"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import multiprocessing
import os
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import pod5

from warpdemux.live_balancing.config_parser import BalancerConfig
from warpdemux.live_balancing.utils import BALANCE_TYPE_WATCHER


@dataclass
class ReadInfo:
    num_minknow_events: int
    end_reason: str
    forced: bool
    channel: int


class Pod5ReadSets:
    """
    stats on reads that were accepted, independ of whether they showed up in pod5 or not:
     - seen_read_ids_per_bc
     - seen_read_count_per_bc
     - seen_read_count (total count accross bcs)

    stats on reads that were accepted and showed up in pod5:
     - balancer.stats_per_bc
     - matched_read_count_per_bc
     - matched_kbases_per_bc
     - matched_read_count (total count accross bcs)
     - matched_kbases (total count accross bcs)
     - matched_estimated_read_length (total count accross bcs)
     - matched_estimated_read_length_per_bc

    stats on reads that were not accepted/seen but showed up in pod5:
     - pod5_unmatched_read_count
     - pod5_unmatched_read_ids
     - pod5_unmatched_kbases

    stats on reads that were accepted but did not show up in pod5:
     - unmatched_seen_read_ids_per_bc
     - unmatched_seen_read_ids
     - unmatched_seen_read_count_per_bc
     - unmatched_seen_read_count
    """

    def __init__(self, num_bcs: int, name: str):
        """sets_seen is of length num_bcs, each element is a set of read_ids"""

        self.num_bcs = num_bcs
        self._sets_seen = []
        for _ in range(self.num_bcs):
            self._sets_seen.append(set())

        self._pod5_reads = set()
        self._pod5_reads_to_match = (
            set()
        )  # reads that showed up in pod5 but were not accepted in run
        self._reads_matched_per_bc = [set() for _ in range(self.num_bcs)]

        self._pod5_kbases = 0
        self._matched_kbases = 0
        self._matched_kbases_per_bc = [0.0 for _ in range(self.num_bcs)]

        self.name = name

    def check_remaining_to_match(self):
        matched = set()

        for read_id, num_kbases in self._pod5_reads_to_match:
            for i in range(self.num_bcs):
                if read_id in self._sets_seen[i]:
                    self._reads_matched_per_bc[i].add((read_id, num_kbases))
                    self._matched_kbases += num_kbases
                    self._matched_kbases_per_bc[i] += num_kbases
                    matched.add((read_id, num_kbases))
                    break

        self._pod5_reads_to_match -= matched
        return len(matched)

    def add_seen_in_run(self, index: int, read_id: str):
        self._sets_seen[index].add(read_id)

    def add_seen_in_pod5(self, read_id: str, num_kbases: float, rejected: bool = False):
        if rejected:
            return 0  # ignore reads that we rejected

        if num_kbases < 20 / 1000:
            return 0  # ignore reads shorter than 20 bases

        self._pod5_reads.add((read_id, num_kbases))
        self._pod5_kbases += num_kbases

        for i in range(self.num_bcs):
            if read_id in self._sets_seen[i]:
                self._reads_matched_per_bc[i].add((read_id, num_kbases))
                self._matched_kbases += num_kbases
                self._matched_kbases_per_bc[i] += num_kbases
                return 1

        # not matched
        self._pod5_reads_to_match.add((read_id, num_kbases))
        return 0

    def add_seen_in_pod5_read_dict(self, read_dict: Dict[str, ReadInfo]):
        n_matched = 0
        # TODO: clean up unused fields and make into named tuple for readability

        for read_id, read_info in read_dict.items():
            n_matched += self.add_seen_in_pod5(
                read_id,
                max(
                    [
                        0,
                        (read_info.num_minknow_events - 100)  # subtract adapter events
                        / 1000,  # convert to kbases
                    ]
                ),
                read_info.forced,
            )

        return n_matched

    @property
    def seen_read_ids_per_bc(self):
        return self._sets_seen

    @property
    def seen_read_count_per_bc(self):
        return np.array([len(s) for s in self._sets_seen])

    @property
    def seen_read_count(self):
        return sum([len(s) for s in self._sets_seen])

    @property
    def pod5_read_count(self):
        return len(self._pod5_reads)

    @property
    def pod5_kbases(self):
        return self._pod5_kbases

    @property
    def pod5_estimated_read_length(self):
        return self.pod5_kbases / self.pod5_read_count

    @property
    def matched_read_count_per_bc(self):
        return np.array([len(s) for s in self._reads_matched_per_bc])

    @property
    def matched_read_ids_per_bc(self):
        return [
            set([read_id for read_id, _ in self._reads_matched_per_bc[i]])
            for i in range(self.num_bcs)
        ]

    @property
    def matched_read_count(self):
        return sum([len(s) for s in self._reads_matched_per_bc])

    @property
    def matched_kbases(self):
        return self._matched_kbases

    @property
    def matched_kbases_per_bc(self):
        return np.array(self._matched_kbases_per_bc, dtype=float)

    @property
    def matched_estimated_read_length(self):
        return (
            self.matched_kbases / self.matched_read_count
            if self.matched_read_count
            else 0
        )

    @property
    def matched_estimated_read_length_per_bc(self):
        return np.divide(
            self.matched_kbases_per_bc,
            self.matched_read_count_per_bc.astype(float),
            out=np.zeros_like(self.matched_kbases_per_bc),
            where=self.matched_read_count_per_bc != 0,
        )

    @property
    def matches_df(self):
        try:
            res = []
            for i in range(self.num_bcs):
                df = pd.DataFrame(
                    self._reads_matched_per_bc[i], columns=["read_id", "kbases"]
                )
                df["barcode"] = i
                res.append(df)
            return pd.concat(res)
        except:
            print("Exception in matches_df")
            traceback.print_exc()
            return

    @property
    def unmatched_seen_read_ids_per_bc(self):
        return [
            self._sets_seen[i] - self.matched_read_ids_per_bc[i]
            for i in range(self.num_bcs)
        ]

    @property
    def unmatched_seen_read_ids(self):
        return set().union(*self.unmatched_seen_read_ids_per_bc)

    @property
    def unmatched_seen_read_count_per_bc(self):
        return np.array([len(s) for s in self.unmatched_seen_read_ids_per_bc])

    @property
    def unmatched_seen_read_count(self):
        return sum([len(s) for s in self.unmatched_seen_read_ids_per_bc])

    @property
    def unmatched_seen_df(self):
        try:
            res = []
            for i in range(self.num_bcs):
                df = pd.DataFrame(
                    self.unmatched_seen_read_ids_per_bc[i],
                    columns=["read_id"],
                )
                df["barcode"] = i
                res.append(df)
            return pd.concat(res)
        except:
            print("Exception in matches_df")
            traceback.print_exc()
            return

    @property
    def unmatched_pod5_read_count(self):
        return len(self._pod5_reads_to_match)

    @property
    def unmatched_pod5_kbases(self):
        return sum([kbases for _, kbases in self._pod5_reads_to_match])

    @property
    def unmatched_pod5_read_ids(self):
        return [read_id for read_id, _ in self._pod5_reads_to_match]


class BarcodeBalancer:
    pod5_watcher: ThreadPoolExecutor = None
    pod5_watcher_queue: multiprocessing.Queue = None
    pod5_read_sets: Pod5ReadSets = None

    def __init__(
        self,
        num_bcs: int,
        balance_threshold: float = 0.01,
        min_stat: int = 100,
        name: str = "adapter_count",
        balance_type: str = "adapter_count",
        pod5_watch_dir: Optional[str] = None,
        channels: np.ndarray = None,
        pod5_check_interval: int = None,
        wait_to_see: int = None,
        watch_for_missing: bool = None,
        blacklist: np.ndarray = None,
        ignorelist: np.ndarray = None,
        max_stats: np.ndarray = None,
        reject_duration: float = None,
        # TODO: replace init args with BalancerConfig object
    ):
        self.start_time = time.time()

        self.config = {}
        self.name = name
        self.num_bcs = num_bcs
        self.balance_threshold = balance_threshold
        self.min_stat = min_stat
        self.reject_duration = reject_duration

        self._counts_per_bc = np.zeros(num_bcs, dtype=int)  # used for non-pod5 modes

        self.balance_type = balance_type
        self.is_watcher = balance_type in BALANCE_TYPE_WATCHER

        if self.is_watcher:
            if pod5_check_interval is None:
                raise ValueError(
                    f"pod5_check_interval is required for balancer {self.name} of type {self.balance_type}"
                )
            if not os.path.isdir(pod5_watch_dir):
                raise ValueError(
                    f"pod5_watch_dir of balancer {self.name} is not a directory"
                )

        self.channels = channels

        self.blacklist = blacklist
        self.ignorelist = ignorelist
        self.max_stats = max_stats

        self.watch_for_missing = watch_for_missing
        self.wait_to_see = wait_to_see

        if self.balance_type == "block_safety":
            raise NotImplementedError("Block safety not implemented yet.")

        elif self.balance_type == "none":
            self.decide = self._decide_none
            self.update = self._update_none

        elif self.balance_type == "adapter_count":
            self.decide = self._decide_on_balance
            self.update = self._update_count

        elif self.balance_type == "read_count":
            self.decide = self._decide_on_balance
            self.update = self._update_watcher
            self._initialize_pod5_watcher(pod5_watch_dir, pod5_check_interval)

        elif self.balance_type == "base_normalization":
            self.decide = self._decide_on_balance
            self.update = self._update_watcher
            self._initialize_pod5_watcher(pod5_watch_dir, pod5_check_interval)

        elif self.balance_type == "reject_all":
            self.decide = self._reject_all
            self.update = self._update_count
        else:
            raise ValueError(f"Unknown type of balance {balance_type}.")

    @staticmethod
    def pod5_watch_worker(
        input_queue, pod5_dir, pod5_read_sets: Pod5ReadSets, pod5_check_interval=1
    ):
        files_seen_set = set()

        def get_pod5_paths():
            """fpath should be abs path to directory were pod5 file is"""
            return [
                f
                for f in os.listdir(pod5_dir)
                if (os.path.splitext(f)[1] == ".pod5" and f not in files_seen_set)
            ]

        def get_reads_from_pod5(pod5_path) -> Union[None, Dict[str, ReadInfo]]:
            read_dict = {}

            try:
                with pod5.Reader(pod5_path) as reader:
                    for read in reader.reads():
                        read_dict[str(read.read_id)] = ReadInfo(
                            num_minknow_events=int(read.num_minknow_events),
                            end_reason=str(read.end_reason.reason),
                            forced=read.end_reason.forced,
                            channel=read.pore.channel,
                        )

                return read_dict
            except IOError:  # file is still being written
                return None

        def run_watcher():
            try:
                for pod5_path in get_pod5_paths():
                    read_dict = get_reads_from_pod5(os.path.join(pod5_dir, pod5_path))

                    if read_dict is None:
                        continue

                    n_matched = pod5_read_sets.add_seen_in_pod5_read_dict(read_dict)
                    files_seen_set.add(pod5_path)

                    print(
                        f"[pod5_watch_worker] Found {n_matched} read_ids from"
                        f" {len(read_dict)} total in file"
                        f" {os.path.basename(pod5_path)}."
                    )

                n_remaining = (
                    pod5_read_sets.unmatched_pod5_read_count
                )  # this should only happen during dummy testing
                n_matched = pod5_read_sets.check_remaining_to_match()
                if n_matched:
                    print(
                        f"[pod5_watch_worker] Found {n_matched} read_ids from"
                        f" {n_remaining} total in remaining reads."
                    )

            except:
                print("Exception in pod5 watch worker")
                traceback.print_exc()
                return

        while True:
            try:
                item = input_queue.get(block=False)
                if item is None:
                    break
            except Exception:
                pass

            try:
                run_watcher()
                time.sleep(
                    pod5_check_interval
                )  # wait before retrying to give files time to be written
            except:
                print("Exception in pod5 watch worker")
                traceback.print_exc()
                return

        # after stop signal, run watcher once more to catch any remaining reads
        run_watcher()

    def _initialize_pod5_watcher(self, pod5_dir, pod5_check_interval):
        self.pod5_watcher = ThreadPoolExecutor(1)
        self.pod5_watcher_queue = multiprocessing.Queue()
        self.pod5_read_sets = Pod5ReadSets(self.num_bcs, self.name)

        self.pod5_watcher.submit(
            self.pod5_watch_worker,
            self.pod5_watcher_queue,
            pod5_dir,
            self.pod5_read_sets,
            pod5_check_interval,
        )

    @property
    def stats_per_bc(self):
        if self.balance_type == "read_count":
            return self.pod5_read_sets.matched_read_count_per_bc
        elif self.balance_type == "base_normalization":
            return self.pod5_read_sets.matched_kbases_per_bc
        elif self.balance_type in ["adapter_count", "none", "reject_all"]:
            return self._counts_per_bc
        else:
            raise NotImplementedError

    def _reject_all(self, index: int, read_id: str, channel: int) -> Tuple[bool, str]:
        if not (np.isin(channel, self.channels, assume_unique=True)):
            return (
                True,
                "untracked_channel",
            )  # Do not balance on barcodes from other channels

        if self.ignorelist[index]:
            return True, "ignored"

        return False, "reject_all"

    def _decide_none(self, index: int, read_id: str, channel: int) -> Tuple[bool, str]:
        if self.blacklist[index]:
            return False, "blacklisted"

        return True, "no_balancing"

    def _update_none(self, index: int, read_id: str, channel: int) -> None:
        pass

    def _decide_on_balance(
        self, index: int, read_id: str, channel: int
    ) -> Tuple[bool, str]:
        if not (np.isin(channel, self.channels, assume_unique=True)):
            return (
                True,
                "untracked_channel",
            )  # Do not balance on barcodes from other channels

        if self.blacklist[index]:
            return False, "blacklisted"

        if self.ignorelist[index]:
            return True, "ignored"

        if self.stats_per_bc[index] < self.min_stat:
            return True, "min_stat"

        if self.stats_per_bc[index] > self.max_stats[index]:
            self._put_on_blacklist(index)
            return False, "max_stat"

        indexer = np.logical_or(self.blacklist, self.ignorelist)
        valid_barcode_sum = self.stats_per_bc.sum() - self.stats_per_bc[indexer].sum()
        valid_num_bcs = self.num_bcs - indexer.sum()

        if valid_barcode_sum <= 0:
            return True, "valid_sum_seen_zero"

        offset = self.stats_per_bc[index] - valid_barcode_sum / valid_num_bcs
        threshold = self.balance_threshold * valid_barcode_sum / valid_num_bcs

        if offset > threshold:
            return False, "unbalanced"

        return True, "balanced"

    def _update_watcher(self, index: int, read_id: str, channel: int):
        self.pod5_read_sets.add_seen_in_run(index, read_id)

    def _update_count(self, index: int, read_id: str, channel: int):
        self.stats_per_bc[index] += 1

    def _put_on_blacklist(self, index: int):
        self.blacklist[index] = True

        if self.ignorelist[index]:
            self.ignorelist[index] = False

    def _put_on_ignorelist(self, index: int):
        self.ignorelist[index] = True

        if self.blacklist[index]:
            self.blacklist[index] = False

    def _handle_missing_barcodes(self):
        if self.watch_for_missing:
            for idx in np.argwhere(self.stats_per_bc < self.min_stat).ravel():
                self._put_on_ignorelist(idx)

    def balance(self, assigned_bc: int, read_id: str, channel: int):
        """
        Decide whether to accept the barcode based on the balance criterion.

        Returns:
        bool: True if the barcode is accepted, False otherwise.
        """
        if self.watch_for_missing and time.time() - self.start_time > self.wait_to_see:
            self._handle_missing_barcodes()
            self.watch_for_missing = False

        decision, reason = self.decide(assigned_bc, read_id, channel)
        if decision:
            self.update(
                assigned_bc, read_id, channel
            )  # watchers are only updated for accepted reads
            return True, reason
        else:
            return False, reason

    def stop(self):
        if self.is_watcher:
            self.pod5_watcher_queue.put(None)
            self.pod5_watcher.shutdown()
            self.pod5_watcher_queue.close()


class BarcodeBalancers:
    def __init__(self, config_list: List[BalancerConfig]):
        self.config_list = config_list

        self._initialize_balancers()
        self._initialize_channel_dict()

    def _initialize_balancers(self):
        self.balancers_dict = {}

        for config in self.config_list:
            self.balancers_dict[config.name] = BarcodeBalancer(
                name=config.name,
                num_bcs=config.num_bcs,
                balance_threshold=config.balance_threshold,
                min_stat=config.min_stat,
                balance_type=config.balance_type,
                channels=config.channels,
                pod5_watch_dir=config.pod5_watch_dir,
                pod5_check_interval=config.pod5_check_interval,
                watch_for_missing=config.watch_for_missing,
                wait_to_see=config.wait_to_see,
                blacklist=config.blacklist,
                ignorelist=config.ignorelist,
                max_stats=config.max_stats,
                reject_duration=config.reject_duration,
            )

    def _initialize_channel_dict(self):
        self.channel_dict = {}
        for balancer in self.balancers:
            for channel in balancer.channels:
                self.channel_dict[channel] = balancer.name

    @property
    def balancers(self) -> List[BarcodeBalancer]:
        return list(self.balancers_dict.values())

    def _get_balancer_with_channel(self, channel: int) -> BarcodeBalancer:
        name = self.channel_dict[channel]
        return self._get_balancer_with_name(name)

    def _get_balancer_with_name(self, name: str) -> Optional[BarcodeBalancer]:
        return self.balancers_dict.get(name, None)

    def get_balancer(
        self, name_or_channel: Union[str, int]
    ) -> Optional[BarcodeBalancer]:
        if isinstance(name_or_channel, int):
            return self._get_balancer_with_channel(name_or_channel)
        else:
            return self._get_balancer_with_name(name_or_channel)

    @property
    def names(self) -> List[str]:
        return list(self.balancers_dict.keys())

    def balance_type_for_channel(self, channel: int) -> str:
        return self.get_balancer(channel).balance_type

    def balance(self, assigned_bc: int, read_id: str, channel: int) -> Tuple[bool, str]:
        """
        Decide whether to accept the barcode based on the balance criterion.

        Returns:
        bool: True if the barcode is accepted, False otherwise.
        """

        balancer = self.get_balancer(channel)
        if balancer is None:
            raise ValueError(f"No balancer found for channel {channel}")

        return balancer.balance(assigned_bc, read_id, channel)

    def stop(self) -> None:
        for balancer in self.balancers:
            balancer.stop()
