"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import multiprocessing
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import toml

from adapted.detect.real_range import real_range_check
from adapted.detect.mvs import mean_var_shift_polyA_detect
from warpdemux.live_balancing.balancer import BarcodeBalancers
from warpdemux.live_balancing.config_parser import MainConfig
from warpdemux.live_balancing.reporting import report_worker
from warpdemux.live_balancing.utils import FailedRead, ReadObject
from warpdemux.live_balancing.worker import (
    balance_worker,
    classification_worker,
    segmentation_worker,
)

from warpdemux import read_until

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


class Logger(object):
    def __init__(self, save_path=None, runid=None):
        self.terminal = sys.stdout

        self.save_path = save_path
        if self.save_path is None:
            self.save_path = os.getcwd()

        if self.save_path is not None:
            if not os.path.isabs(self.save_path):
                self.save_path = os.path.join(REPO_ROOT, self.save_path)

            if not os.path.isdir(self.save_path):
                os.makedirs(self.save_path)

        self.log = open(os.path.join(self.save_path, f"run_{runid}.log"), "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        pass


class ChannelRepeatedUnblockDuration:
    """

    Can't use ReadData number:
    number:
        The minknow assigned number of this read  Read numbers always
        increment throughout the experiment, and are unique per
        channel - however they are not necessarily contiguous.

    Instead, we use time stamps.

    The first unblock is always the default_unblock duration.
    If a second unblock is requested within time_window seconds from the first block,
     the unblock_duration is increased.
    If, in addition, a third unblock is requested within time_window seconds from the second block,
        the unblock_duration is increased again.

    As soon as the time_window is exceeded, the unblock_count is reset to 1.

    For default settings, this automatically resets the unblock_count to 1 after the 3rd (long) unblock.

    """

    def __init__(
        self,
        channels,
        time_window=1.5,
        unblock_duration_1=0.1,
        unblock_duration_2=0.5,
        unblock_duration_3=2.0,
    ):
        self.channels = channels
        self.unblock_time = {ch: 0.1 for ch in channels}
        self.unblock_count = {ch: 0 for ch in channels}
        self.unblock_last = {ch: 0 for ch in channels}

        self.unblock_duration_1 = unblock_duration_1
        self.unblock_duration_2 = unblock_duration_2
        self.unblock_duration_3 = unblock_duration_3
        self.time_window = time_window

    def update(self, channel: int):
        t = time.time()

        diff = t - self.unblock_last[channel]
        if diff < self.time_window:
            self.unblock_count[channel] += 1
            self.unblock_last[channel] = t

            if self.unblock_count[channel] == 2:
                return self.unblock_duration_2
            elif self.unblock_count[channel] == 3:
                return self.unblock_duration_3
            else:
                raise RuntimeError(
                    "unblock_count should be no more than 3; check your logic."
                )

        else:  # more than time_window seconds passed
            # (also happens after 3rd unblock within 3s for default settings)
            self.unblock_last[channel] = t
            self.unblock_count[channel] = 1
            return self.unblock_duration_1


class Session:
    def __init__(
        self,
        client: read_until.ReadUntilClient,
        config_path: str,
    ):
        self.client = client
        self._check_client()

        self.config: MainConfig = self._parse_and_validate_config(config_path)
        self._setup_logger()
        toml.dump(self.config.to_dict(), sys.stdout)

        self.skip_stats = {
            "missed_obs": [],
            "missed_reads": 0,
            "too_long_reads": 0,
            "not_real_read": 0,
        }

        self.CRUD = ChannelRepeatedUnblockDuration(
            channels=self.config.channels,
            time_window=self.config.acquisition.repeated_unblock_time_window,
            unblock_duration_1=self.config.balancing.reject_duration,
            unblock_duration_2=self.config.acquisition.repeated_unblock_duration_2,
            unblock_duration_3=self.config.acquisition.repeated_unblock_duration_3,
        )

        self.balancers = BarcodeBalancers(self.config.balancers)

        self.segmentation_queue = multiprocessing.Queue()
        self.classification_queue = multiprocessing.Queue()
        self.barcode_queue = multiprocessing.Queue()
        self.result_queue = multiprocessing.Queue()

        self.segmentation_executor = ThreadPoolExecutor(
            max_workers=self.config.processing.nproc_segmentation
        )
        self.classification_executor = ThreadPoolExecutor(
            max_workers=self.config.processing.nproc_classification
        )
        self.balancing_executor = ThreadPoolExecutor(1)
        self.reporting_executor = ThreadPoolExecutor(1)

    def _parse_and_validate_config(self, config_path: str):
        return MainConfig(config_path)

    def _setup_logger(self):
        self.logger = Logger(
            self.config.reporting.save_path, self.config.reporting.runid
        )
        sys.stdout = self.logger
        sys.stderr = self.logger

    def _initializer(self):
        for _ in range(self.config.processing.nproc_segmentation):
            self.segmentation_executor.submit(
                segmentation_worker,
                self.segmentation_queue,
                self.classification_queue,
                self.config.model,
            )

        for _ in range(self.config.processing.nproc_classification):
            self.classification_executor.submit(
                classification_worker,
                self.classification_queue,
                self.barcode_queue,
                self.config.model,
            )

        self.balancing_executor.submit(
            balance_worker,
            self.barcode_queue,
            self.result_queue,
            self.config.balancing,
            self.client,
            self.balancers,
        )

        self.reporting_executor.submit(
            report_worker,
            self.result_queue,
            self.config.reporting,
            self.skip_stats,
            self.balancers,
        )

    def _check_client(self):
        if not self.client.signal_dtype == np.float32:
            raise ValueError("Client calibrated_signal should be True")
        if self.client.one_chunk:
            raise ValueError("Client should be set to one_chunk=False")
        if not self.client.cache_type.__name__ in [
            "AccumulatingCache",
            "PreallocAccumulatingCache",
        ]:
            print(
                self.client.cache_type.__name__,
            )
            raise ValueError("Client cache_type should be AccumulatingCache")

    def start(self):
        self._initializer()

        # Having a running client will cause very weird shenenigans in the signal buffer.
        # Trust me, don't do it.
        if self.client.is_running:
            self.client.reset()

        self.client.run()

    def stop(self):
        self.client.reset()

        for _ in range(self.config.processing.nproc_segmentation):
            self.segmentation_queue.put(None)
        self.segmentation_executor.shutdown()

        for _ in range(self.config.processing.nproc_classification):
            self.classification_queue.put(None)
        self.classification_executor.shutdown()

        self.barcode_queue.put(None)
        self.balancing_executor.shutdown()
        self.balancers.stop()

        self.result_queue.put(None)
        self.reporting_executor.shutdown()

        time.sleep(0.5)

        assert self.segmentation_queue.empty()
        self.segmentation_queue.close()

        assert self.classification_queue.empty()
        self.classification_queue.close()

        assert self.barcode_queue.empty()
        self.barcode_queue.close()

        assert self.result_queue.empty()
        self.result_queue.close()

    def run(self, batch_size=100):
        """Run the analysis on the data from the client"""
        self.start()

        while self.client.is_running:
            if self.client.queue_length > 0:
                for channel, read in self.client.get_read_chunks(
                    batch_size=batch_size,
                    last=False,  # False -> FIFO, True -> LIFO,
                    min_chunk_length=self.config.acquisition.min_chunk_size,  # we don't expect poly A before this
                ):
                    start_time = time.time()
                    calibrated_signal = np.frombuffer(
                        read.raw_data, self.client.signal_dtype
                    )

                    missed_obs = (
                        read.chunk_start_sample - read.start_sample
                    )  # if negative then start_sample is later than chunk_start, so within first chunk
                    self.skip_stats["missed_obs"].append(missed_obs)

                    if missed_obs > self.config.acquisition.max_missed_start_offset:
                        self.client.stop_receiving_read(channel, read.number)
                        self.skip_stats["missed_reads"] += 1

                        self.result_queue.put(
                            FailedRead(
                                reason="missed_obs",
                                decision="retain",
                                channel=channel,
                                read_start_sample=read.start_sample,
                                chunk_start_sample=read.chunk_start_sample,
                                read_number=read.number,
                                read_id=read.id,
                                chunk_length=calibrated_signal.size,
                                polya_start=0,
                                data_arr=np.array([]),
                                start_time=start_time,
                                time_per_step=[time.time() - start_time],
                            ).to_result(self.config.model.num_bcs)
                        )
                        continue

                    # trimmning start of signal to read start
                    # NOTE: polya_start is thus reported relative to the start of the read!
                    if missed_obs < 0:
                        calibrated_signal = calibrated_signal[-missed_obs:]

                    chunk_size = calibrated_signal.size

                    if chunk_size > self.config.acquisition.max_chunk_size:
                        self.skip_stats["too_long_reads"] += 1
                        self.client.stop_receiving_read(channel, read.number)
                        # if read is too long dont decide on it
                        self.result_queue.put(
                            FailedRead(
                                reason="too_long_read",
                                decision="retain",
                                channel=channel,
                                read_start_sample=read.start_sample,
                                chunk_start_sample=read.chunk_start_sample,
                                read_number=read.number,
                                read_id=read.id,
                                chunk_length=chunk_size,
                                polya_start=0,
                                data_arr=np.array([]),
                                start_time=start_time,
                                time_per_step=[time.time() - start_time],
                            ).to_result(self.config.model.num_bcs)
                        )
                        continue

                    polya_start = mean_var_shift_polyA_detect(
                        calibrated_signal,
                        params=self.config.model.sdc.streaming,
                    )

                    if (
                        polya_start == 0
                        and chunk_size <= self.config.acquisition.max_chunk_size
                    ):
                        continue

                    # print(f"Channel: {channel}, Read id: {read.id}, Signal length:{calibrated_signal.size}, Chunk size: {read.chunk_length}, Poly A position: {polya_start}, missed_obs: {missed_obs}")
                    # np.save(f"debug/{read.id}", calibrated_signal)
                    # np.save(f"debug/copied_{read.id}", calibrated_signal)

                    # arriving here means a poly A tail was detected
                    self.client.stop_receiving_read(channel, read.number)

                    # check that basic signal characteristics match
                    if not real_range_check(
                        calibrated_signal[:polya_start],
                        self.config.model.sdc.real_range,
                    ):
                        self.skip_stats["not_real_read"] += 1
                        self.result_queue.put(
                            FailedRead(
                                reason="not_real_read",
                                decision="reject",
                                channel=channel,
                                read_start_sample=read.start_sample,
                                chunk_start_sample=read.chunk_start_sample,
                                read_number=read.number,
                                read_id=read.id,
                                chunk_length=chunk_size,
                                polya_start=0,
                                data_arr=np.array([]),
                                start_time=start_time,
                                time_per_step=[time.time() - start_time],
                            ).to_result(self.config.model.num_bcs)
                        )
                        self.client.unblock_read(
                            channel, read.number, duration=self.CRUD.update(channel)
                        )

                        continue

                    clip = min(
                        polya_start + self.config.model.spc.sig_extract.padding,
                        chunk_size,
                    )

                    self.segmentation_queue.put(
                        ReadObject(
                            channel=channel,
                            read_start_sample=read.start_sample,
                            chunk_start_sample=read.chunk_start_sample,
                            read_number=read.number,
                            read_id=read.id,
                            chunk_length=chunk_size,
                            polya_start=polya_start,
                            data_arr=calibrated_signal[:clip],
                            start_time=start_time,
                            time_per_step=[time.time() - start_time],
                        )
                    )
                time.sleep(0.005)

        self.stop()

    def reject_all(self, batch_size=100, wait_between_requests=0.001, duration=None):
        """Reject all reads immediately (non-batched) when they appear in queue"""
        self.start()

        start = time.time()
        report_duration = 0
        while self.client.is_running:
            if self.client.queue_length > 0:
                for channel, read in self.client.get_read_chunks(
                    batch_size=batch_size, last=False  # False -> FIFO, True -> LIFO
                ):
                    self.client.unblock_read(channel, read.number)

            if time.time() - report_duration >= 10:
                report_duration = time.time()
                # print(missed_offsets)
                print("Total read events missed:", self.client.missed_reads)
            time.sleep(wait_between_requests)
            if duration is not None:
                if time.time() - start > duration:
                    break

        self.stop()
        return self.client.missed_reads
