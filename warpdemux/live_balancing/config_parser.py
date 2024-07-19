"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import os
import time
from typing import List

import numpy as np
import toml

from warpdemux.config.utils import (
    get_model_num_bcs,
    get_model_spc_config,
    get_model_spc_live_config,
)
from warpdemux.live_balancing._defaults import *
from warpdemux.live_balancing.utils import BALANCE_TYPE_WATCHER

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

CHANNEL_MIN = 1


class ConfigParserBase:
    section = None

    def __init__(self, config: dict):
        pass

    def _check_for_unknown_keys(self, config: dict, ignore_keys: List[str] = None):
        ignore_keys = ignore_keys or []
        try:
            for key in config[self.section].keys():
                if key in ignore_keys:
                    pass
                elif key not in self.__dict__:
                    raise ValueError(f"Unknown key {key} in config [{self.section}].")
        except KeyError:
            for key in config.keys():
                if key in ignore_keys:
                    pass
                elif key not in self.__dict__:
                    raise ValueError(f"Unknown key {key} in config.")

    def to_dict(self):
        d = self.__dict__.copy()
        d.pop("section", None)
        return d


class ModelConfig(ConfigParserBase):
    """
    model_name: str
        Name of the model to use for classification.
    num_bcs: int
        Number of barcodes to report.
    spc: SignalProcessingConfig
        signal processing config
    """

    section = "model"

    def __init__(self, config: dict):
        self.model_name = config[self.section].get("model_name", DEFAULT_MODEL_NAME)

        self.spc = get_model_spc_config(self.model_name)
        self.sdc = get_model_spc_live_config(self.model_name)
        self.num_bcs = get_model_num_bcs(self.model_name)

        self._check_for_unknown_keys(
            config,
        )

        print(self.sdc)


class ProcessingConfig(ConfigParserBase):
    """
    nproc_segmentation: int
        Number of processes to use for segmentation.
    nproc_classification: int
        Number of processes to use for classification.
    """

    section = "processing"

    def __init__(self, config: dict):
        if self.section not in config.keys():
            config[self.section] = {}

        self.nproc_segmentation = config[self.section].get(
            "nproc_segmentation", DEFAULT_NPROC_SEGMENTATION
        )
        self.nproc_classification = config[self.section].get(
            "nproc_classification", DEFAULT_NPROC_CLASSIFICATION
        )

        self._check_for_unknown_keys(config)


class AcquisitionConfig(ConfigParserBase):
    """
    max_missed_start_offset: int
        Maximum number of samples to miss the start offset by.
    max_chunk_size: int
        Maximum number of samples to extract from each read.
    min_chunk_size: int
        Minimum length of accumulated chunk to be retrieved.
    min_adapter_length: int
        Minimum length of adapter to be detected.
    repeated_unblock_time_window: float
        Time window in seconds to consider for repeated unblock detection.
    repeated_unblock_duration_2: float
        Rejecetion durations in seconds for first repeated unblock detection.
    repeated_unblock_duration_3: float
        Rejecetion durations in seconds for second repeated unblock detection.
    """

    section = "acquisition"

    def __init__(self, config: dict):
        if self.section not in config.keys():
            config[self.section] = {}

        self.max_missed_start_offset = config[self.section].get(
            "max_missed_start_offset", DEFAULT_MAX_MISSED_START_OFFSET
        )
        self.max_chunk_size = config[self.section].get(
            "max_chunk_size", DEFAULT_MAX_CHUNK_SIZE
        )
        self.min_chunk_size = config[self.section].get(
            "min_chunk_size", DEFAULT_MIN_CHUNK_SIZE
        )

        self.min_adapter_length = config[self.section].get(
            "min_adapter_length", self.min_chunk_size
        )

        if self.min_chunk_size > self.max_chunk_size:
            raise ValueError(
                f"min_chunk_size {self.min_chunk_size} can't be larger than"
                f" max_chunk_size {self.max_chunk_size}. Please check your config."
            )

        # TODO: refactor these params to their own section?
        self.repeated_unblock_time_window = config[self.section].get(
            "repeated_unblock_time_window", DEAFULT_REPEATED_UNBLOCK_TIME_WINDOW
        )
        self.repeated_unblock_duration_2 = config[self.section].get(
            "repeated_unblock_duration_2", DEFAULT_REPEATED_UNBLOCK_DURATION_2
        )
        self.repeated_unblock_duration_3 = config[self.section].get(
            "repeated_unblock_duration_3", DEFAULT_REPEATED_UNBLOCK_DURATION_3
        )
        self._check_for_unknown_keys(config)


class BalancingConfig(ConfigParserBase):
    """
    pred_conf_threshold: float
        Confidence threshold for the classifier. If the confidence is below this threshold, the read is considered unclassified.
    reject_duration: float
        Duration in seconds to reject a channel after a read is rejected.
    max_signal_after_polya: int
        Number of signal allowed after detection of poly A start. Used to suppresss rejection commands when read is too far translocated.
    """

    section = "balancing"

    def __init__(self, config: dict):
        if self.section not in config.keys():
            config[self.section] = {}
        self.pred_conf_threshold = config[self.section].get(
            "pred_conf_threshold", DEFAULT_PRED_CONF_THRESHOLD
        )
        self.reject_duration = config["balancing"].get(
            "reject_duration", DEFAULT_REJECT_DURATION
        )
        self.max_signal_after_polya = config["balancing"].get(
            "max_signal_after_polya", DEFAULT_MAX_SIGNAL_AFTER_POLYA
        )
        self._check_for_unknown_keys(config)


class ReportingConfig(ConfigParserBase):
    """
    save_path: str
        Path to save the report and results to.
    save_every_sec: int
        Save the report every save_every_sec seconds.
    runid: str
        Unique identifier for the run based on the current time.
    """

    section = "reporting"

    def __init__(self, config: dict):
        if self.section not in config.keys():
            config[self.section] = {}
        self.save_path = config[self.section].get("save_path", None)
        self.runid = time.strftime("%Y-%m-%d_%H-%M-%S")
        self.save_every_sec = config[self.section].get(
            "save_every_sec", DEFAULT_SAVE_EVERY_SEC
        )

        self._validate_save_path()

        self._check_for_unknown_keys(config)

    def _validate_save_path(self):
        if self.save_path is None:
            self.save_path = os.path.join(os.getcwd(), "results")
        if not os.path.isabs(self.save_path):
            self.save_path = os.path.join(REPO_ROOT, self.save_path)
        if not os.path.isdir(self.save_path):
            os.makedirs(self.save_path)


class BalancerConfig(ConfigParserBase):
    """
    name: str
        Name of the balancer. Defaults to `balance_type`.
    balance_threshold: float
        Threshold for the balance statistic. If the balance statistic is below this threshold, the read is rejected.
    min_stat: float
        Minimum number of reads or kbases to have seen for a barcode before balancing kicks in for that barcode.
    balance_type: str
        Type of balancer to use. One of ['none', 'reject_all', 'adapter_count', 'read_count', 'base_normalization']
    pred_conf_threshold: float
        Confidence threshold for the classifier. If the confidence is below this threshold, the read is considered unclassified.
    channel_frac: float
        Fraction of channels to use for this balancer.
    channel_num: int
        Number of channels to use for this balancer.
    pod5_watch_dir: str
        Path to the directory to watch for pod5 files. Required for balance_type 'read_count' and 'base_normalization'.
    pod5_check_interval: float
        Interval in seconds to check for new pod5 files. Required for balance_type 'read_count' and 'base_normalization'.
    channels: list
        List of channels to use for this balancer.
    watch_for_missing: bool
        Boolean flag to watch for missing barcodes. Missing barcodes are set on ignorelist.
    wait_to_see: int
        Time in seconds to wait before deciding if a barcode is missing.
    blacklist: np.ndarray[bool]
        Boolean array of size num_bcs. If the entry with index `i` is `True`, the barcode with index `i` is blacklisted.
    ignorelist: np.ndarray[bool]
        Boolean array of size num_bcs. If the entry with index `i` is `True`, the barcode with index `i` is ignored.
    max_stats: np.ndarray[float]
        Float array of size num_bcs. If the number of reads or kbases for a barcode exceeds this value, the barcode is blacklisted.
    reject_duration: float
        Balancer specific duration in seconds to reverse the channel charge when rejecting a read.
    """

    section = "balancers_item"

    def __init__(
        self,
        config: dict,
        num_bcs: int,
        n_all_channels: int,
        available_channels: list,
    ):
        if self.section not in config.keys():
            config[self.section] = {}

        self.num_bcs = num_bcs
        self.balance_threshold = config.get(
            "balance_threshold", DEFAULT_BALANCE_THRESHOLD
        )
        self.min_stat = config.get("min_stat", DEFAULT_MIN_STAT)
        self.balance_type = config.get("balance_type", DEFAULT_BALANCE_TYPE)
        self.name = config.get("name", self.balance_type)
        self.pred_conf_threshold = config.get(
            "pred_conf_threshold", DEFAULT_PRED_CONF_THRESHOLD
        )
        self.channel_num = config.get("channel_num", None)
        if self.channel_num is not None:
            if config.get("channel_frac", None) is not None:
                raise ValueError(
                    "Only one of channel_frac and channel_num can be specified."
                )
            self.channel_frac = None
        else:
            self.channel_frac = config.get("channel_frac", DEFAULT_CHANNEL_FRAC)

        self.reject_duration = config.get("reject_duration", None)
        self.pod5_watch_dir = config.get("pod5_watch_dir", None)
        self.pod5_check_interval = None

        self.watch_for_missing = config.get(
            "watch_for_missing", DEFAULT_WATCH_FOR_MISSING
        )
        self.wait_to_see = config.get("wait_to_see", DEFAULT_WAIT_TO_SEE)

        self.blacklist = np.zeros(self.num_bcs, dtype=bool)
        for i in range(self.num_bcs):
            if config.get(f"blacklist_barcode{i:02d}", False):
                self.blacklist[i] = True
                config.pop(f"blacklist_barcode{i:02d}", None)

        self.ignorelist = np.zeros(self.num_bcs, dtype=bool)
        for i in range(self.num_bcs):
            if not config.get(f"watch_barcode{i:02d}", True):
                self.ignorelist[i] = True
                config.pop(f"watch_barcode{i:02d}", None)

        self.max_stats = np.full(self.num_bcs, np.inf)
        for i in range(self.num_bcs):
            if config.get(f"max_barcode{i:02d}", None) is not None:
                self.max_stats[i] = config.get(f"max_barcode{i:02d}")
                config.pop(f"max_barcode{i:02d}", None)

        # sanity check, barcode can't be both blacklisted and ignored
        for i in range(self.num_bcs):
            if self.blacklist[i] and self.ignorelist[i]:
                raise ValueError(
                    f"Barcode {i} can't be both blacklisted and ignored. Please check"
                    " your config."
                )

        if self.balance_type in BALANCE_TYPE_WATCHER:
            if self.pod5_watch_dir in [None, ""]:
                raise ValueError(
                    f"pod5_watch_dir is required for mode {self.balance_type}"
                )
            if not os.path.isabs(self.pod5_watch_dir):
                self.pod5_watch_dir = os.path.join(REPO_ROOT, self.pod5_watch_dir)

            self.pod5_check_interval = config.get(
                "pod5_check_interval", DEFAULT_POD5_CHECK_INTERVAL
            )

        n_channels = (
            self.channel_num
            if self.channel_num is not None
            else int(self.channel_frac * n_all_channels)
        )
        try:
            self.channels = np.array(
                sorted([available_channels.pop() for _ in range(n_channels)])
            )
        except StopIteration:
            raise ValueError(
                f"channel_frac {self.channel_frac} is too high, only"
                f" {len(available_channels)} channels are available. This can happen"
                " when you did not specify channel_frac for each individual balancer."
            )

        self._check_for_unknown_keys(config)

    def to_dict(self):
        d = super(self.__class__, self).to_dict()
        d["channels"] = self.channels.tolist()
        d["blacklist"] = np.argwhere(self.blacklist).ravel().tolist()
        d["ignorelist"] = np.argwhere(self.ignorelist).ravel().tolist()
        return d


class FlowcellConfig(ConfigParserBase):
    """
    flowcell_type: str
        One of ['flongle', 'minion', 'promethion']
    channel_num: int
        Number of channels on the flowcell.
    min_channel: int
        Minimum channel number to use.
    max_channel: int
        Maximum channel number to use.
    """

    section = "flowcell"

    channel_num_dict = {
        "flongle": 126,
        "minion": 512,
        "promethion": 2675,
    }

    def __init__(self, config: dict):
        if self.section not in config.keys():
            raise ValueError(f"Flowcell section missing in config.")

        self.flowcell_type = config[self.section].get("flowcell_type", None)
        if self.flowcell_type is None:
            raise ValueError(f"Flowcell type missing in config.")
        if self.flowcell_type not in self.channel_num_dict.keys():
            raise ValueError(
                f"Unknown flowcell type {self.flowcell_type}. Supported types are {self.channel_num_dict.keys()}."
            )

        self.channel_num = self.channel_num_dict[self.flowcell_type]
        self.min_channel = config[self.section].get("min_channel", DEFAULT_MIN_CHANNEL)
        self.max_channel = config[self.section].get("max_channel", self.channel_num)
        if self.min_channel < CHANNEL_MIN:
            raise ValueError(
                f"min_channel {self.min_channel} can't be smaller than"
                f" {CHANNEL_MIN} (flowcell {self.flowcell_type}). Please check your config."
            )
        if self.max_channel > self.channel_num:
            raise ValueError(
                f"max_channel {self.max_channel} can't be larger than"
                f" channel_num {self.channel_num} (flowcell {self.flowcell_type}). "
                "Please check your config."
            )
        self._check_for_unknown_keys(config)


class MainConfig:
    def __init__(self, config_toml_path: str):
        with open(config_toml_path, "r") as f:
            run_config = toml.load(f)

        # make sure all keys are present
        for key in [
            "processing",
            "acquisition",
            "balancing",
            "reporting",
            "model",
        ]:
            if key == "balancers":
                run_config[key] = run_config.get(key, [{}])
            else:
                run_config[key] = run_config.get(key, {})

        self.model: "ModelConfig" = ModelConfig(run_config)
        self.flowcell: "FlowcellConfig" = FlowcellConfig(run_config)
        self.processing: "ProcessingConfig" = ProcessingConfig(run_config)
        self.acquisition: "AcquisitionConfig" = AcquisitionConfig(run_config)
        self.balancing: "BalancingConfig" = BalancingConfig(run_config)
        self.reporting: "ReportingConfig " = ReportingConfig(run_config)

        self.channels = np.arange(
            self.flowcell.min_channel, self.flowcell.max_channel + 1
        ).tolist()

        self.balancers: List["BalancerConfig"] = self._create_balancers(run_config)

    def _create_balancers(self, run_config: dict):
        channels = (
            np.random.permutation(self.channels).copy().tolist()
        )  # will get popped from
        n_all_channels = len(channels)  # save number of channels for later

        balancers = [
            BalancerConfig(
                balancer_config,
                num_bcs=self.model.num_bcs,
                n_all_channels=n_all_channels,
                available_channels=channels,
            )
            for balancer_config in run_config["balancers"]
        ]

        # check for remaining channels, if any
        #   if there are, and there is a 'none' balancer, add them to the 'none' balancer
        #   if there are, and there is no 'none' balancer, create a 'none' balancer with the remaining channels
        if len(channels):
            channel_frac = len(channels) / n_all_channels
            none_balancer = None
            for balancer_config in balancers:
                # find first none balancer
                if balancer_config.balance_type == "none":
                    none_balancer = balancer_config
                    break

            if none_balancer is None:
                balancers.append(
                    BalancerConfig(
                        {
                            "balance_threshold": DEFAULT_BALANCE_THRESHOLD,
                            "min_stat": DEFAULT_MIN_STAT,
                            "balance_type": "none",
                            "name": "unused_channels",
                            "pred_conf_threshold": DEFAULT_PRED_CONF_THRESHOLD,
                            "channel_frac": channel_frac,
                            "pod5_watch_dir": None,
                            "pod5_check_interval": None,
                            "channels": channels,
                        },
                        self.model.num_bcs,
                        n_all_channels,
                        channels,
                    )
                )
            else:
                none_balancer.channels = np.concatenate(
                    (none_balancer.channels, channels)
                )
                none_balancer.channel_frac += channel_frac

        # Having duplicate balance types is not allowed
        balance_names = [balancer_config.name for balancer_config in balancers]
        if len(balance_names) != len(set(balance_names)):
            raise ValueError(
                f"Duplicate balancer found in config: {balance_names}. Please check"
                " your config. When using multiple balancers of the same balance_type,"
                " make sure to give each an unique name!"
            )
        return balancers

    def to_dict(self):
        # TODO: partially move this to sub config classes
        return {
            "processing": self.processing.to_dict(),
            "acquisition": self.acquisition.to_dict(),
            "model": self.model.to_dict(),
            "balancing": self.balancing.to_dict(),
            "reporting": self.reporting.to_dict(),
            "balancers": [balancer.to_dict() for balancer in self.balancers],
        }
