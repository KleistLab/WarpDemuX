"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import time

import attrs
import numpy as np

BALANCE_TYPE = [
    "none",
    "adapter_count",
    "read_count",
    "base_normalization",
    "block_safety",
    "reject_all",
]

BALANCE_TYPE_WATCHER = ["read_count", "base_normalization", "block_safety"]


@attrs.define
class ReadObject:
    """Read object passed between the different workers.
    data_arr is updated at each step.
    """

    channel: int
    read_start_sample: int
    chunk_start_sample: int
    read_number: int
    read_id: str
    chunk_length: int
    polya_start: int
    data_arr: (
        np.ndarray
    )  # first: raw signal, then segment means, then FPT, then prob_dist
    start_time: float
    time_per_step: (
        []
    )  # time_polya_detect, time_segmentation, time_classification, time_balancing
    is_outlier: bool = attrs.field(default=False, init=False)


@attrs.define
class Result:
    """Result holder

    barcode_label is 'unclassified' if the read is unassigned based on the conf_threshold
    """

    channel: int
    read_start_sample: int
    chunk_start_sample: int
    read_number: int
    read_id: str
    chunk_length: int
    polya_start: int
    prob_dist: np.ndarray
    is_outlier: bool
    barcode_label: str
    decision: str
    reason: str
    balance_type: str
    balancer_name: str
    time_polya_detect: float
    time_segmentation: float
    time_classification: float
    time_balancing: float
    total_compute_time: float
    total_time: float

    def __init__(
        self,
        read_obj: ReadObject,
        label: str,
        decision: str,
        reason: str,
        balancer_name: str,
        balance_type: str,
    ):
        self.channel = read_obj.channel
        self.read_start_sample = read_obj.read_start_sample
        self.chunk_start_sample = read_obj.chunk_start_sample
        self.read_number = read_obj.read_number
        self.read_id = read_obj.read_id
        self.chunk_length = read_obj.chunk_length
        self.polya_start = read_obj.polya_start
        self.prob_dist = read_obj.data_arr
        self.is_outlier = read_obj.is_outlier
        self.barcode_label = label
        self.decision = decision
        self.reason = reason
        self.balancer_name = balancer_name
        self.balance_type = balance_type
        self.time_polya_detect = read_obj.time_per_step[0]
        self.time_segmentation = read_obj.time_per_step[1]
        self.time_classification = read_obj.time_per_step[2]
        self.time_balancing = read_obj.time_per_step[3]
        self.total_compute_time = sum(read_obj.time_per_step)
        self.total_time = time.time() - read_obj.start_time

    def to_dict(self):
        return {
            "channel": self.channel,
            "read_start_sample": self.read_start_sample,
            "chunk_start_sample": self.chunk_start_sample,
            "read_number": self.read_number,
            "read_id": self.read_id,
            "chunk_length": self.chunk_length,
            "polya_start": self.polya_start,
            "barcode_label": self.barcode_label,
            "is_outlier": self.is_outlier,
            "decision": self.decision,
            "reason": self.reason,
            "balancer_name": self.balancer_name,
            "balance_type": self.balance_type,
            **{f"p{i:02d}": p for i, p in enumerate(self.prob_dist.squeeze()[:-1])},
            "p-1": self.prob_dist.squeeze()[-1],
            "time_polya_detect": self.time_polya_detect,
            "time_segmentation": self.time_segmentation,
            "time_classification": self.time_classification,
            "time_balancing": self.time_balancing,
            "total_compute_time": self.total_compute_time,
            "total_time": self.total_time,
        }


@attrs.define
class FailedRead(ReadObject):
    """Read object passed between the different workers.
    data_arr is updated at each step.
    """

    channel: int
    read_start_sample: int
    chunk_start_sample: int
    read_number: int
    read_id: str
    chunk_length: int
    polya_start: int
    data_arr: np.ndarray
    start_time: float
    time_per_step: []
    decision: str
    reason: str

    def __init__(self, reason, decision, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reason = reason
        self.decision = decision

        self.time_per_step = self.time_per_step + [0] * (4 - len(self.time_per_step))

    def to_dict(self):
        return {
            "channel": self.channel,
            "read_start_sample": self.read_start_sample,
            "chunk_start_sample": self.chunk_start_sample,
            "read_number": self.read_number,
            "read_id": self.read_id,
            "chunk_length": self.chunk_length,
            "polya_start": self.polya_start,
            "decision": self.decision,
            "reason": self.reason,
        }

    def to_result(self, num_bcs):
        self.data_arr = np.zeros((num_bcs, 1))
        return Result(
            self,
            label="failed_read",
            reason=self.reason,
            decision=self.decision,
            balance_type="failed_read",
            balancer_name="failed_read",  # failed reads get their own 'balancer';
            # these reads are not classified, the're kicked out before balancing
        )
