"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import os
from dataclasses import dataclass, field
from typing import List, Union

import numpy as np
from adapted.config.base import BaseConfig


@dataclass
class OutputConfig(BaseConfig):
    output_dir: str = ""

    save_dwell_time: bool = False
    save_fpts: bool = False
    save_boundaries: bool = True

    output_subdir_pred: str = "predictions"
    output_subdir_fail: str = "failed_reads"
    output_subdir_fpts: str = "fingerprints"
    output_subdir_boundaries: str = "boundaries"

    def __post_init__(self):
        self.output_dir_pred = os.path.join(self.output_dir, self.output_subdir_pred)
        self.output_dir_fail = os.path.join(self.output_dir, self.output_subdir_fail)
        self.output_dir_fpts = os.path.join(self.output_dir, self.output_subdir_fpts)
        self.output_dir_boundaries = os.path.join(
            self.output_dir, self.output_subdir_boundaries
        )

        # create output directories if valid
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)
            os.makedirs(self.output_dir_pred, exist_ok=True)
            os.makedirs(self.output_dir_fail, exist_ok=True)

            if self.save_boundaries:
                os.makedirs(self.output_dir_boundaries, exist_ok=True)
            if self.save_fpts:
                os.makedirs(self.output_dir_fpts, exist_ok=True)


@dataclass
class InputConfig(BaseConfig):
    files: List[str] = field(default_factory=list)
    read_ids: Union[List[str], np.ndarray] = field(default_factory=list)
    continue_from: str = ""
    retry_from: str = ""
    n_reads: int = -1
    preprocessed: bool = False


@dataclass
class BatchConfig(BaseConfig):
    num_proc: int = -1  # default to number of cores
    batch_size_output: int = 4000
    minibatch_size: int = 1000
    bidx_pass: int = 0
    bidx_fail: int = 0


@dataclass
class TaskConfig(BaseConfig):
    pass
