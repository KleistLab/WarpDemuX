"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Union

import numpy as np
from adapted.config.base import BaseConfig


@dataclass
class OutputConfig(BaseConfig):
    output_dir: str = ""

    save_llr_trace: bool = False
    save_dwell_time: bool = False
    save_boundaries: bool = True


@dataclass
class InputConfig(BaseConfig):
    files: List[str] = field(default_factory=list)
    read_ids: Union[List[str], np.ndarray] = field(default_factory=list)

    preprocessed: bool = False


@dataclass
class ResegmentTaskConfig(BaseConfig):
    detect_result_dict: Dict[str, Dict[str, Any]] = field(default_factory=dict)
