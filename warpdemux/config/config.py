"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass
from typing import Optional, Union

from adapted.config.config import Config as AdaptedConfig
from adapted.config.file_proc import (
    BatchConfig,
    TaskConfig as DetectTaskConfig,
)

from warpdemux.config.classification import ClassificationConfig
from warpdemux.config.file_proc import (
    InputConfig,
    OutputConfig,
    ResegmentTaskConfig,
)
from warpdemux.config.sig_proc import SigProcConfig


@dataclass
class Config(AdaptedConfig):
    input: InputConfig = InputConfig()
    output: OutputConfig = OutputConfig()
    task: Union[DetectTaskConfig, ResegmentTaskConfig] = DetectTaskConfig()
    batch: BatchConfig = BatchConfig()
    sig_proc: SigProcConfig = SigProcConfig()
    classif: Optional[ClassificationConfig] = None
