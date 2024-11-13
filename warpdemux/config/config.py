"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass
from typing import Optional, Union

from adapted.config.config import Config as AdaptedConfig

from warpdemux.config.classification import ClassificationConfig
from warpdemux.config.file_proc import (
    BatchConfig,
    InputConfig,
    OutputConfig,
    TaskConfig,
)
from warpdemux.config.sig_proc import SigProcConfig


@dataclass
class Config(AdaptedConfig):
    input: InputConfig = InputConfig()
    output: OutputConfig = OutputConfig()
    batch: BatchConfig = BatchConfig()
    task: TaskConfig = TaskConfig()
    sig_proc: SigProcConfig = SigProcConfig()
    classif: ClassificationConfig = ClassificationConfig()
