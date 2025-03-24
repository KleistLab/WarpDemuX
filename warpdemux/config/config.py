"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass, field

from adapted.config.config import Config as AdaptedConfig

from warpdemux.config.classification import ClassificationConfig
from warpdemux.config.file_proc import (BatchConfig, InputConfig, OutputConfig,
                                        TaskConfig)
from warpdemux.config.sig_proc import SigProcConfig


@dataclass
class Config(AdaptedConfig):
    input: InputConfig = field(default_factory=InputConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    batch: BatchConfig = field(default_factory=BatchConfig)
    task: TaskConfig = field(default_factory=TaskConfig)
    sig_proc: SigProcConfig = field(default_factory=SigProcConfig)
    classif: ClassificationConfig = field(default_factory=ClassificationConfig)
