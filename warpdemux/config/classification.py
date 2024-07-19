"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass

from adapted.config.base import BaseConfig


@dataclass
class ClassificationConfig(BaseConfig):
    model_name: str = "WDX12"
    block_size: int = 1000
    num_proc: int = -1
