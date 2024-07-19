"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass

from adapted.config.base import BaseConfig as AdaptedBaseConfig
from adapted.config.sig_proc import (
    SigProcConfig as AdaptedSigProcConfig,
)


@dataclass
class SigNormConfig(AdaptedBaseConfig):
    winsor_window: int = 5
    outlier_thresh: float = 5.0


@dataclass
class SigExtractConfig(AdaptedBaseConfig):
    padding: int = 5
    normalization: str = "none"


@dataclass
class SegmentationConfig(AdaptedBaseConfig):
    min_obs_per_base: int = 15
    running_stat_width: int = 30
    num_events: int = 110
    accept_less_cpts: bool = False
    pad_if_less: bool = True

    normalization: str = "none"
    barcode_num_events: int = 25


@dataclass
class SigProcConfig(AdaptedSigProcConfig):
    sig_norm: SigNormConfig = SigNormConfig()
    sig_extract: SigExtractConfig = SigExtractConfig()
    segmentation: SegmentationConfig = SegmentationConfig()
