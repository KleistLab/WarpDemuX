"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass
from typing import Tuple, Union

from adapted.config.base import BaseConfig as AdaptedBaseConfig
from adapted.config.sig_proc import SigProcConfig as AdaptedSigProcConfig


@dataclass
class SigExtractConfig(AdaptedBaseConfig):
    padding: int = 5
    normalization: str = "none"


@dataclass
class SegmentationConfig(AdaptedBaseConfig):
    """Segmentation parameters for barcode detection.

    Parameters
    ----------
    min_obs_per_base: int
        Minimum number of observations per base.
    running_stat_width: int
        Width of the running statistic window.
    num_events: int
        Number of events to detect.
    accept_less_cpts: bool
        Whether to accept less cpts than num_events.

    consensus_refinement: bool
        Whether to use consensus barcode refinement.
    consensus_model: str
        The consensus model to use. See warpdemux.config._consensus
    refinement_optimal_cpts: bool
        Whether to use the optimal segmentation for the barcode segmentation.

    normalization: str
        Normalization method. Can be "none", "mean", or "median".
    barcode_num_events: Union[int, Tuple[int, int]]
        Number of events to keep for barcode detection.
        If a tuple is provided, the first element is the number of barcode events to detect,
        and the second element is the number of barcode events to keep.
    """

    min_obs_per_base: int = 15
    running_stat_width: int = 30
    num_events: int = 110
    accept_less_cpts: bool = False

    consensus_refinement: bool = False
    consensus_model: str = ""

    refinement_optimal_cpts: bool = False

    normalization: str = "none"
    barcode_num_events: Union[int, Tuple[int, int]] = 25


@dataclass
class SigProcConfig(AdaptedSigProcConfig):
    sig_extract: SigExtractConfig = SigExtractConfig()
    segmentation: SegmentationConfig = SegmentationConfig()
