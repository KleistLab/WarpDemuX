"""
---------------------------------------------------------------------------------
DISCLAIMER:

The code contained within this file is based on, or directly taken from, the
Tombo software. Tombo is a software package provided by ONT Research and the
original code can be found at their official repository:

https://github.com/nanoporetech/tombo

All rights and acknowledgments go to the original authors and contributors of Tombo.
Any modifications made to the code from its original version, if any, are the
responsibility of the current file's maintainer and do not reflect the views or
practices of the original Tombo developers.

Use and distribution of this code should respect the licensing terms and conditions
set by the original Tombo developers.

This file is licensed under the Mozilla Public License 2.0 (MPL 2.0).
---------------------------------------------------------------------------------
"""

import logging

import numpy as np
import pyximport

pyximport.install(setup_args={"include_dirs": np.get_include()})
from ._c_segmentation import c_new_means, c_windowed_t_test


def windowed_t_test(
    raw_signal: np.ndarray,
    running_stat_width: int = 30,
) -> np.ndarray:
    try:
        t_scores = c_windowed_t_test(
            raw_signal.astype(np.float64),
            running_stat_width,
        )
        return t_scores

    except Exception as e:
        logging.error(e)
        return np.zeros(0, dtype=int)


def compute_base_means(raw_signal: np.ndarray, base_starts: np.ndarray) -> np.ndarray:
    """
    Compute mean base values from a raw signal based on start positions of bases.

    This is an updated version of the Tombo implementation. It ensures that segments
    at the beginning and end of the signal are also included.

    Parameters:
    ----------
    raw_signal : np.ndarray
        Raw nanopore signal observation values.

    base_starts : np.ndarray
        Array containing 0-based starting positions of bases within the raw signal.

    Returns:
    -------
    np.ndarray
        Array containing the mean values of bases.
    """

    if base_starts[0] != 0:
        base_starts = np.insert(base_starts, 0, 0)
    if base_starts[-1] != raw_signal.size:
        base_starts = np.append(base_starts, raw_signal.size)

    return c_new_means(raw_signal.astype(np.float64), base_starts)
