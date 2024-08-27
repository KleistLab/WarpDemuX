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

import numpy as np
import pyximport

pyximport.install(setup_args={"include_dirs": np.get_include()})
from ._c_segmentation import c_new_means, c_windowed_t_test


def windowed_t_test(
    raw_signal: np.ndarray,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
) -> np.ndarray:
    try:
        t_scores = c_windowed_t_test(
            raw_signal.astype(np.float64),
            min_obs_per_base,
            running_stat_width,
        )
        return t_scores

    except Exception as e:
        print(e)
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


def identify_stalls(
    all_raw_signal,
    window_size=5 * 25,
    threshold=2,
    edge_buffer=100,
    min_consecutive_obs=200,
    n_windows=5,
    mini_window_size=25,
    return_metric=False,
):
    """Identify locations where bases have stalled in the pore.
    source: tombo/tombo_stats.py
    """

    def compute_running_mean_diffs():
        """Compute average difference between n_window neighboring window means
        each of size window_size.
        """
        moving_average = np.cumsum(all_raw_signal)
        moving_average[mini_window_size:] = (
            moving_average[mini_window_size:] - moving_average[:-mini_window_size]
        )
        moving_average = moving_average[mini_window_size - 1 :] / mini_window_size

        # extract moving window averages at n_window offsets
        offsets = [
            moving_average[
                int(mini_window_size * offset) : int(
                    -mini_window_size * (n_windows - offset - 1)
                )
            ]
            for offset in range(n_windows - 1)
        ] + [
            moving_average[int(mini_window_size * (n_windows - 1)) :],
        ]
        # compute difference between all pairwise offset
        diffs = [
            np.abs(offsets[i] - offsets[j])
            for i in range(n_windows)
            for j in range(i + 1, n_windows)
        ]

        # compute average over offset differences at each valid position
        diff_sums = diffs[0].copy()
        for diff_i in diffs:
            diff_sums += diff_i
        return diff_sums / len(diffs)

    # if the raw signal is too short to compute stall metrics
    if all_raw_signal.shape[0] < window_size:
        if return_metric:
            return [], np.repeat(np.NAN, all_raw_signal.shape[0])
        return []

    # identify potentially stalled signal from either running window means
    # or running percentile difference methods
    stall_metric = np.empty(all_raw_signal.shape, all_raw_signal.dtype)
    stall_metric[:] = np.NAN
    start_offset = int(window_size * 0.5)
    end_offset = all_raw_signal.shape[0] - window_size + start_offset + 1

    assert window_size == mini_window_size * n_windows
    stall_metric[start_offset:end_offset] = compute_running_mean_diffs()

    # identify contiguous windows over threshold for minimal stretches
    with np.errstate(invalid="ignore"):
        stall_locs = np.where(
            np.diff(np.concatenate([[False], stall_metric <= threshold]))
        )[0]
    if stall_metric[-1] <= threshold:
        stall_locs = np.concatenate([stall_locs, [stall_metric.shape[0]]])
    stall_locs = stall_locs.reshape(-1, 2)
    stall_locs = stall_locs[(np.diff(stall_locs) > min_consecutive_obs).flatten()]
    if stall_locs.shape[0] == 0:
        if return_metric:
            return [], stall_metric
        return []

    # expand windows out to region that gave result below threshold
    # since windows are centered (minus edge buffer)
    expand_width = (window_size // 2) - edge_buffer
    if expand_width > 0:
        stall_locs[:, 0] -= expand_width
        stall_locs[:, 1] += expand_width
        # collapse intervals that now overlap
        merged_stall_locs = []
        prev_int = stall_locs[0]
        for curr_int in stall_locs:
            if curr_int[0] > prev_int[1]:
                # add previous interval to all intervals
                merged_stall_locs.append(prev_int)
                prev_int = curr_int
            else:
                # extend previous interval since these overlap
                prev_int[1] = curr_int[1]
        merged_stall_locs.append(prev_int)
        stall_locs = merged_stall_locs

    if return_metric:
        return stall_locs, stall_metric
    return stall_locs


def remove_stall_cpts(stall_ints, valid_cpts):
    """Remove stall points from valid points

    source: tombo/tombo_stats.py
    """
    if len(stall_ints) == 0:
        return valid_cpts

    # RNA data contains stall regions that can cause problems for
    # banded dynamic programming so they are removed here
    stall_int_iter = iter(stall_ints)
    curr_stall_int = next(stall_int_iter)
    non_stall_cpts = []
    # loop over valid cpts
    for i, cpt in enumerate(valid_cpts):
        # iterate through stall intervals until the current interval end
        # is greater than the cpt to check against
        while cpt > curr_stall_int[1]:
            try:
                curr_stall_int = next(stall_int_iter)
            except StopIteration:
                break
        if not (curr_stall_int[0] < cpt < curr_stall_int[1]):
            non_stall_cpts.append(i)

    return valid_cpts[non_stall_cpts]
