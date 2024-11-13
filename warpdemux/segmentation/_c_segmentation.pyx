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
cimport numpy as np

cdef bint boolean_variable = True

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

DTYPE_INT = np.int64
ctypedef np.int64_t DTYPE_INT_t

DTYPE_INT16 = np.int16
ctypedef np.int16_t DTYPE_INT16_t

cdef extern from "math.h":
    double sqrt(double m)


def c_new_means(
        np.ndarray[DTYPE_t] norm_signal not None,
        np.ndarray[DTYPE_INT_t] new_segs not None):
    cdef DTYPE_INT_t n_segs = new_segs.shape[0] - 1
    cdef np.ndarray[DTYPE_t] means_arr = np.empty(n_segs, dtype=DTYPE)
    cdef DTYPE_t curr_sum
    cdef DTYPE_INT_t idx, seg_idx
    for idx in range(n_segs):
        curr_sum = 0
        for seg_idx in range(new_segs[idx], new_segs[idx + 1]):
            curr_sum += norm_signal[seg_idx]
        means_arr[idx] = curr_sum / (new_segs[idx + 1] - new_segs[idx])
    return means_arr

def c_valid_cpts_w_cap_t_test(
        np.ndarray[DTYPE_t] raw_signal, DTYPE_INT_t min_base_obs,
        DTYPE_INT_t running_stat_width, DTYPE_INT_t num_cpts,
        bint accept_less_cpts,
        bint return_scores):
    cdef DTYPE_INT_t pos, idx
    cdef DTYPE_t pos_diff, m1, m2, var1, var2
    cdef DTYPE_INT_t num_cands = max( raw_signal.shape[0] - (running_stat_width * 2), 0)

    if num_cands == 0:
        return np.empty(0, dtype=DTYPE_INT)
        
    # note these will not actually be t-scores, but will be a monotonic transform
    # so the rank order will be the same
    cdef np.ndarray[DTYPE_t] t_scores = np.empty(num_cands, dtype=DTYPE)
    for pos in range(num_cands):
        # compute means
        m1 = 0
        for idx in range(running_stat_width):
            m1 += raw_signal[pos + idx]
        m1 /= running_stat_width
        m2 = 0
        for idx in range(running_stat_width):
            m2 += raw_signal[pos + running_stat_width + idx]
        m2 /= running_stat_width

        # compute sum of variances
        var1 = 0
        for idx in range(running_stat_width):
            pos_diff = raw_signal[pos + idx] - m1
            var1 += pos_diff * pos_diff
        var2 = 0
        for idx in range(running_stat_width):
            pos_diff = raw_signal[pos + running_stat_width + idx] - m2
            var2 += pos_diff * pos_diff

        if var1 + var2 == 0:
            t_scores[pos] = 0.0
        elif m1 > m2:
            t_scores[pos] = (m1 - m2) / sqrt(var1 + var2)
        else:
            t_scores[pos] = (m2 - m1) / sqrt(var1 + var2)

    cdef np.ndarray[DTYPE_INT_t] candidate_poss = np.argsort(
        t_scores).astype(DTYPE_INT)[::-1]

    cdef np.ndarray[DTYPE_INT_t] cpts = np.empty(num_cpts, dtype=DTYPE_INT)
    cpts[0] = candidate_poss[0] + running_stat_width
    blacklist_pos = set(range(
        candidate_poss[0] - min_base_obs + 1, candidate_poss[0] + min_base_obs))
    cdef DTYPE_INT_t cand_pos
    cdef DTYPE_INT_t cand_idx = 1
    cdef DTYPE_INT_t added_cpts = 1
    while added_cpts < num_cpts:
        cand_pos = candidate_poss[cand_idx]
        if cand_pos not in blacklist_pos:
            cpts[added_cpts] = cand_pos + running_stat_width
            added_cpts += 1
            blacklist_pos.update(range(
                cand_pos - min_base_obs + 1, cand_pos + min_base_obs))
        cand_idx += 1
        if cand_idx >= num_cands:
            if not accept_less_cpts:
                raise NotImplementedError('Fewer changepoints found than requested')
            else:
                return cpts[:added_cpts] if not return_scores else (cpts[:added_cpts], t_scores)
    return cpts if not return_scores else (cpts, t_scores)


def c_windowed_t_test(
        np.ndarray[DTYPE_t] raw_signal,
        DTYPE_INT_t running_stat_width):
    cdef DTYPE_INT_t pos, idx
    cdef DTYPE_t pos_diff, m1, m2, var1, var2
    cdef DTYPE_INT_t num_cands = raw_signal.shape[0] - (running_stat_width * 2)
    # note these will not actually be t-scores, but will be a monotonic transform
    # so the rank order will be the same
    cdef np.ndarray[DTYPE_t] t_scores = np.empty(num_cands, dtype=DTYPE)
    for pos in range(num_cands):
        # compute means
        m1 = 0
        for idx in range(running_stat_width):
            m1 += raw_signal[pos + idx]
        m1 /= running_stat_width
        m2 = 0
        for idx in range(running_stat_width):
            m2 += raw_signal[pos + running_stat_width + idx]
        m2 /= running_stat_width

        # compute sum of variances
        var1 = 0
        for idx in range(running_stat_width):
            pos_diff = raw_signal[pos + idx] - m1
            var1 += pos_diff * pos_diff
        var2 = 0
        for idx in range(running_stat_width):
            pos_diff = raw_signal[pos + running_stat_width + idx] - m2
            var2 += pos_diff * pos_diff

        if var1 + var2 == 0:
            t_scores[pos] = 0.0
        elif m1 > m2:
            t_scores[pos] = (m1 - m2) / sqrt(var1 + var2)
        else:
            t_scores[pos] = (m2 - m1) / sqrt(var1 + var2)

    return t_scores

def c_windowed_diff(np.ndarray[DTYPE_t] raw_signal, DTYPE_INT_t min_base_obs,
        DTYPE_INT_t running_stat_width, ):
    cdef np.ndarray[DTYPE_t] raw_cumsum = np.cumsum(
        np.concatenate([[0.0], raw_signal]))
    # get difference between all neighboring running_stat_width regions
    return np.abs(
        (2 * raw_cumsum[running_stat_width:-running_stat_width]) -
        raw_cumsum[:-2*running_stat_width] -
        raw_cumsum[2*running_stat_width:])