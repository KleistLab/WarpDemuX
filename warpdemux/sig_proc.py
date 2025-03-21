"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, TypeVar, Union, overload

import numpy as np
from adapted.container_types import DetectResults
from adapted.container_types import ReadResult as AdaptedReadResult
from dtaidistance.dtw import warping_paths_fast
from dtaidistance.subsequence import SubsequenceAlignment
from ruptures import KernelCPD
from scipy.signal import find_peaks

from warpdemux.config.sig_proc import SegmentationConfig, SigProcConfig
from warpdemux.segmentation.segmentation import (compute_base_means,
                                                 windowed_t_test)


@dataclass
class ReadResult(AdaptedReadResult):
    # read_id
    # success
    # fail_reason
    # detect_results
    barcode_fpt: Optional[np.ndarray] = None
    dwell_times: Optional[np.ndarray] = None
    adapter_dt_med: Optional[float] = None
    adapter_dt_mad: Optional[float] = None
    adapter_event_mean: Optional[float] = None
    adapter_event_std: Optional[float] = None
    adapter_event_med: Optional[float] = None
    adapter_event_mad: Optional[float] = None
    seg_cons_query_start: Optional[int] = None
    seg_cons_query_end: Optional[int] = None
    sig_barcode_start: Optional[int] = None

    def to_summary_dict(self) -> Dict[str, Any]:
        summary_dict = super().to_summary_dict()
        summary_dict.update(
            {
                "adapter_dt_med": self.adapter_dt_med,
                "adapter_dt_mad": self.adapter_dt_mad,
                "adapter_event_mean": self.adapter_event_mean,
                "adapter_event_std": self.adapter_event_std,
                "adapter_event_med": self.adapter_event_med,
                "adapter_event_mad": self.adapter_event_mad,
                "seg_cons_query_start": self.seg_cons_query_start,
                "seg_cons_query_end": self.seg_cons_query_end,
                "sig_barcode_start": self.sig_barcode_start,
            }
        )
        return summary_dict

    def set_read_id(self, read_id: str):
        self.read_id = read_id


###############
# Normalization
###############


def mad_normalize(signal: np.ndarray, accept_nan: bool = False) -> np.ndarray:
    """
    Normalize a signal using Median Absolute Deviation (MAD).

    Parameters:
    ----------
    signal : np.ndarray
        The input signal array to be normalized.

    Returns:
    -------
    np.ndarray
        The normalized signal array.
    """

    if np.isnan(signal).any() and not accept_nan:
        msg = "Signal contains NaN values."
        logging.error(msg)
        raise ValueError(msg)
    elif np.isnan(signal).any():
        shift = np.nanmedian(signal, axis=-1, keepdims=True)
        scale = np.nanmedian(np.abs(signal - shift), axis=-1, keepdims=True)
    else:
        shift = np.median(signal, axis=-1, keepdims=True)
        scale = np.median(np.abs(signal - shift), axis=-1, keepdims=True)

    return (signal - shift) / scale


def mean_normalize(signal: np.ndarray, accept_nan: bool = False) -> np.ndarray:
    if np.isnan(signal).any() and not accept_nan:
        msg = "Signal contains NaN values."
        logging.error(msg)
        raise ValueError(msg)
    elif np.isnan(signal).any():
        shift = np.nanmean(signal, axis=-1, keepdims=True)
        scale = np.nanstd(signal, axis=-1, keepdims=True)
    else:
        shift = np.mean(signal, axis=-1, keepdims=True)
        scale = np.std(signal, axis=-1, keepdims=True)

    return (signal - shift) / scale


def normalize(
    signal: np.ndarray, method: str = "mean", accept_nan: bool = False
) -> np.ndarray:
    if signal.size == 0:
        return signal

    if np.isnan(signal).any() and not accept_nan:
        msg = "Signal contains NaN values."
        logging.error(msg)
        raise ValueError(msg)

    if method == "mean":
        norm_signal = mean_normalize(signal, accept_nan=accept_nan)
    elif method == "median":
        norm_signal = mad_normalize(signal, accept_nan=accept_nan)
    elif method == "none":
        norm_signal = signal
    else:
        msg = f"Normalization method {method} not recognized."
        logging.error(msg)
        raise ValueError(msg)

    return norm_signal


def normalize_wrt(
    to_norm: np.ndarray, ref: np.ndarray, method: str = "mean"
) -> np.ndarray:
    """
    Normalize a signal using the mean and standard deviation (or median and MAD) computed from a reference signal.
    The normalization is applied over the last dimension of the input array.

    Args:
        to_norm: Signal to normalize, can be 1D (M1,) or 2D (N, M1)
        ref: Reference signal, must be 1D (M2,)
        method: Normalization method, either "mean" or "median"

    Returns:
        Normalized signal with same shape as input to_norm
    """
    if to_norm.ndim == 1:
        to_norm = to_norm[:, None]

    if method == "mean":
        norm_signal = ((to_norm - np.mean(ref)) / np.std(ref)).squeeze(-1)
    elif method == "median":
        med = np.median(ref)
        mad = np.median(np.abs(ref - med))
        norm_signal = ((to_norm - med) / mad).squeeze(-1)
    else:
        msg = f"Normalization method {method} not recognized."
        logging.error(msg)
        raise ValueError(msg)

    return norm_signal


###############
# Adapter
###############


def discrepenacy_curve_to_cpts(
    scores: np.ndarray,
    num_events: int = 110,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    accept_less_cpts=False,
) -> np.ndarray:
    peaks, _ = find_peaks(scores, distance=min_obs_per_base)

    if peaks.size < num_events and not accept_less_cpts:
        return np.array([])

    valid_cpts = peaks[np.argsort(scores[peaks])[-num_events:]] + running_stat_width
    valid_cpts.sort()

    signal_len = scores.size + 2 * running_stat_width

    if valid_cpts[0] != 0:
        valid_cpts = np.insert(valid_cpts, 0, 0)
    if valid_cpts[-1] != signal_len:
        valid_cpts = np.append(valid_cpts, signal_len)

    return valid_cpts


def segment_signal(
    raw_signal: np.ndarray,
    num_events: int = 110,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    accept_less_cpts=False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Segment the given raw signal using the `c_windowed_t_test` function.

    If segmentation is not possible due to a NotImplementedError, an empty array is returned.

    Parameters:
    ----------
    raw_signal : np.ndarray
        The raw signal array to be segmented.

    min_obs_per_base : int, optional
        Minimum observations per base.

    running_stat_width : int, optional
        Width for the running statistic.

    num_events : int, optional
        The desired number of events.

    Returns:
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        - Array containing segmented, normalized signal values.
        - Array containing dwell times.
        - Array containing windowed t-test scores.
    """

    scores = windowed_t_test(
        raw_signal,
        running_stat_width=running_stat_width,
    )

    valid_cpts = discrepenacy_curve_to_cpts(
        scores,
        num_events=num_events,
        min_obs_per_base=min_obs_per_base,
        running_stat_width=running_stat_width,
        accept_less_cpts=accept_less_cpts,
    )

    if valid_cpts.size == 0:
        return np.array([]), np.array([]), np.array([])

    dwell_times = valid_cpts[1:] - valid_cpts[:-1]
    event_means = compute_base_means(raw_signal, valid_cpts)

    return event_means, dwell_times, scores


def segment_signal_with_consensus_guided_barcode_refinement(
    raw_signal: np.ndarray,
    consensus_query: np.ndarray,
    config: SegmentationConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int, int]:
    """
    Segment the given raw signal using the `c_windowed_t_test` function.

    Parameters
    ----------
    raw_signal : np.ndarray
        The raw signal array to be segmented.
    consensus_query : np.ndarray
        The consensus signal of the constant adapter sequence to be used for determining the barcode start position in the segmented signal.
        The consensus signal should be normalized to the same scale as the raw signal.
    config : SegmentationConfig
        The segmentation configuration to be used for the barcode refinement.


    Returns
    -------
    Tuple[np.ndarray, np.ndarray, int, int, int]
        - Normalized, consensus-guided barcode-refined, segmented signal values (mean value per event)
        - Dwell times
        - Normalized adapter segmentation means
        - Adapter dwell times
        - Segment query start index
        - Segment query end index
        - Signal barcode start index
    """

    def _get_subseq_match(query, series):

        norm_series = normalize(series, config.consensus_subseq_match_normalization)
        sa = SubsequenceAlignment(
            query,
            norm_series,
            use_c=True,
        )
        # NOTE: align outside of class to control psi and penalty
        _, sa.paths = warping_paths_fast(
            query,
            norm_series,
            penalty=config.consensus_subseq_match_penalty,
            psi=config.consensus_subseq_match_psi,
            compact=False,
            psi_neg=False,
        )
        sa._compute_matching()
        match = sa.best_match()
        startidx, endidx = match.segment
        return startidx, endidx

    adapter_event_means, adapter_dwell_times, adapter_scores = segment_signal(
        raw_signal,
        num_events=config.num_events,
        min_obs_per_base=min(
            config.min_obs_per_base,
            int(round(raw_signal.size / config.num_events / 2)),
        ),
        running_stat_width=min(
            config.running_stat_width,
            int(round(raw_signal.size / config.num_events)),
        ),
        accept_less_cpts=config.accept_less_cpts,
    )

    if adapter_event_means.size == 0:
        return np.array([]), np.array([]), np.array([]), np.array([]), 0, 0, 0

    # seg_barcode_start
    seg_query_start, seg_query_end = _get_subseq_match(
        consensus_query, adapter_event_means
    )
    sig_barcode_start = int(np.sum(adapter_dwell_times[:seg_query_end]))

    barcode_scores = adapter_scores[sig_barcode_start:]

    # if tuple, first entry is segmentation num_events, second is num_events to keep
    if hasattr(config.barcode_num_events, '__len__'):
        if not isinstance(config.barcode_num_events, (tuple, list, np.ndarray)):
            msg = "barcode_num_events must be a tuple, list or numpy array when using multiple values"
            logging.error(msg)
            raise TypeError(msg)
        barcode_segm_num_events = int(config.barcode_num_events[0])
    else:
        assert isinstance(config.barcode_num_events, int)
        barcode_segm_num_events = int(config.barcode_num_events)

    # logging.info(f"seg_query_end: {seg_query_end}")

    if config.refinement_optimal_cpts:
        algo = KernelCPD(kernel="linear", min_size=config.min_obs_per_base)
        algo.fit(barcode_scores.reshape(-1, 1))
        valid_cpts = np.array(
            algo.predict(n_bkps=barcode_segm_num_events)
        )  # size num_events_barcode+1, includes signal length
        valid_cpts = np.insert(valid_cpts, 0, 0)
    else:
        valid_cpts = discrepenacy_curve_to_cpts(
            barcode_scores,
            num_events=barcode_segm_num_events,
            min_obs_per_base=config.min_obs_per_base,
            running_stat_width=config.running_stat_width,
            accept_less_cpts=False,
        )
    
    if valid_cpts.size == 0:
        return np.array([]), np.array([]), np.array([]), np.array([]), 0, 0, 0

    barcode_dwell_times = valid_cpts[1:] - valid_cpts[:-1]
    barcode_event_means = compute_base_means(raw_signal[sig_barcode_start:], valid_cpts)

    return (
        barcode_event_means,
        barcode_dwell_times,
        adapter_event_means,
        adapter_dwell_times,
        seg_query_start,
        seg_query_end,
        sig_barcode_start,
    )


# TODO: define end and begin extract_padding?
def extract_adapter(
    signal,
    adapter_start: int,
    adapter_end: int,
    extract_padding: int,
) -> np.ndarray:
    start = max(0, adapter_start - extract_padding)
    stop = min(signal.size, adapter_end + extract_padding)

    return signal[start:stop]


def detect_results_to_fpt(
    calibrated_signal: np.ndarray,
    spc: "SigProcConfig",
    detect_results: DetectResults,
    consensus_query: np.ndarray = np.array([]),
) -> ReadResult:
    if not detect_results.success:
        return ReadResult(
            success=False,
            fail_reason=detect_results.fail_reason,
            barcode_fpt=np.array([]),
            dwell_times=np.array([]),
            detect_results=detect_results,
        )

    assert (
        detect_results.adapter_start is not None
        and detect_results.adapter_end is not None
    )  # make type checker happy

    adapter_sig = extract_adapter(
        calibrated_signal,
        adapter_start=detect_results.adapter_start,
        adapter_end=detect_results.adapter_end,
        extract_padding=spc.sig_extract.padding,
    )

    med = np.nanmedian(adapter_sig)
    mad = np.nanmedian(np.abs(adapter_sig - med))

    # impute mad outliers, required: signal flickers strongly affect segmentation results
    # in place
    np.clip(
        adapter_sig,
        med - spc.core.sig_norm_outlier_thresh * mad,
        med + spc.core.sig_norm_outlier_thresh * mad,
        out=adapter_sig,
    )

    try:
        adapter_sig = normalize(
            adapter_sig,
            spc.sig_extract.normalization,
            accept_nan=True,
        )  # can be used to normalize by adapter mean or median prior to segmentation
    except Exception as e:
        return ReadResult(
            success=False,
            fail_reason=f"signal normalization failed: {e}",
            barcode_fpt=np.array([]),
            dwell_times=np.array([]),
            detect_results=detect_results,
        )

    seg_cons_query_start = None
    seg_cons_query_end = None
    sig_barcode_start = None

    if spc.segmentation.consensus_refinement:
        # TODO: move to config validation
        if type(spc.segmentation.barcode_num_events) == int:
            raise ValueError(
                "barcode_num_events is an integer in consensus refinement mode, "
                "use a tuple instead"
            )

        (
            barcode_event_means,
            barcode_dwell_times,
            adapter_event_means,
            adapter_dwell_times,
            seg_cons_query_start,
            seg_cons_query_end,
            sig_barcode_start,
        ) = segment_signal_with_consensus_guided_barcode_refinement(
            adapter_sig,
            consensus_query,
            config=spc.segmentation,
        )
        if not barcode_event_means.size:
            return ReadResult(
                success=False,
                fail_reason="event segmentation failed",
                barcode_fpt=np.array([]),
                dwell_times=np.array([]),
                detect_results=detect_results,
            )

        norm_barcode_event_means = normalize_wrt(
            barcode_event_means, adapter_event_means, spc.segmentation.normalization
        )

        # rename for consistency with non-refined segmentation
        norm_segment_avgs = norm_barcode_event_means
        dwell_times = barcode_dwell_times

        adapter_dt_med = float(np.median(adapter_dwell_times))
        adapter_dt_mad = float(np.median(np.abs(adapter_dwell_times - adapter_dt_med)))

        adapter_event_mean = float(adapter_event_means.mean())
        adapter_event_std = float(adapter_event_means.std())
        adapter_event_med = float(np.median(adapter_event_means))
        adapter_event_mad = float(
            np.median(np.abs(adapter_event_means - adapter_event_med))
        )

        # filter out consensus query outliers
        if (
            seg_cons_query_start > spc.segmentation.consensus_subseq_match_ub_start
            or seg_cons_query_end < spc.segmentation.consensus_subseq_match_lb_end
            or seg_cons_query_end > spc.segmentation.consensus_subseq_match_ub_end
        ):
            return ReadResult(
                success=False,
                fail_reason="consensus query outlier",
                barcode_fpt=np.array([]),
                dwell_times=np.array([]),
                detect_results=detect_results,
                adapter_dt_med=adapter_dt_med,
                adapter_dt_mad=adapter_dt_mad,
                adapter_event_mean=adapter_event_mean,
                adapter_event_std=adapter_event_std,
                adapter_event_med=adapter_event_med,
                adapter_event_mad=adapter_event_mad,
                seg_cons_query_start=seg_cons_query_start,
                seg_cons_query_end=seg_cons_query_end,
                sig_barcode_start=sig_barcode_start,
            )
    else:
        segment_avgs, dwell_times, _ = segment_signal(
            adapter_sig,
            num_events=spc.segmentation.num_events,
            min_obs_per_base=min(
                spc.segmentation.min_obs_per_base,
                round(adapter_sig.size / spc.segmentation.num_events / 2),
            ),
            running_stat_width=min(
                spc.segmentation.running_stat_width,
                round(adapter_sig.size / spc.segmentation.num_events),
            ),
            accept_less_cpts=spc.segmentation.accept_less_cpts,
        )

        if not segment_avgs.size:
            return ReadResult(
                success=False,
                fail_reason="event segmentation failed",
                barcode_fpt=np.array([]),
                dwell_times=np.array([]),
                detect_results=detect_results,
            )

        try:
            # normalize the segment signal
            norm_segment_avgs = normalize(
                segment_avgs,
                spc.segmentation.normalization,
                accept_nan=False,
            )
        except Exception as e:
            return ReadResult(
                success=False,
                fail_reason=f"segment normalization failed: {e}",
                barcode_fpt=np.array([]),
                dwell_times=np.array([]),
                detect_results=detect_results,
            )

        adapter_dt_med = float(np.median(dwell_times))
        adapter_dt_mad = float(np.median(np.abs(dwell_times - adapter_dt_med)))
        adapter_event_mean = float(segment_avgs.mean())
        adapter_event_std = float(segment_avgs.std())
        adapter_event_med = float(np.median(segment_avgs))
        adapter_event_mad = float(np.median(np.abs(segment_avgs - adapter_event_med)))

    nsegm_retain = (
        spc.segmentation.barcode_num_events
        if not spc.segmentation.consensus_refinement
        else spc.segmentation.barcode_num_events[1]  # already confirmed as tuple
    )

    if norm_segment_avgs.size < nsegm_retain:
        # pad 3' end with NaNs
        norm_segment_avgs = np.pad(
            norm_segment_avgs,
            (spc.segmentation.barcode_num_events - norm_segment_avgs.size, 0),
            mode="constant",
            constant_values=np.nan,
        )
        dwell_times = np.pad(
            dwell_times,
            (spc.segmentation.barcode_num_events - dwell_times.size, 0),
            mode="constant",
            constant_values=np.nan,
        )

    return ReadResult(
        success=True,
        fail_reason="",
        barcode_fpt=norm_segment_avgs[-min(nsegm_retain, norm_segment_avgs.size) :],
        dwell_times=dwell_times[-min(nsegm_retain, dwell_times.size) :],
        detect_results=detect_results,
        adapter_dt_med=adapter_dt_med,
        adapter_dt_mad=adapter_dt_mad,
        adapter_event_mean=adapter_event_mean,
        adapter_event_std=adapter_event_std,
        adapter_event_med=adapter_event_med,
        adapter_event_mad=adapter_event_mad,
        seg_cons_query_start=seg_cons_query_start,
        seg_cons_query_end=seg_cons_query_end,
        sig_barcode_start=sig_barcode_start,
    )
