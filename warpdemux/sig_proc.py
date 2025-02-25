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
from dtaidistance.dtw import warping_path_fast
from scipy.signal import find_peaks

from warpdemux.config.sig_proc import SigProcConfig
from warpdemux.segmentation.segmentation import compute_base_means, windowed_t_test


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
    seg_barcode_start: Optional[int] = None
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
                "seg_barcode_start": self.seg_barcode_start,
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
        return ((to_norm - np.mean(ref)) / np.std(ref)).squeeze(-1)
    elif method == "median":
        med = np.median(ref)
        mad = np.median(np.abs(ref - med))
        return ((to_norm - med) / mad).squeeze(-1)
    else:
        raise ValueError(f"Normalization method {method} not recognized.")


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


np_ndarray = TypeVar("np_ndarray", bound=np.ndarray)


@overload
def segment_signal(
    raw_signal: np_ndarray,
    num_events: int = 110,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    accept_less_cpts: bool = False,
    return_scores: bool = False,
) -> Tuple[np_ndarray, np_ndarray]: ...


@overload
def segment_signal(
    raw_signal: np_ndarray,
    num_events: int = 110,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    accept_less_cpts: bool = False,
    return_scores: bool = True,
) -> Tuple[np_ndarray, np_ndarray, np_ndarray]: ...


def segment_signal(
    raw_signal: np_ndarray,
    num_events: int = 110,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    accept_less_cpts=False,
    return_scores=False,
) -> Union[Tuple[np_ndarray, np_ndarray], Tuple[np_ndarray, np_ndarray, np_ndarray]]:
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
    np.ndarray
        Array containing segmented, normalized signal values.
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

    dwell_times = valid_cpts[1:] - valid_cpts[:-1]
    event_means = compute_base_means(raw_signal, valid_cpts)

    if return_scores:
        return event_means, dwell_times, scores
    else:
        return event_means, dwell_times


@overload
def segment_signal_with_consensus_guided_barcode_refinement(
    raw_signal: np_ndarray,
    consensus_signal: np_ndarray,
    consensus_divergence_index: int = 83,
    num_events_signal: int = 110,
    num_events_barcode: int = 25,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    normalize_wrt_method: str = "mean",
    return_adapter_segmentation_results: bool = False,
) -> Tuple[np_ndarray, np_ndarray, int, int]: ...


@overload
def segment_signal_with_consensus_guided_barcode_refinement(
    raw_signal: np_ndarray,
    consensus_signal: np_ndarray,
    consensus_divergence_index: int = 83,
    num_events_signal: int = 110,
    num_events_barcode: int = 25,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    normalize_wrt_method: str = "mean",
    return_adapter_segmentation_results: bool = True,
) -> Tuple[np_ndarray, np_ndarray, np_ndarray, np_ndarray, int, int]: ...


def segment_signal_with_consensus_guided_barcode_refinement(
    raw_signal: np_ndarray,
    consensus_signal: np_ndarray,
    consensus_divergence_index: int = 83,
    num_events_signal: int = 110,
    num_events_barcode: int = 25,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    normalize_wrt_method: str = "mean",
    return_adapter_segmentation_results: bool = False,
) -> Union[
    Tuple[np_ndarray, np_ndarray, int, int],
    Tuple[np_ndarray, np_ndarray, np_ndarray, np_ndarray, int, int],
]:
    """
    Segment the given raw signal using the `c_windowed_t_test` function.

    If segmentation is not possible (e.g. because there is too little signal for the desired number of events under
    the running_stat_width and/or min_obs_per_base), an empty array is returned.

    Parameters:
    ----------
    raw_signal : np.ndarray
        The raw signal array to be segmented.

    consensus_signal : np.ndarray
        The consensus signal array to be used for determining the barcode start position in the segmented signal.
        Multiple consensus signals can be provided by passing a 2D array with shape (n_signals, signal_length).
        In this case, majority voting across the consensus signals is used to determine the barcode start position.

    consensus_divergence_index : int, optional
        The divergence index for the consensus signal.

    num_events_signal : int, optional
        The desired number of events in the signal.

    num_events_barcode : int, optional
        The desired number of events in the barcode.

    min_obs_per_base : int, optional
        Minimum observations per base.

    running_stat_width : int, optional
        Width for the running statistic used for segmentation.

    accept_less_cpts : bool, optional
        Whether to accept less cpts than num_events_signal.

    Returns:
    -------
    np.ndarray
        Array containing normalized, consensus-guided barcode-refined, segmented signal values (mean value per event).
    """

    def _get_barcode_start(
        segment_signal: np_ndarray,
        consensus_signal: np_ndarray,
        consensus_divergence_index: int,
    ) -> int:
        if consensus_signal.ndim == 1:
            consensus_signal = consensus_signal.reshape(1, -1)

        positions = []

        norm_segment_signal = normalize(segment_signal, "mean")
        for cons_signal in consensus_signal:
            path = np.array(
                warping_path_fast(
                    norm_segment_signal,
                    cons_signal,
                    include_distance=False,
                    window=30,
                    penalty=1.0,
                )
            )  # [(query_idx, ref_idx), query_idx, ref_idx), ...]

            path_index = np.argmax(path[:, 1] == consensus_divergence_index)
            positions.append(path[path_index, 0])

        pos = np.argmax(np.bincount(positions))
        return int(pos)

    adapter_event_means, adapter_dwell_times, adapter_scores = segment_signal(
        raw_signal,
        num_events=num_events_signal,
        min_obs_per_base=min_obs_per_base,
        running_stat_width=running_stat_width,
        accept_less_cpts=False,
        return_scores=True,
    )

    if adapter_event_means.size == 0:
        return np.array([])  # no refinement possible TODO: handle this case

    seg_barcode_start = _get_barcode_start(
        adapter_event_means, consensus_signal, consensus_divergence_index
    )
    sig_barcode_start = int(np.sum(adapter_dwell_times[:seg_barcode_start]))

    barcode_scores = adapter_scores[sig_barcode_start:]

    valid_cpts = discrepenacy_curve_to_cpts(
        barcode_scores,
        num_events=num_events_barcode,
        min_obs_per_base=min_obs_per_base,
        running_stat_width=running_stat_width,
        accept_less_cpts=False,
    )

    barcode_dwell_times = valid_cpts[1:] - valid_cpts[:-1]
    barcode_event_means = compute_base_means(raw_signal[sig_barcode_start:], valid_cpts)

    norm_barcode_event_means = normalize_wrt(
        barcode_event_means, adapter_event_means, normalize_wrt_method
    )

    if return_adapter_segmentation_results:
        return (
            norm_barcode_event_means,
            barcode_dwell_times,
            adapter_event_means,
            adapter_dwell_times,
            seg_barcode_start,
            sig_barcode_start,
        )
    else:
        return (
            norm_barcode_event_means,
            barcode_dwell_times,
            seg_barcode_start,
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
    consensus_signal: np.ndarray = np.array([]),
    consensus_divergence_index: int = 0,
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

    seg_barcode_start = None
    sig_barcode_start = None
    if spc.segmentation.consensus_refinement:
        if type(spc.segmentation.barcode_num_events) == int:
            raise ValueError(
                "barcode_num_events is an integer in consensus refinement mode, "
                "use a tuple instead"
            )

        (
            norm_segment_avgs,
            dwell_times,
            adapter_event_means,
            adapter_dwell_times,
            seg_barcode_start,
            sig_barcode_start,
        ) = segment_signal_with_consensus_guided_barcode_refinement(
            adapter_sig,
            consensus_signal,
            consensus_divergence_index=consensus_divergence_index,
            num_events_signal=spc.segmentation.num_events,
            num_events_barcode=spc.segmentation.barcode_num_events[0],
            min_obs_per_base=min(
                spc.segmentation.min_obs_per_base,
                round(adapter_sig.size / spc.segmentation.num_events / 2),
            ),
            running_stat_width=min(
                spc.segmentation.running_stat_width,
                round(adapter_sig.size / spc.segmentation.num_events),
            ),
            normalize_wrt_method=spc.segmentation.normalization,
            return_adapter_segmentation_results=True,
        )
        if not norm_segment_avgs.size:
            return ReadResult(
                success=False,
                fail_reason="event segmentation failed",
                barcode_fpt=np.array([]),
                dwell_times=np.array([]),
                detect_results=detect_results,
            )

        adapter_dt_med = float(np.median(adapter_dwell_times))
        adapter_dt_mad = float(np.median(np.abs(adapter_dwell_times - adapter_dt_med)))

        adapter_event_mean = float(adapter_event_means.mean())
        adapter_event_std = float(adapter_event_means.std())
        adapter_event_med = float(np.median(adapter_event_means))
        adapter_event_mad = float(
            np.median(np.abs(adapter_event_means - adapter_event_med))
        )

    else:

        segment_avgs, dwell_times = segment_signal(
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
        else spc.segmentation.barcode_num_events[1]
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
        seg_barcode_start=seg_barcode_start,
        sig_barcode_start=sig_barcode_start,
    )
