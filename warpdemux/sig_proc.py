"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

from dataclasses import dataclass
from typing import Any, Dict, Optional
import logging
import numpy as np
from scipy.signal import find_peaks

from adapted.detect.combined import DetectResults, combined_detect
from adapted.file_proc.file_proc import (
    ReadResult as AdaptedReadResult,
)

from warpdemux.config.sig_proc import SigProcConfig
from warpdemux.segmentation.segmentation import (
    compute_base_means,
    identify_stalls,
    remove_stall_cpts,
    windowed_t_test,
)


@dataclass
class ReadResult(AdaptedReadResult):
    # read_id
    # success
    # fail_reason
    # detect_results
    barcode_fpt: Optional[np.ndarray] = None
    dwell_times: Optional[np.ndarray] = None
    adapter_med_dt: Optional[float] = None
    adapter_mad_dt: Optional[float] = None

    def to_summary_dict(self) -> Dict[str, Any]:
        summary_dict = super().to_summary_dict()
        summary_dict.update(
            {
                "adapter_med_dt": self.adapter_med_dt,
                "adapter_mad_dt": self.adapter_mad_dt,
            }
        )
        return summary_dict

    def set_read_id(self, read_id: str):
        self.read_id = read_id


###############
# Normalization
###############


def impute_window_median(
    signal: np.ndarray, indices: np.ndarray, window_size: int = 5
) -> np.ndarray:
    """Impute the given indices of the signal with window medians.

    Parameters:
    -----------
    signal : np.ndarray
        The input signal array.

    indices : np.ndarray
        Indices of the signal array to be imputed.

    window_size : int, optional
        The size of the window used for calculating the median for imputation.

    Returns:
    --------
    np.ndarray
        The signal array with imputed values.
    """
    # check if window size is at least 3
    if window_size < 3:
        msg = "window_size should be at least 3"
        logging.error(msg)
        raise ValueError(msg)

    half_window = window_size // 2

    signal_copy = signal.copy()
    for i_spike in indices:
        if i_spike < 2:
            signal_copy[i_spike] = np.median(signal_copy[: i_spike + window_size])
        elif i_spike > (len(signal_copy) - (half_window + 1)):
            signal_copy[i_spike] = np.median(signal_copy[i_spike - window_size :])
        else:
            signal_copy[i_spike] = np.median(
                signal_copy[i_spike - half_window : i_spike + half_window]
            )
    return signal_copy


def mad_outlier_indices(signal: np.ndarray, outlier_thresh: float = 5) -> np.ndarray:
    """
    Identify the indices of values in the signal that are outliers based on MAD.

    Parameters:
    ----------
    signal : np.ndarray
        The input signal array.

    outlier_thresh : float, optional
        The threshold multiplier for MAD to determine outliers. Default is 5.

    Returns:
    -------
    np.ndarray
        Array of indices that correspond to outlier values in the signal.
    """

    med = np.median(signal)
    mad = np.median(np.abs(signal - med))
    lower_lim = med - (mad * outlier_thresh)
    upper_lim = med + (mad * outlier_thresh)

    return ((signal < lower_lim) | (signal > upper_lim)).nonzero()[0]


def mad_winsor(
    signal: np.ndarray, outlier_thresh: float = 5, window_size: int = 5
) -> np.ndarray:
    """
    Apply winsorization to a signal using MAD to limit the effect of outliers.

    Parameters
    ----------
    signal : np.ndarray
        The input signal array to be winsorized.
    outlier_thresh : float, optional
        The threshold multiplier for MAD to determine outliers. Default is 5.
    window_size : int, optional
        The size of the window used for imputation of outliers. Default is 5.

    Returns
    -------
    np.ndarray
        The winsorized signal array with outliers handled according to MAD.
    """

    outlier_indices = mad_outlier_indices(signal, outlier_thresh)
    wins_signal = impute_window_median(signal, outlier_indices, window_size)
    return wins_signal


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

    if np.isnan(signal).any():
        if not accept_nan:
            msg = "Signal contains NaN values."
            logging.error(msg)
            raise ValueError(msg)
        else:
            shift = np.nanmedian(signal)
            scale = np.nanmedian(np.abs(signal - shift))
    else:
        shift = np.median(signal)
        scale = np.median(np.abs(signal - shift))

    norm_signal = (signal - shift) / scale
    return norm_signal


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
        if accept_nan:
            norm_signal = (signal - np.nanmean(signal)) / np.nanstd(signal)
        else:
            norm_signal = (signal - np.mean(signal)) / np.std(signal)
    elif method == "none":
        norm_signal = signal
    elif method == "median":
        norm_signal = mad_normalize(signal, accept_nan=accept_nan)
    else:
        msg = f"Normalization method {method} not recognized."
        logging.error(msg)
        raise ValueError(msg)

    return norm_signal


###############
# Adapter
###############


def segment_signal(
    raw_signal: np.ndarray,
    num_events: int = 110,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
    remove_stalls: bool = False,
    stall_window_size=5 * 25,
    stall_threshold=2,
    stall_edge_buffer=100,
    stall_min_consecutive_obs=200,
    stall_n_windows=5,
    stall_mini_window_size=25,
    accept_less_cpts=False,
):
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

    remove_stalls : bool, optional
        Whether to remove stall points from the valid points.

    stall_window_size : int, optional
        The size of the window used for imputation of outliers.

    stall_threshold : int, optional
        The threshold for stall detection.

    stall_edge_buffer : int, optional
        The buffer to add to the edge of the window for stall detection.

    stall_min_consecutive_obs : int, optional
        The minimum number of consecutive observations to consider a stall.

    stall_n_windows : int, optional
        The number of windows to use for stall detection.

    stall_mini_window_size : int, optional
        The size of the mini window used for stall detection.

    Returns:
    -------
    np.ndarray
        Array containing segmented, normalized signal values.
    """

    scores = windowed_t_test(
        raw_signal,
        min_obs_per_base=min_obs_per_base,
        running_stat_width=running_stat_width,
    )
    peaks, _ = find_peaks(scores, distance=min_obs_per_base)

    if peaks.size < num_events and not accept_less_cpts:
        return np.array([])

    valid_cpts = peaks[np.argsort(scores[peaks])[-num_events:]] + running_stat_width
    valid_cpts.sort()

    if remove_stalls:
        stall_cpts = identify_stalls(
            raw_signal,
            window_size=stall_window_size,
            threshold=stall_threshold,
            edge_buffer=stall_edge_buffer,
            min_consecutive_obs=stall_min_consecutive_obs,
            n_windows=stall_n_windows,
            mini_window_size=stall_mini_window_size,
        )
        valid_cpts = remove_stall_cpts(stall_cpts, valid_cpts)

    if valid_cpts[0] != 0:
        valid_cpts = np.insert(valid_cpts, 0, 0)
    if valid_cpts[-1] != raw_signal.size:
        valid_cpts = np.append(valid_cpts, raw_signal.size)

    dwell_times = valid_cpts[1:] - valid_cpts[:-1]
    event_means = compute_base_means(raw_signal, valid_cpts)

    return event_means, dwell_times


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

    # impute mad outliers, required: signal flickers strongly affect segmentation results
    norm_adapter_sig = mad_winsor(
        adapter_sig,
        outlier_thresh=spc.sig_norm.outlier_thresh,
        window_size=spc.sig_norm.winsor_window,
    )

    try:
        norm_adapter_sig = normalize(
            norm_adapter_sig,
            spc.sig_extract.normalization,
            accept_nan=False,
        )  # can be used to normalize by adapter mean or median prior to segmentation
    except Exception as e:
        return ReadResult(
            success=False,
            fail_reason=f"signal normalization failed: {e}",
            barcode_fpt=np.array([]),
            dwell_times=np.array([]),
            detect_results=detect_results,
        )

    segment_avgs, dwell_times = segment_signal(
        norm_adapter_sig,
        num_events=spc.segmentation.num_events,
        min_obs_per_base=min(
            spc.segmentation.min_obs_per_base,
            round(norm_adapter_sig.size / spc.segmentation.num_events / 2),
        ),
        running_stat_width=min(
            spc.segmentation.running_stat_width,
            round(norm_adapter_sig.size / spc.segmentation.num_events),
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

    if (
        norm_segment_avgs.size < spc.segmentation.barcode_num_events
        and spc.segmentation.pad_if_less
    ):
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

    adapter_med_dt = float(np.median(dwell_times))
    adapter_mad_dt = float(np.median(np.abs(dwell_times - adapter_med_dt)))

    return ReadResult(
        success=True,
        fail_reason="",
        barcode_fpt=norm_segment_avgs[
            -min(spc.segmentation.barcode_num_events, norm_segment_avgs.size) :
        ],
        dwell_times=dwell_times[
            -min(spc.segmentation.barcode_num_events, dwell_times.size) :
        ],
        detect_results=detect_results,
        adapter_med_dt=adapter_med_dt,
        adapter_mad_dt=adapter_mad_dt,
    )


def barcode_fpt(
    calibrated_signal: np.ndarray,
    full_sig_len: int,
    spc: "SigProcConfig",
    llr_return_trace: bool = False,
) -> "ReadResult":
    """Process the raw signal, extract the adapter segment, preprocess, and segment it.

    mvs_median_after_window = 5000: since we are not limited by the speed of the live balancing,
    we can afford to use a larger window for the mvs polyA detection.
    Even if we have a very short polyA tail with a signal drop after, increasing the median window will
    still detect the boundary as the median will be pulled up by the subsequent RNA signal.

    Default parameters were fitted for non-live data."""

    detect_results = combined_detect(
        calibrated_signal,
        full_signal_len=full_sig_len,
        spc=spc,
        llr_return_trace=llr_return_trace,
    )

    return detect_results_to_fpt(calibrated_signal, spc, detect_results)
