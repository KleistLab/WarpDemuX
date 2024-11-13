"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np
from adapted.detect.combined import DetectResults
from adapted.file_proc.file_proc import ReadResult as AdaptedReadResult
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


###############
# Adapter
###############


def segment_signal(
    raw_signal: np.ndarray,
    num_events: int = 110,
    min_obs_per_base: int = 15,
    running_stat_width: int = 30,
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

    Returns:
    -------
    np.ndarray
        Array containing segmented, normalized signal values.
    """

    scores = windowed_t_test(
        raw_signal,
        running_stat_width=running_stat_width,
    )
    peaks, _ = find_peaks(scores, distance=min_obs_per_base)

    if peaks.size < num_events and not accept_less_cpts:
        return np.array([])

    valid_cpts = peaks[np.argsort(scores[peaks])[-num_events:]] + running_stat_width
    valid_cpts.sort()

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

    if norm_segment_avgs.size < spc.segmentation.barcode_num_events:
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
