import importlib.resources as pkg_resources
import os
import traceback
from functools import partial
from typing import Any, Callable, Dict, List, Tuple, Union

import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm

from adapted.detect.combined import DetectResults
from adapted.output import save_detected_boundaries, save_traces

from warpdemux.config.classification import ClassificationConfig
from warpdemux.config.config import Config
from warpdemux.config.sig_proc import SigProcConfig
from warpdemux.models import model_files
from warpdemux.models.dtw_svm import DTW_SVM_Model
from warpdemux.sig_proc import ReadResult, barcode_fpt, detect_results_to_fpt


def resegment_wrapper(
    signal: np.ndarray,
    signal_len: int,
    full_sig_len: int,  # required for matching adapted.file_proc.file_proc.process_in_batches signature
    read_id: str,
    spc: "SigProcConfig",
    detect_result_dict: Dict[str, Dict[str, Any]],
) -> "ReadResult":
    try:
        _vals = detect_result_dict.get(read_id, {})
        _vals.pop("adapter_med_dt", None)
        _vals.pop("adapter_mad_dt", None)
        _vals.update(
            {
                "success": True,
            }
        )
        detect_results = DetectResults(**_vals)

        proc_res = detect_results_to_fpt(
            calibrated_signal=signal[:signal_len],
            spc=spc,
            detect_results=detect_results,
        )
        proc_res.set_read_id(read_id)
        return proc_res
    except Exception as e:
        # shouldn't happen, errors should be caught earlier and result in empty fpt
        print(f"Failed on read {read_id}: {e}")
        traceback.print_exc()
        return ReadResult(read_id=read_id, success=False, fail_reason="unknown")


def barcode_fpt_wrapper(
    signal: np.ndarray,
    signal_len: int,
    full_sig_len: int,
    read_id: str,
    spc: "SigProcConfig",
    llr_return_trace: bool = False,
) -> "ReadResult":
    try:
        # partial ReadResult, missing read_id
        proc_res = barcode_fpt(
            calibrated_signal=signal[
                :signal_len
            ],  # in case signal_len < preloaded buffer
            full_sig_len=full_sig_len,
            spc=spc,
            llr_return_trace=llr_return_trace,
        )
        proc_res.set_read_id(read_id)
        return proc_res

    except Exception as e:
        # shouldn't happen, errors should be caught earlier and result in empty fpt
        print(f"Failed on read {read_id}: {e}")
        traceback.print_exc()
        return ReadResult(read_id=read_id, success=False, fail_reason="unknown")


def save_processed_signals(
    list_of_processing_results: List["ReadResult"],
    filename: str,
    save_dwell_time: bool = True,
):
    """Save processed reads to a numpy file."""

    read_ids, barcode_fpts, dwell_times = [], [], []
    for res in list_of_processing_results:
        read_ids.append(res.read_id)
        barcode_fpts.append(res.barcode_fpt)
        dwell_times.append(res.dwell_times)
    read_ids = np.array(read_ids)
    barcode_fpts = np.array(barcode_fpts)
    dwell_times = np.array(dwell_times)

    if save_dwell_time:
        np.savez(
            filename, read_ids=read_ids, signals=barcode_fpts, dwell_times=dwell_times
        )
    else:
        np.savez(filename, read_ids=read_ids, signals=barcode_fpts)

    return read_ids, barcode_fpts, dwell_times


def save_predictions(
    predictions: pd.DataFrame,
    read_ids: Union[List[str], np.ndarray],
    output_dir: str,
    suffix: str = "",
) -> None:
    os.makedirs(output_dir, exist_ok=True)

    cols = predictions.columns.tolist()
    predictions["read_id"] = read_ids
    predictions[["read_id", *cols]].to_csv(
        os.path.join(output_dir, f"barcode_predictions{suffix}.csv"),
        index=False,
    )


def load_model(model_name: str) -> DTW_SVM_Model:
    with pkg_resources.path(model_files, f"{model_name}.joblib") as model_path:
        return joblib.load(model_path)


def predict_batch(
    model: DTW_SVM_Model,
    signals: np.ndarray,
    read_ids: np.ndarray,
    batch_idx: int,
    output_dir: str,
    classif_config: ClassificationConfig,
) -> None:
    pbar = tqdm(
        total=1, desc=f"Predicting on fpts batch {batch_idx}", position=1, leave=False
    )
    predictions = model.predict(
        signals,
        pbar=True,  # computing kernel matrix pbar
        pbar_kwargs=dict(position=2, leave=False),
        block_size=classif_config.block_size,
        nproc=classif_config.num_proc,
        return_df=True,
    )
    pbar.update(1)

    save_predictions(
        predictions=predictions,
        read_ids=read_ids,
        output_dir=output_dir,
        suffix=f"_{batch_idx}",
    )
    pbar.close()


def _save_detect_fpts_batch(
    pass_or_fail: str,
    results: List[ReadResult],
    batch_idx: int,
    output_dir: str,
    save_llr_trace: bool = False,
    save_dwell_time: bool = False,
    save_boundaries: bool = True,
) -> Union[Tuple[np.ndarray, np.ndarray], None]:
    """Save detected boundaries and return read_ids and signals if pass_or_fail == 'pass'."""

    if pass_or_fail == "pass":
        fn = "detected_boundaries"
    elif pass_or_fail == "fail":
        fn = "failed_reads"
    else:
        raise ValueError(
            f"Invalid pass_or_fail: {pass_or_fail}. Must be 'pass' or 'fail'."
        )

    if save_boundaries:
        save_detected_boundaries(
            results,  # type: ignore
            os.path.join(output_dir, f"{fn}_{batch_idx}.csv"),
            save_fail_reasons=pass_or_fail == "fail",
        )

    if save_llr_trace:
        save_traces(
            results,  # type: ignore
            os.path.join(output_dir, f"{fn}_{batch_idx}_llr_trace.npz"),
        )

    if pass_or_fail == "pass":
        read_ids, signals, _ = save_processed_signals(
            results,
            os.path.join(output_dir, f"barcode_fpts_{batch_idx}.npz"),
            save_dwell_time=save_dwell_time,
        )

        return read_ids, signals

    return None


def save_detect_fpts_batch(
    pass_or_fail: str,
    results: List[ReadResult],
    batch_idx: int,
    output_dir: str,
    save_llr_trace: bool = False,
    save_dwell_time: bool = False,
    save_boundaries: bool = True,
) -> Union[Tuple[np.ndarray, np.ndarray], None]:
    """Save detected boundaries and return read_ids and signals if pass_or_fail == 'pass'."""
    return _save_detect_fpts_batch(
        pass_or_fail=pass_or_fail,
        results=results,
        batch_idx=batch_idx,
        output_dir=output_dir,
        save_llr_trace=save_llr_trace,
        save_dwell_time=save_dwell_time,
        save_boundaries=save_boundaries,
    )


def save_detect_fpts_and_predict_batch(
    pass_or_fail: str,
    results: List[ReadResult],
    batch_idx: int,
    output_dir: str,
    model: DTW_SVM_Model,
    classif_config: ClassificationConfig,
    save_llr_trace: bool = False,
    save_dwell_time: bool = False,
    save_boundaries: bool = True,
) -> None:
    res = _save_detect_fpts_batch(
        pass_or_fail=pass_or_fail,
        results=results,
        batch_idx=batch_idx,
        output_dir=output_dir,
        save_llr_trace=save_llr_trace,
        save_dwell_time=save_dwell_time,
        save_boundaries=save_boundaries,
    )

    if res is not None:  # 'pass' batch
        read_ids, signals = res
        predict_batch(
            model=model,
            signals=signals,
            read_ids=read_ids,
            batch_idx=batch_idx,
            output_dir=output_dir,
            classif_config=classif_config,
        )


def get_save_detect_fpts_and_predict_batch_fn(
    config: Config,
) -> Callable:
    print("Loading model...")

    assert (
        config.classif is not None
    ), "Classification config is required for classification."
    model = load_model(config.classif.model_name)

    return partial(
        save_detect_fpts_and_predict_batch, model=model, classif_config=config.classif
    )
