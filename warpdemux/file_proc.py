"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import concurrent.futures
import importlib.resources as pkg_resources
import inspect
import logging
import multiprocessing
import os
import queue
import re
import subprocess
import sys
import threading
import time
import traceback
from typing import Any, Callable, Dict, Generator, List, Set, Tuple, Union

import joblib
import numpy as np
import pandas as pd
from adapted.container_types import DetectResults
from adapted.detect.cnn import load_cnn_model
from adapted.detect.combined import (combined_detect_cnn, combined_detect_llr2,
                                     combined_detect_start_peak)
from adapted.output import save_detected_boundaries
from pod5.reader import Reader
from tqdm import tqdm

from warpdemux._consensus import ALL as CONSENSUS_ALL
from warpdemux.config.config import Config
from warpdemux.config.sig_proc import SigProcConfig
from warpdemux.models import model_files
from warpdemux.models.dtw_svm import DTW_SVM_Model
from warpdemux.sig_proc import ReadResult, detect_results_to_fpt

_STOP_SIGNAL = threading.Event()

# Initialize the shared batch index and a threading lock
_BIDX_LOCK_PASS = threading.Lock()
_BIDX_LOCK_FAIL = threading.Lock()
_BIDX_LOCK_PREDICT = threading.Lock()

_BIDX_PASS = 0
_BIDX_FAIL = 0
_BIDX_PREDICT = 0

_TOTAL_READS_LOCK = threading.Lock()
_TOTAL_READS = -1


def set_global_var(name, value):
    if name == "total_reads":
        with _TOTAL_READS_LOCK:
            global _TOTAL_READS
            _TOTAL_READS = value
    elif name == "bidx_pass":
        with _BIDX_LOCK_PASS:
            global _BIDX_PASS
            _BIDX_PASS = value
    elif name == "bidx_fail":
        with _BIDX_LOCK_FAIL:
            global _BIDX_FAIL
            _BIDX_FAIL = value
    elif name == "bidx_predict":
        with _BIDX_LOCK_PREDICT:
            global _BIDX_PREDICT
            _BIDX_PREDICT = value
    else:
        raise ValueError(f"Invalid name: {name}")


def get_global_var(name: str) -> int:
    if name == "total_reads":
        with _TOTAL_READS_LOCK:
            return _TOTAL_READS
    elif name == "bidx_pass":
        with _BIDX_LOCK_PASS:
            return _BIDX_PASS
    elif name == "bidx_fail":
        with _BIDX_LOCK_FAIL:
            return _BIDX_FAIL
    elif name == "bidx_predict":
        with _BIDX_LOCK_PREDICT:
            return _BIDX_PREDICT
    else:
        raise ValueError(f"Invalid name: {name}")


def increment_ridx(ridx_dict, ridx_name, ridx_lock, value=1):
    with ridx_lock:
        ridx_dict[ridx_name] += value


def increment_bidx(bidx_dict, bidx_name, bidx_lock, value=1):
    with bidx_lock:
        bidx_dict[bidx_name] += value


def bidx_name_from_ridx_name(ridx_name: str) -> str:
    if "pass" in ridx_name:
        return "pass"
    elif "fail" in ridx_name:
        return "fail"
    elif "predict" in ridx_name:
        return "predict"
    else:
        raise ValueError(f"Invalid ridx_name: {ridx_name}")


def determine_bidx_from_file(file: str) -> int:
    return int(file.split("_")[-1].split(".")[0])


def determine_max_bidx_from_dir(dir_path: str) -> int:

    return max(determine_bidx_from_file(f) for f in os.listdir(dir_path))


def scan_processed_reads(
    continue_from_path: str, scan_failed: bool = False, result_type: str = "predictions"
) -> Tuple[Set[str], int, int]:
    processed_reads = set()
    max_pass_bidx = -1
    max_fail_bidx = -1

    if result_type not in ["predictions", "fingerprints"]:
        raise ValueError(
            f"Invalid result_type: {result_type}. Must be 'predictions' or 'fingerprints'."
        )

    if scan_failed:
        fail_sub = os.path.join(continue_from_path, "failed_reads")
        for file in os.listdir(fail_sub):
            if file.startswith("failed_reads_") and file.endswith(".csv"):
                bidx = int(file.split("_")[-1].split(".")[0])
            max_fail_bidx = max(max_fail_bidx, bidx)
            with open(os.path.join(fail_sub, file), "r") as f:
                processed_reads.update(line.split(",")[0] for line in f.readlines()[1:])

    if result_type == "predictions":
        pass_sub = os.path.join(continue_from_path, "predictions")
        starts_with = "barcode_predictions_"
        extension = "csv"
    elif result_type == "fingerprints":
        pass_sub = os.path.join(continue_from_path, "fingerprints")
        starts_with = "barcode_fpts_"
        extension = "npz"

    for file in os.listdir(pass_sub):
        if file.startswith(starts_with) and file.endswith(extension):
            bidx = int(file.split("_")[-1].split(".")[0])
            max_pass_bidx = max(max_pass_bidx, bidx)
            if extension == "csv":
                with open(os.path.join(pass_sub, file), "r") as f:
                    processed_reads.update(
                        line.split(",")[0] for line in f.readlines()[1:]
                    )
            elif extension == "npz":
                npz = np.load(os.path.join(pass_sub, file))
                processed_reads.update(npz["read_ids"])
    return processed_reads, max_pass_bidx, max_fail_bidx


def handle_previous_results(
    config: Config,
) -> Set[str]:

    scan_failed = config.task.preprocess
    result_type = "predictions" if config.task.predict else "fingerprints"

    processed_reads, max_pass_bidx, max_fail_bidx = scan_processed_reads(
        config.input.continue_from, scan_failed, result_type
    )
    config.batch.bidx_pass = max_pass_bidx + 1
    config.batch.bidx_fail = max_fail_bidx + 1
    return processed_reads


def barcode_fpt_wrapper(
    signal: np.ndarray,
    read_id: str,
    detect_results: DetectResults,
    spc: "SigProcConfig",
) -> "ReadResult":
    """Note: signal should not contain nans"""

    if spc.segmentation.consensus_refinement:
        if spc.segmentation.consensus_model == "":
            raise ValueError(
                "consensus_model must be specified when consensus_refinement is True"
            )
        if spc.segmentation.consensus_model not in CONSENSUS_ALL:
            raise ValueError(
                f"Invalid consensus model: {spc.segmentation.consensus_model}"
            )
        consensus_query = CONSENSUS_ALL[spc.segmentation.consensus_model]
    else:
        consensus_query = np.array([])

    try:
        # partial ReadResult, missing read_id
        proc_res = detect_results_to_fpt(
            signal,
            spc,
            detect_results,
            consensus_query,
        )
        proc_res.set_read_id(read_id)
        return proc_res

    except Exception as e:
        # shouldn't happen, errors should be caught earlier and result in empty fpt
        logging.debug(f"Failed on read {read_id}: {e}")
        logging.debug(traceback.format_exc())
        return ReadResult(read_id=read_id, success=False, fail_reason="unknown")


def yield_signals_from_pod5(
    pod5_files: Union[List[str], Set[str]],
    read_ids_incl: Set[str],
    read_ids_excl: Set[str],
    batch_size: int,
    preload_size: int,
) -> Generator[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], None, None]:
    if read_ids_incl and read_ids_excl:
        read_ids_incl = read_ids_incl.difference(read_ids_excl)
        read_ids_excl = set()

    selection = list(read_ids_incl) if read_ids_incl else None

    N = batch_size
    m = preload_size
    i_id = 0

    signals = np.empty((N, m), dtype=np.float32)
    full_lengths = np.empty(N, dtype=np.int32)
    in_arr_lengths = np.empty(N, dtype=np.int32)
    read_ids = np.empty(N, dtype=object)

    for filename in pod5_files:
        with Reader(filename) as file_obj:
            for read_record in file_obj.reads(selection=selection, missing_ok=True):
                if str(read_record.read_id) in read_ids_excl:
                    continue

                _m = min(m, read_record.num_samples)
                full_lengths[i_id] = read_record.num_samples
                in_arr_lengths[i_id] = _m
                signals[i_id, :_m] = read_record.signal_pa[:_m]
                if _m < m:
                    signals[i_id, _m:] = np.nan
                read_ids[i_id] = str(read_record.read_id)

                if i_id == N - 1:
                    yield signals, in_arr_lengths, full_lengths, read_ids
                    signals = np.empty(
                        (N, m),
                        dtype=np.float32,
                    )
                    full_lengths = np.empty(N, dtype=np.int32)
                    in_arr_lengths = np.empty(N, dtype=np.int32)
                    read_ids = np.empty(N, dtype=object)
                    i_id = 0
                else:
                    i_id += 1

    if i_id > 0:
        yield signals[:i_id], in_arr_lengths[:i_id], full_lengths[:i_id], read_ids[
            :i_id
        ]


def yield_fpts_from_npz(
    npz_files: Union[List[str], Set[str]],
    read_ids_incl: Set[str],
    read_ids_excl: Set[str],
    batch_size: int,
) -> Generator[Tuple[np.ndarray, np.ndarray], None, None]:
    if read_ids_incl and read_ids_excl:
        read_ids_incl = read_ids_incl.difference(read_ids_excl)
        read_ids_excl = set()

    N = batch_size

    fpts = np.empty((0, 0), dtype=np.float32)
    read_ids = np.empty(0, dtype=object)

    for filename in npz_files:
        with np.load(filename) as data:
            file_fpts = data["signals"]
            file_read_ids = data["read_ids"]
            # Get indices of excluded read_ids
            if read_ids_excl:
                exclude_idx = np.array(
                    [i for i, rid in enumerate(file_read_ids) if rid in read_ids_excl]
                )
                file_fpts = np.delete(file_fpts, exclude_idx, axis=0)
                file_read_ids = np.delete(file_read_ids, exclude_idx)
            elif read_ids_incl:
                include_idx = np.array(
                    [i for i, rid in enumerate(file_read_ids) if rid in read_ids_incl]
                )
                file_fpts = file_fpts[include_idx]
                file_read_ids = file_read_ids[include_idx]

            if fpts.size > 0:
                fpts = np.concatenate((fpts, file_fpts), axis=0)
                read_ids = np.concatenate((read_ids, file_read_ids), axis=0)
            else:
                fpts = file_fpts
                read_ids = file_read_ids

            while len(fpts) >= N:
                batch_fpts = fpts[:N].copy()
                batch_read_ids = read_ids[:N].copy()
                yield batch_fpts, batch_read_ids
                fpts = fpts[N:]
                read_ids = read_ids[N:]

    if fpts.size > 0:
        yield fpts, read_ids


def worker_enqueue_minibatches_raw(
    files: Set[str],
    read_ids_incl: Set[str],
    read_ids_excl: Set[str],
    config: Config,
    preloaded_minibatch_queue: multiprocessing.Queue,
    ridx_dict: dict,
    ridx_lock: multiprocessing.Lock,  # type: ignore
):

    for minibatch in yield_signals_from_pod5(
        pod5_files=files,
        read_ids_incl=read_ids_incl,
        read_ids_excl=read_ids_excl,
        batch_size=config.batch.minibatch_size,
        preload_size=config.sig_proc.sig_preload_size,
    ):
        preloaded_minibatch_queue.put(minibatch)  # has maxsize
        increment_ridx(ridx_dict, "enqueued", ridx_lock, len(minibatch[0]))

    preloaded_minibatch_queue.put(None)  # signal end of work
    logging.debug("Worker enqueue_minibatches finished")


def worker_enqueue_minibatches_fpts(
    npz_files: Set[str],
    read_ids_incl: Set[str],
    read_ids_excl: Set[str],
    config: Config,
    preloaded_minibatch_queue: multiprocessing.Queue,
    ridx_dict: dict,
    ridx_lock: multiprocessing.Lock,  # type: ignore
):

    for minibatch in yield_fpts_from_npz(
        npz_files=npz_files,
        read_ids_incl=read_ids_incl,
        read_ids_excl=read_ids_excl,
        batch_size=config.batch.minibatch_size,
    ):
        preloaded_minibatch_queue.put(minibatch)  # has maxsize
        increment_ridx(ridx_dict, "enqueued", ridx_lock, len(minibatch[0]))

    preloaded_minibatch_queue.put(None)  # signal end of work
    logging.debug("Worker enqueue_minibatches finished")


def worker_detect_and_predict_on_preloaded_signals(
    preloaded_minibatch,
    model_predict,
    model_detect,
    config,
    ridx_dict,
    ridx_lock,
    save_pass_queue,
    save_fail_queue,
    save_predict_queue,
):
    success = []
    fail = []
    signals, _, full_lengths, read_ids = preloaded_minibatch

    if config.sig_proc.primary_method == "llr":
        detect_results = combined_detect_llr2(
            batch_of_signals=signals,
            full_signal_lens=full_lengths,
            spc=config.sig_proc,
        )
    elif config.sig_proc.primary_method == "cnn":
        detect_results = combined_detect_cnn(
            batch_of_signals=signals,
            full_signal_lens=full_lengths,
            model=model_detect,
            spc=config.sig_proc,
        )
    elif config.sig_proc.primary_method == "start_peak":
        detect_results = combined_detect_start_peak(
            batch_of_signals=signals,
            full_signal_lens=full_lengths,
            spc=config.sig_proc,
        )

    del model_detect
    assert isinstance(detect_results, list)

    for signal, read_id, detect_result in zip(signals, read_ids, detect_results):
        result = barcode_fpt_wrapper(
            signal=signal,
            read_id=read_id,
            detect_results=detect_result,
            spc=config.sig_proc,
        )
        if result.success:
            success.append(result)
        else:
            fail.append(result)

    del preloaded_minibatch, signals, full_lengths, read_ids

    if fail:
        save_fail_queue.put(fail)
        increment_ridx(ridx_dict, "done_fail", ridx_lock, len(fail))

    if success:
        save_output_pass = config.output.save_boundaries or config.output.save_fpts
        if save_output_pass:
            save_pass_queue.put(success)
        increment_ridx(ridx_dict, "done_pass", ridx_lock, len(success))

        if config.task.predict:
            read_ids, fpts = zip(*[(res.read_id, res.barcode_fpt) for res in success])

            predictions = model_predict.predict(
                np.vstack(fpts),
                pbar=False,
                nproc=1,
                return_df=True,
            )
            predictions = add_read_id_col_to_predictions(predictions, read_ids)

            save_predict_queue.put(predictions)
            increment_ridx(ridx_dict, "done_predict", ridx_lock, len(predictions))


def worker_predict_on_preloaded_fpts(
    preloaded_minibatch,
    model_predict,
    model_detect,  # unused
    config,  # unused
    ridx_dict,
    ridx_lock,
    save_pass_queue,  # unused
    save_fail_queue,  # unused
    save_predict_queue,
):
    """Process preloaded fingerprints through prediction model.

    Args:
        preloaded_minibatch: Tuple of (fingerprints array, read_ids list)
        model_predict: Prediction model to use
        model_detect: Unused, for compatibility with worker_detect_predict_on_preloaded_signals
        config: Unused, for compatibility with worker_detect_predict_on_preloaded_signals
        ridx_dict: Shared dictionary for read indices
        ridx_lock: Lock for ridx_dict access
        save_pass_queue: Unused, for compatibility with worker_detect_predict_on_preloaded_signals
        save_fail_queue: Unused, for compatibility with worker_detect_predict_on_preloaded_signals
        save_predict_queue: Queue to save prediction results to
    """
    assert model_detect is None, "model_detect should be None"
    assert save_pass_queue is None, "save_pass_queue should be None"
    assert save_fail_queue is None, "save_fail_queue should be None"
    assert config is not None, "config should be provided"

    fpts, read_ids = preloaded_minibatch

    predictions = model_predict.predict(
        fpts,
        pbar=False,
        nproc=1,
        return_df=True,
    )
    predictions = add_read_id_col_to_predictions(predictions, read_ids)

    save_predict_queue.put(predictions)
    increment_ridx(ridx_dict, "done_predict", ridx_lock, len(predictions))


def _queue_batch_processor_df(
    queue, process_fn, ridx_dict, ridx_name, ridx_lock, bidx_dict, bidx_lock, config
):
    """Collects pd.DataFrame results of variable length from the queue and processes them in batches of batch_size rows."""
    batch: List[pd.DataFrame] = []
    batch_read_count = 0

    while True:
        # Wait for results from the queue
        result = queue.get()
        if result is None:
            break
        if not isinstance(result, pd.DataFrame):
            raise ValueError(f"Invalid result type: {type(result)}")
        batch.append(result)
        batch_read_count += len(result)

        # If the batch reaches the desired size, process it
        if batch_read_count >= config.batch.batch_size_output:
            current_batch = pd.concat(batch, axis=0)
            process_fn(
                current_batch.iloc[: config.batch.batch_size_output].copy(),
                config=config,
                bidx_dict=bidx_dict,
                bidx_lock=bidx_lock,
            )
            batch = [current_batch.iloc[config.batch.batch_size_output :]]
            batch_read_count -= config.batch.batch_size_output
            increment_ridx(
                ridx_dict, ridx_name, ridx_lock, config.batch.batch_size_output
            )
            increment_bidx(bidx_dict, bidx_name_from_ridx_name(ridx_name), bidx_lock, 1)

    # save remainder
    if batch_read_count > 0:
        current_batch = pd.concat(batch, axis=0)
        process_fn(
            current_batch, config=config, bidx_dict=bidx_dict, bidx_lock=bidx_lock
        )
        increment_ridx(ridx_dict, ridx_name, ridx_lock, batch_read_count)
        increment_bidx(bidx_dict, bidx_name_from_ridx_name(ridx_name), bidx_lock, 1)


def _queue_batch_processor_list(
    queue, process_fn, ridx_dict, ridx_name, ridx_lock, bidx_dict, bidx_lock, config
):
    """Collects ReadResult results from the queue and processes them in batches."""
    batch: List[ReadResult] = []
    batch_read_count = 0

    while True:
        # Wait for results from the queue
        result = queue.get()

        if result is None:
            break

        if not isinstance(result, list):
            raise ValueError(f"Invalid result type: {type(result)}")

        batch.extend(result)
        batch_read_count += len(result)

        if batch_read_count >= config.batch.batch_size_output:
            process_fn(
                batch[: config.batch.batch_size_output],
                config=config,
                bidx_dict=bidx_dict,
                bidx_lock=bidx_lock,
            )
            batch = batch[config.batch.batch_size_output :]
            batch_read_count -= config.batch.batch_size_output
            increment_ridx(
                ridx_dict, ridx_name, ridx_lock, config.batch.batch_size_output
            )
            increment_bidx(bidx_dict, bidx_name_from_ridx_name(ridx_name), bidx_lock, 1)

    # save remainder
    if batch_read_count > 0:
        process_fn(batch, config=config, bidx_dict=bidx_dict, bidx_lock=bidx_lock)
        increment_ridx(ridx_dict, ridx_name, ridx_lock, batch_read_count)
        increment_bidx(bidx_dict, bidx_name_from_ridx_name(ridx_name), bidx_lock, 1)


def queue_batch_processor(
    queue: multiprocessing.Queue,
    process_fn: Callable,
    ridx_dict: Dict[str, int],
    ridx_name: str,
    ridx_lock: multiprocessing.Lock,  # type: ignore
    bidx_dict: Dict[str, int],
    bidx_lock: multiprocessing.Lock,  # type: ignore
    config: Config,
):
    """Collects results from the queue and processes them in batches."""

    # Infer the type based on the process_fn signature
    fn_sig = inspect.signature(process_fn)
    first_param_type = fn_sig.parameters[next(iter(fn_sig.parameters))].annotation

    try:
        if first_param_type == pd.DataFrame:
            _queue_batch_processor_df(
                queue,
                process_fn,
                ridx_dict,
                ridx_name,
                ridx_lock,
                bidx_dict,
                bidx_lock,
                config,
            )
        elif first_param_type == List[ReadResult]:
            _queue_batch_processor_list(
                queue,
                process_fn,
                ridx_dict,
                ridx_name,
                ridx_lock,
                bidx_dict,
                bidx_lock,
                config,
            )
        else:
            raise ValueError(
                f"Unsupported process_fn first argument type: {first_param_type}"
            )
    except Exception as e:
        logging.error(f"Error in queue_batch_processor: {e}")
        logging.error(traceback.format_exc())
        raise e


def save_batch_outputs_pass(
    read_results: List[ReadResult], bidx_dict, bidx_lock, config: Config
):
    with bidx_lock:
        bidx = bidx_dict["pass"]
    save_detect_results(
        "pass",
        results=read_results,
        batch_idx=bidx,
        save_fpts=config.output.save_fpts,
        save_dwell_time=config.output.save_dwell_time,
        save_boundaries=config.output.save_boundaries,
        output_dir_boundaries=config.output.output_dir_boundaries,
        output_dir_fpts=config.output.output_dir_fpts,
    )


def save_batch_outputs_fail(
    read_results: List[ReadResult], bidx_dict, bidx_lock, config: Config
):
    with bidx_lock:
        bidx = bidx_dict["fail"]
    save_detect_results(
        "fail",
        results=read_results,
        batch_idx=bidx,
        output_dir_fail=config.output.output_dir_fail,
        save_fpts=False,
        save_dwell_time=False,
        save_boundaries=True,
    )
    del read_results


def save_batch_predictions(
    predictions: pd.DataFrame, bidx_dict, bidx_lock, config: Config
):
    with bidx_lock:
        bidx = bidx_dict["predict"]
    save_predictions(
        predictions,
        os.path.join(
            config.output.output_dir_pred,
            f"barcode_predictions_{bidx}.csv",
        ),
    )
    del predictions


def save_detect_results(
    pass_or_fail: str,
    results: List[ReadResult],
    batch_idx: int,
    save_fpts: bool = True,
    save_dwell_time: bool = False,
    save_boundaries: bool = True,
    output_dir_boundaries: str = "",
    output_dir_fpts: str = "",
    output_dir_fail: str = "",
    **kwargs,
) -> Union[Tuple[np.ndarray, np.ndarray], None]:
    """Save detected boundaries and return read_ids and signals if pass_or_fail == 'pass'."""

    if pass_or_fail == "pass":
        fn = "detected_boundaries"
        dirn = output_dir_boundaries
    elif pass_or_fail == "fail":
        fn = "failed_reads"
        dirn = output_dir_fail
    else:
        msg = f"Invalid pass_or_fail: {pass_or_fail}. Must be 'pass' or 'fail'."
        logging.error(msg)
        raise ValueError(msg)

    if save_boundaries:
        save_detected_boundaries(
            results,  # type: ignore
            os.path.join(dirn, f"{fn}_{batch_idx}.csv"),
            save_fail_reasons=pass_or_fail == "fail",
        )

    if pass_or_fail == "pass" and save_fpts:
        read_ids, signals, _ = save_fpts_signals(
            results,
            os.path.join(output_dir_fpts, f"barcode_fpts_{batch_idx}.npz"),
            save_dwell_time=save_dwell_time,
        )

        return read_ids, signals

    return None


def save_fpts_signals(
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
    num_reads = len(read_ids)

    if save_dwell_time:
        np.savez(
            filename,
            num_reads=num_reads,
            read_ids=read_ids,
            signals=barcode_fpts,
            dwell_times=dwell_times,
        )
    else:
        np.savez(filename, num_reads=num_reads, read_ids=read_ids, signals=barcode_fpts)

    return read_ids, barcode_fpts, dwell_times


def save_predictions(
    predictions: pd.DataFrame,
    filename: str,
) -> None:

    predictions.to_csv(
        filename,
        index=False,
    )


def add_read_id_col_to_predictions(
    predictions: pd.DataFrame,
    read_ids: Union[List[str], np.ndarray],
) -> pd.DataFrame:

    cols = predictions.columns.tolist()

    if "read_id" in cols:
        raise ValueError("'read_id' already in dataframe")

    predictions["read_id"] = read_ids
    return predictions[["read_id", *cols]]


def worker_progress_report(ridx_dict, ridx_lock, save_predictions: bool = True):
    total_reads = get_global_var("total_reads")
    total_reads_set = total_reads != -1

    success_str = "predict" if save_predictions else "pass"

    # Create progress bars for total, failed, and predicted reads
    pbar_total = tqdm(
        desc="Total progress",
        unit="reads",
        position=0,
        total=total_reads if total_reads_set else None,
    )
    pbar_fail = tqdm(
        desc="Failed reads",
        position=1,
        total=total_reads if total_reads_set else None,
        bar_format="{desc}",  # Custom format without speed
    )
    pbar_success = tqdm(
        desc=f"{success_str.capitalize()}ed reads",
        position=2,
        total=total_reads if total_reads_set else None,
        bar_format="{desc}",  # Custom format without speed
    )

    last_fail = 0
    last_success = 0
    time_passed = 0

    while not _STOP_SIGNAL.is_set():
        if time_passed >= 1:  # update every second
            if not total_reads_set:
                total_reads = get_global_var("total_reads")
                total_reads_set = total_reads != -1
                if total_reads_set:
                    for pbar in [pbar_total, pbar_fail, pbar_success]:
                        pbar.total = total_reads

            with ridx_lock:
                n_fail = ridx_dict["done_fail"]
                n_success = ridx_dict[f"done_{success_str}"]
                n_total = n_fail + n_success

                new_fail = n_fail - last_fail
                new_success = n_success - last_success
                new_total = new_fail + new_success

                # Update progress bars
                if new_fail > 0:
                    pbar_fail.update(new_fail)
                if new_success > 0:
                    pbar_success.update(new_success)
                pbar_total.update(new_total)

                # Update descriptions with percentages
                if n_total > 0:
                    fail_pct = (n_fail / n_total) * 100
                    success_pct = (n_success / n_total) * 100
                    pbar_fail.set_description_str(
                        f"Failed reads   {n_fail:,} | {fail_pct:.1f}%"
                    )
                    pbar_success.set_description_str(
                        f"{success_str.capitalize()}ed reads   {n_success:,} | {success_pct:.1f}%"
                    )

                last_fail = n_fail
                last_success = n_success

            time_passed = 0
        time.sleep(0.1)
        time_passed += 0.1

    # Final update
    with ridx_lock:
        n_fail = ridx_dict["done_fail"]
        n_success = ridx_dict[f"done_{success_str}"]
        n_total = n_fail + n_success

        # Update progress bars one last time
        pbar_fail.update(n_fail - last_fail)
        pbar_success.update(n_success - last_success)
        pbar_total.update((n_fail - last_fail) + (n_success - last_success))

        if n_total > 0:
            fail_pct = (n_fail / n_total) * 100
            success_pct = (n_success / n_total) * 100
            pbar_fail.set_description_str(
                f"Failed reads   {n_fail:,} | {fail_pct:.1f}%"
            )
            pbar_success.set_description_str(
                f"{success_str.capitalize()}ed reads   {n_success:,} | {success_pct:.1f}%"
            )

    # Close all progress bars
    pbar_total.close()
    pbar_fail.close()
    pbar_success.close()


def check_total_num_reads(pod5_files: Set[str]) -> int:
    """Count total number of reads across all POD5 files using pod5 inspect."""

    total_reads = 0

    pod5_exec = os.path.join(sys.exec_prefix, "bin", "pod5")
    logging.debug(f"pod5_exec: {pod5_exec}")

    for file in pod5_files:
        try:
            # Run pod5 inspect command and capture output
            result = subprocess.run(
                [pod5_exec, "inspect", "summary", file],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running pod5 inspect on {file}: {e}")
            raise

        try:
            # Get the last line of output
            last_line = result.stdout.strip().split("\n")[-1]

            # Extract number of reads using regex
            match = re.search(r"Found\s+\d+\s+batches,\s+(\d+)\s+reads", last_line)
            if match:
                reads = int(match.group(1))
                total_reads += reads
            else:
                logging.warning(f"Could not parse read count from summary for {file}")

        except subprocess.CalledProcessError as e:
            logging.error(f"Error parsing pod5 inspect output on {file}: {e}")
            raise

    return total_reads


def check_total_num_fpts(npz_files: Set[str]) -> int:
    """Count total number of fpts across all npz files."""

    total_fpts = 0

    python_exec = os.path.join(sys.exec_prefix, "bin", "python")
    logging.debug(f"python_exec: {python_exec}")

    cmd = (
        lambda file: f"""import numpy as np
with np.load("{file}", allow_pickle=True) as npz:
    try:
        print(npz["num_reads"])
    except:
        print(npz["read_ids"].size)
"""
    )

    for file in npz_files:
        try:
            # Run pod5 inspect command and capture output
            result = subprocess.run(
                [python_exec, "-c", cmd(file)],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running python script on {file}: {e}")
            logging.error(f"stdout: {e.stdout}")
            logging.error(f"stderr: {e.stderr}")
            raise

        try:
            num_fpts = int(result.stdout.strip())
            total_fpts += num_fpts

        except subprocess.CalledProcessError as e:
            logging.error(f"Error parsing python script output on {file}: {e}")
            raise

    return total_fpts


# TODO: if read_ids_incl contains read_ids that are not in the files, total_reads will be incorrect
# TODO: if read_ids_excl contains read_ids that are not in the files, total_reads will be incorrect
def set_total_num_reads(
    files: Set[str], read_ids_incl: Set[str], read_ids_excl: Set[str]
) -> None:

    if any(file.endswith(".pod5") for file in files):
        set_total_num_reads_pod5(files, read_ids_incl, read_ids_excl)
    elif any(file.endswith(".npz") for file in files):
        set_total_num_reads_npz(files, read_ids_incl, read_ids_excl)
    else:
        raise ValueError("No POD5 or NPZ files found in the input")


def set_total_num_reads_pod5(
    pod5_files: Set[str], read_ids_incl: Set[str], read_ids_excl: Set[str]
) -> None:
    if read_ids_incl:
        total_reads = len(read_ids_incl)
    else:
        total_reads_pod5 = check_total_num_reads(pod5_files)
        total_reads = total_reads_pod5 - len(read_ids_excl)
    set_global_var("total_reads", total_reads)
    logging.debug(f"total_reads: {total_reads}")


def set_total_num_reads_npz(
    npz_files: Set[str], read_ids_incl: Set[str], read_ids_excl: Set[str]
) -> None:
    if read_ids_incl:
        total_fpts = len(read_ids_incl)
    else:
        total_fpts_npz = check_total_num_fpts(npz_files)
        total_fpts = total_fpts_npz - len(read_ids_excl)
    set_global_var("total_reads", total_fpts)
    logging.debug(f"total_reads: {total_fpts}")


# TODO: move models
def load_model(model_name: str) -> DTW_SVM_Model:
    with pkg_resources.path(model_files, f"{model_name}.joblib") as model_path:
        return joblib.load(model_path)


def run_demux(
    files: Set[str], read_ids_incl: Set[str], read_ids_excl: Set[str], config: Config
) -> Dict[str, Any]:
    """enqueue_minibatches_thread takes care of preloading the signals in the queue.
    Maxsize of queue is set to num_proc. Detect is the most time consuming step, we want to limit the signal preloading to be in sync with the speed of detect.
    """

    save_output_pass = config.output.save_boundaries or config.output.save_fpts
    save_predictions = config.output.save_predictions

    save_pass_thread = None
    save_fail_thread = None
    save_predict_thread = None

    # set global variables for batch indices for continue mode
    set_global_var("bidx_pass", config.batch.bidx_pass)
    set_global_var("bidx_fail", config.batch.bidx_fail)
    set_global_var(
        "bidx_predict", config.batch.bidx_pass
    )  # batch index for predict is the same as for pass

    logging.debug(f"starting mp manager")
    with multiprocessing.Manager() as manager:
        if config.task.predict:
            assert (
                config.classif is not None
            ), "Classification configuration not provided"  # make linter happy
            model_predict = load_model(config.classif.model_name)
        else:
            model_predict = None

        if config.task.preprocess:
            if config.sig_proc.primary_method != "cnn":
                model_detect = None
            else:
                model_detect = load_cnn_model(config.sig_proc.cnn_boundaries.model_name)
        else:
            model_detect = None

        preloaded_minibatch_queue = manager.Queue(maxsize=config.batch.num_proc)
        save_pass_queue = manager.Queue() if save_output_pass else None
        save_fail_queue = manager.Queue() if config.task.preprocess else None
        save_predict_queue = manager.Queue() if save_predictions else None

        ridx_lock = manager.Lock()
        bidx_lock = manager.Lock()

        # Create a shared dictionary that will be shared across processes
        ridx_dict = manager.dict()
        ridx_dict["enqueued"] = 0
        ridx_dict["done_pass"] = 0
        ridx_dict["done_fail"] = 0
        ridx_dict["done_predict"] = 0
        ridx_dict["saved_pass"] = 0
        ridx_dict["saved_fail"] = 0
        ridx_dict["saved_predict"] = 0

        bidx_dict = manager.dict()
        bidx_dict["pass"] = get_global_var("bidx_pass")
        bidx_dict["fail"] = get_global_var("bidx_fail")
        bidx_dict["predict"] = get_global_var("bidx_predict")

        # puts None into the preloaded_minibatch_queue when having preloaded all minibatches
        # and then returns
        # None is the sentinel value that allows the worker processes to exit
        total_reads_thread = threading.Thread(
            target=set_total_num_reads,
            args=(files, read_ids_incl, read_ids_excl),
        )

        target_enqueue_minibatches = (
            worker_enqueue_minibatches_raw
            if config.task.preprocess
            else worker_enqueue_minibatches_fpts
        )

        enqueue_minibatches_thread = threading.Thread(
            target=target_enqueue_minibatches,
            args=(
                files,
                read_ids_incl,
                read_ids_excl,
                config,
                preloaded_minibatch_queue,
                ridx_dict,
                ridx_lock,
            ),
        )

        save_pass_thread = (
            threading.Thread(
                target=queue_batch_processor,
                args=(
                    save_pass_queue,
                    save_batch_outputs_pass,
                    ridx_dict,
                    "saved_pass",
                    ridx_lock,
                    bidx_dict,
                    bidx_lock,
                    config,
                ),
            )
            if save_output_pass
            else None
        )

        save_fail_thread = (
            threading.Thread(
                target=queue_batch_processor,
                args=(
                    save_fail_queue,
                    save_batch_outputs_fail,
                    ridx_dict,
                    "saved_fail",
                    ridx_lock,
                    bidx_dict,
                    bidx_lock,
                    config,
                ),
            )
            if config.task.preprocess
            else None
        )

        # save predictions
        save_predict_thread = (
            threading.Thread(
                target=queue_batch_processor,
                args=(
                    save_predict_queue,
                    save_batch_predictions,
                    ridx_dict,
                    "saved_predict",
                    ridx_lock,
                    bidx_dict,
                    bidx_lock,
                    config,
                ),
            )
            if save_predictions
            else None
        )

        progress_bar_thread = threading.Thread(
            target=worker_progress_report,
            args=(
                ridx_dict,
                ridx_lock,
                save_predictions,
            ),
        )

        logging.debug(f"all threads created")

        total_reads_thread.start()

        if save_output_pass and save_pass_thread:
            save_pass_thread.start()
        if config.task.preprocess and save_fail_thread:
            save_fail_thread.start()
        if save_predictions and save_predict_thread:
            save_predict_thread.start()
        enqueue_minibatches_thread.start()

        progress_bar_thread.start()
        logging.debug("all threads started")

        def handle_completed_future(future):
            try:
                future.result()  # This will raise any exceptions that occurred
            except Exception as e:
                logging.error(f"Error in worker: {e}")
                logging.error(traceback.format_exc())

        # submit jobs to the process pool as minibatches are preloaded
        # the job is finished when the processed signal and predictions have been added to the save queues
        # not done logic is used to manage the number of concurrent jobs and the preloaded queue
        # when a job is started, it takes from the preloaded queue
        # submitting all jobs at once would require a large memory footprint for the preloaded queue

        target_process_pool = (
            worker_detect_and_predict_on_preloaded_signals
            if config.task.preprocess
            else worker_predict_on_preloaded_fpts
        )
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=config.batch.num_proc,
        ) as process_pool:
            not_done_futures = set()
            not_done_count = 0

            all_enqueued = False

            while not all_enqueued:
                # First, check for completed futures
                done_futures, not_done_futures = concurrent.futures.wait(
                    not_done_futures,
                    timeout=0.001,  #  prevent busy waiting
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )

                # Handle completed futures
                for future in done_futures:
                    handle_completed_future(future)
                    not_done_count -= 1

                # Then try to submit new work if we have capacity
                while not_done_count < 1.1 * config.batch.num_proc:
                    try:
                        minibatch = preloaded_minibatch_queue.get()
                    except queue.Empty:
                        if not not_done_futures:  # If no pending work and queue empty
                            all_enqueued = True
                            break
                        continue

                    if minibatch is None:
                        all_enqueued = True
                        break  # all minibatches have been submitted to the worker pool

                    future = process_pool.submit(
                        target_process_pool,
                        minibatch,
                        model_predict,
                        model_detect,
                        config,
                        ridx_dict,
                        ridx_lock,
                        save_pass_queue,
                        save_fail_queue,
                        save_predict_queue,
                    )
                    not_done_futures.add(future)
                    not_done_count += 1

        logging.debug("All detect and predict workers have been launched")
        enqueue_minibatches_thread.join()
        logging.debug("Enqueue minibatches thread finished")

        # Wait for all remaining futures to complete
        concurrent.futures.wait(not_done_futures)

        # Handle any remaining completed futures
        for future in not_done_futures:
            handle_completed_future(future)

        logging.debug("All detect and predict workers have finished")

        # signal the save threads that no more minibatches will be added to the queues
        if save_fail_queue:
            save_fail_queue.put(None)
        if save_pass_queue:
            save_pass_queue.put(None)
        if save_predict_queue:
            save_predict_queue.put(None)

        logging.debug("Signalled save threads to finish")

        logging.debug("Waiting for save threads to finish")

        total_reads_thread.join()
        if save_fail_thread:
            save_fail_thread.join()
        if save_pass_thread:
            save_pass_thread.join()
        if save_predict_thread:
            save_predict_thread.join()

        logging.debug("Save threads finished")

        _STOP_SIGNAL.set()
        progress_bar_thread.join()
        logging.debug("Progress bar thread finished")

        success_str = "predict" if save_predictions else "pass"
        if ridx_dict[f"done_{success_str}"] > 0:
            if config.task.preprocess:
                pretext = "Adapter was successfully detected in"
            else:
                pretext = "Predictions were successfully generated in"
            logging.info(
                f"{pretext} {ridx_dict[f'done_{success_str}']} / "
                f"{ridx_dict['done_fail'] + ridx_dict[f'done_{success_str}']} reads  "
                f"({ridx_dict[f'done_{success_str}'] / ridx_dict['enqueued'] * 100:.2f}%)."
            )
        else:
            logging.info(f"All reads failed.")

        return dict(ridx_dict)
