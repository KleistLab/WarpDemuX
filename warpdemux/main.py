"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import os
import sys
import time
import logging

import numpy as np
from adapted.file_proc.file_proc import get_file_read_id_map, process
from adapted.io_utils import lexsort_num_suffix
from tqdm import tqdm

import warpdemux  # sets adapted submodule path for editable install
from warpdemux.file_proc import (
    barcode_fpt_wrapper,
    load_model,
    predict_batch,
    save_detect_fpts_batch,
    resegment_wrapper,
)
from warpdemux.parser import parse_args
from warpdemux.logger import setup_logger


def main(args=None):

    command, config = parse_args()
    setup_logger(os.path.join(config.output.output_dir, "warpdemux.log"))

    logging.info(f"Command: {' '.join(sys.argv)}")
    logging.info(f"Saving output to: {config.output.output_dir}")

    print_files = (
        config.input.files[: min(3, len(config.input.files))]
        + [f"..."]
        + config.input.files[-min(3, len(config.input.files)) :]
        if len(config.input.files) > 3
        else config.input.files
    )
    print_files_str = "\n".join(print_files)
    logging.info(f"Input filenames:\n{print_files_str}")
    logging.info(f"Total number of input files: {len(config.input.files)}")

    os.makedirs(config.output.output_dir, exist_ok=True)
    fpts_files = config.input.files
    if command == "fpts" or (command == "demux" and not config.input.preprocessed):
        # report config
        logging.info("SigProcConfig:")
        config.sig_proc.pretty_print(file=logging.getLogger().handlers[0].stream)  # type: ignore

        # Preprocess input_read_ids into batches
        logging.info(f"Indexing read IDs...")
        start_time = time.time()

        file_read_id_map = get_file_read_id_map(config)
        logging.info(f"Indexing took: {time.time() - start_time:.2f} seconds")

        config.input.files = []  # no longer needed, save space
        config.input.read_ids = []  # no longer needed, save space

        # save spc that were used
        config.sig_proc.to_toml(os.path.join(config.output.output_dir, "config.toml"))

        # if command == "fpts":  # only preprocess
        process(
            file_read_id_map=file_read_id_map,
            config=config,
            task_fn=barcode_fpt_wrapper,
            results_fn=save_detect_fpts_batch,
        )
        # NOTE: combining is predict and preprocess messes with mp and is slow. Don't do it.
        # elif not config.input.preprocessed:  # command == "demux", preprocess and predict
        #     process(
        #         file_read_id_map=file_read_id_map,
        #         config=config,
        #         task_fn=barcode_fpt_wrapper,
        #         results_fn=get_save_detect_fpts_and_predict_batch_fn(config),
        #     )

        # TODO: check if there are already any fpts in output dir, might lead to unexpected results..
        fpts_files = lexsort_num_suffix(
            [
                os.path.join(config.output.output_dir, f)
                for f in os.listdir(config.output.output_dir)
                if f.startswith("barcode_fpts_") and f.endswith(".npz")
            ]
        )
    elif command == "resegment":
        logging.info("ATTENTION: only using relevant configs for resegmentation. ")

        logging.info("Resegmentation configs:")
        logging.info("Normalization:")
        logging.info(config.sig_proc.sig_norm)
        logging.info("Extraction:")
        logging.info(config.sig_proc.sig_extract)
        logging.info("Segmentation:")
        logging.info(config.sig_proc.segmentation)

        logging.info(f"Indexing read IDs...")
        start_time = time.time()

        file_read_id_map = get_file_read_id_map(config)
        logging.info(f"Indexing took: {time.time() - start_time:.2f} seconds")
        config.input.files = []  # no longer needed, save space
        config.input.read_ids = []  # no longer needed, save space

        # save spc that were used
        config.sig_proc.to_toml(os.path.join(config.output.output_dir, "config.toml"))

        # add a NOTE.txt file to output dir with the following content:
        # "This directory contains resegmented signals. The signals are not raw signals, but rather normalized and extracted signals."
        with open(
            os.path.join(config.output.output_dir, "RESEGMENTATION_NOTE.txt"), "w"
        ) as f:
            f.write(
                "This directory contains resegmented signals. "
                "Detection parameters in the config.toml may not reflect the ones used during detection."
            )

        process(
            file_read_id_map=file_read_id_map,
            config=config,
            task_fn=resegment_wrapper,
            results_fn=save_detect_fpts_batch,
        )

    if command == "demux":
        logging.info("Loading model...")
        assert (
            config.classif is not None
        ), "Classification config is required for classification."
        model = load_model(config.classif.model_name)

        # check for existing predictions in output_dir, skip the corresponding fpts files
        existing_predictions = {
            int(f.split("_")[-1].split(".")[0])
            for f in os.listdir(config.output.output_dir)
            if f.startswith("barcode_predictions_")
        }
        batch_offset = len(existing_predictions) + 1
        fpts_files = [
            f
            for f in fpts_files
            if int(f.split("_")[-1].split(".")[0]) not in existing_predictions
        ]
        remaining_str = "remaining" if batch_offset > 0 else ""
        logging.info(
            f"Predicting barcodes for {len(fpts_files)} {remaining_str} fingerprint files..."
        )
        pbar = tqdm(
            enumerate(fpts_files),
            total=len(fpts_files),
            desc="Predicting barcodes",
            file=sys.stdout,
        )
        for index, fpts_fpath in pbar:
            npz = np.load(fpts_fpath, allow_pickle=True)
            signals = npz["signals"][
                :, -config.sig_proc.segmentation.barcode_num_events :
            ]
            read_ids = npz["read_ids"]

            predict_batch(
                model=model,
                signals=signals,
                read_ids=read_ids,
                batch_idx=index + batch_offset,
                output_dir=config.output.output_dir,
                classif_config=config.classif,
            )
        pbar.close()
        logging.info(f"Predictions completed in {pbar.format_dict['elapsed']:.2f}s")

    logging.info("Done.")


if __name__ == "__main__":
    main()
