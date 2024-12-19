"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import logging
import os
import sys
import time
from typing import Any, Dict

import warpdemux  # sets adapted submodule path for editable install
from warpdemux.file_proc import handle_previous_results, run_demux
from warpdemux.logger import setup_logger
from warpdemux.parser import parse_args


def main(args=None) -> Dict[str, Any]:

    config = parse_args()
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

    read_ids_excl = set()
    read_ids_incl = set()
    if config.input.continue_from:

        # Preprocess input_read_ids into batches
        logging.info(f"Indexing previous results...")
        start_time = time.time()

        read_ids_excl = handle_previous_results(
            config,
        )
        logging.info(f"Indexing took: {time.time() - start_time:.2f} seconds")
        logging.info(f"Found {len(read_ids_excl)} previously processed reads.")

    files = set(config.input.files)
    read_ids_incl = set(config.input.read_ids)
    config.input.files = (
        []
    )  # clear to prevent long lists being copied around to all processes
    config.input.read_ids = (
        []
    )  # clear to prevent long lists being copied around to all processes

    logging.info(f"Config.output: {config.output.dict() }")
    logging.info(f"Config.batch: {config.batch.dict() }")
    logging.info(f"Config.task: {config.task.dict() }")
    if config.task.predict:
        logging.info(f"Config.classif: {config.classif.dict() }")

    if config.task.preprocess:
        # save spc that were used
        config.sig_proc.to_toml(os.path.join(config.output.output_dir, "config.toml"))

    command_verb = {
        "demux": "demultiplexing",
        "prep": "preprocessing",
        "predict": "predicting",
    }[config.task.command]
    logging.info(f"Starting {command_verb}...")
    start_time = time.time()
    ridx_dict = run_demux(files, read_ids_incl, read_ids_excl, config)
    end_time = time.time()

    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.info(
        f"{command_verb.capitalize()} took: {int(hours):02d}:{int(minutes):02d}:{int(seconds):02d} (HH:MM:SS)"
    )

    return ridx_dict


if __name__ == "__main__":
    main()
