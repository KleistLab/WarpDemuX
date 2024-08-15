import os
import sys
import time

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


def main(args=None):
    print("Command executed:")
    print(" ".join(sys.argv))

    command, config = parse_args()

    print("Saving output to:", config.output.output_dir)

    print_files = (
        config.input.files[: min(3, len(config.input.files))]
        + [f"..."]
        + config.input.files[-min(3, len(config.input.files)) :]
        if len(config.input.files) > 3
        else config.input.files
    )
    print_files_str = "\n".join(print_files)
    print(f"Input Filenames:\n{print_files_str}")
    print(f"Total files: {len(config.input.files)}")

    os.makedirs(config.output.output_dir, exist_ok=True)
    fpts_files = config.input.files
    if command == "fpts" or (command == "demux" and not config.input.preprocessed):
        # report config
        print("SigProcConfig:")
        config.sig_proc.pretty_print()

        # Preprocess input_read_ids into batches
        print(f"Preprocessing read IDs for {len(config.input.files)} files")
        start_time = time.time()

        file_read_id_map = get_file_read_id_map(config)
        print(f"Time taken: {time.time() - start_time:.2f} seconds")

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
        print("ATTENTION: only using relevant configs for resegmentation. ")

        print("Resegmentation configs:")
        print("Normalization:")
        print(config.sig_proc.sig_norm)
        print("Extraction:")
        print(config.sig_proc.sig_extract)
        print("Segmentation:")
        print(config.sig_proc.segmentation)

        print(f"Preprocessing read IDs for {len(config.input.files)} files")
        start_time = time.time()

        file_read_id_map = get_file_read_id_map(config)
        print(f"Time taken: {time.time() - start_time:.2f} seconds")
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
        print("Loading model...")
        assert (
            config.classif is not None
        ), "Classification config is required for classification."
        model = load_model(config.classif.model_name)

        for index, fpts_fpath in tqdm(
            enumerate(fpts_files), total=len(fpts_files), desc="Predicting barcodes"
        ):
            npz = np.load(fpts_fpath, allow_pickle=True)
            signals = npz["signals"][
                :, -config.sig_proc.segmentation.barcode_num_events :
            ]
            read_ids = npz["read_ids"]

            predict_batch(
                model=model,
                signals=signals,
                read_ids=read_ids,
                batch_idx=index,
                output_dir=config.output.output_dir,
                classif_config=config.classif,
            )


if __name__ == "__main__":
    main()
