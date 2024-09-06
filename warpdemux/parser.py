import argparse
import os
from datetime import datetime
from typing import Tuple
from adapted.output import np

import pandas as pd
from adapted.config.base import load_nested_config_from_file
from adapted.config.file_proc import (
    BatchConfig,
    TaskConfig as DetectTaskConfig,
)
from adapted.io_utils import input_to_filelist

from warpdemux.config.utils import get_model_spc_config, get_chemistry_specific_config

from warpdemux.config.classification import ClassificationConfig
from warpdemux.config.config import Config
from warpdemux.config.file_proc import InputConfig, OutputConfig, ResegmentTaskConfig
from warpdemux.config.sig_proc import SigProcConfig
from warpdemux._version import __version__


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


parent_parser = argparse.ArgumentParser(
    add_help=False,
)

parent_parser.add_argument(
    "--input",
    type=str,
    nargs="+",
    help=("Input file(s) or directory(s)."),
)
parent_parser.add_argument(
    "--output",
    type=str,
    default=None,
    help="Path to where the run output folder should be created. Default is the current working directory.",
)

parent_parser.add_argument(
    "--save_dwell_times",
    type=str2bool,
    default=False,
    help="Whether to save segment dwell times in the output file. Default is False.",
)


parent_parser.add_argument(
    "-p",
    "--num_proc",
    type=int,
    default=None,
    help=(
        "Number of num_proc to use for parallel processing. If not specified, all"
        " available cores will be used."
    ),
)

parent_parser.add_argument(
    "--batch_size",
    type=int,
    default=4000,
    help=("Number of reads per output file."),
)

parent_parser.add_argument(
    "--minibatch_size",
    type=int,
    default=50,
    help=(
        "Number of reads per worker. These reads are loaded into memory prior to"
        " processing. Choose depending on the max_obs_adapter value and the amount"
        " of memory available."
    ),
)

parent_parser.add_argument(
    "--create_subdir",
    type=str2bool,
    default=True,
    help="Whether to create a subdirectory for the output. Default is True.",
)

parser = argparse.ArgumentParser(description="Process files.", parents=[parent_parser])


subparsers = parser.add_subparsers(title="workflows", dest="command")

fpts_parser = subparsers.add_parser(
    "fpts",
    help="Generate barcode fingerprints from raw signal, no demultiplexing.",
    parents=[parent_parser],
)
reseg_parser = subparsers.add_parser(
    "resegment",
    help="Re-segment raw signal to barcode fingerprints, no adapter end detection or demultiplexing.",
    parents=[parent_parser],
)
demux_parser = subparsers.add_parser(
    "demux",
    help="Demultiplex raw signal or preprocessed barcode fingerprints.",
    parents=[parent_parser],
)

for _parser in [fpts_parser, reseg_parser]:
    _parser.add_argument(
        "--config",
        type=str,
        help="Path to a valid configuration toml to use. See warpdemux/config/config_files/ for options.",
    )

    _parser.add_argument(
        "--chemistry",
        type=str,
        choices=["RNA002", "RNA004"],
        help="Specify the chemistry to use. If provided, --config is not required and will be ignored if both are provided.",
    )

for _parser in [fpts_parser, demux_parser]:
    _parser.add_argument(
        "--read_id_csv",
        type=str,
        default=None,
        help=(
            "Path to a csv file containing read IDs to be processed. Should contain a"
            " 'read_id' column."
        ),
    )

    _parser.add_argument(
        "--read_id_csv_colname",
        type=str,
        default="read_id",
        help=(
            "Column name in 'read_id_csv' containing the read IDs to be processed. Defaults"
            " to 'read_id'. This argument is ignored if '--preprocessed' is set."
        ),
    )
    _parser.add_argument(
        "--DEBUG",
        action="store_true",
        help="Save ADAPTed traces for adapter detection debug.",
    )


reseg_parser.add_argument(
    "--prev_results",
    type=str,
    nargs="+",
    required=True,
    help=(
        "Directory or path(s) to csv file(s) with previous detection results to be processed. "
        "This should be the 'detected_boundaries_[batch_id].csv' files from a previous run. "
    ),
)


demux_parser.add_argument(
    "--preprocessed",
    action="store_true",
    help=(
        "If specified, the input files are assumed to be preprocessed and the"
        " preprocessing step is skipped. Input files should be `barcode_fpts.npz`"
        " files obtained through `read_data_handler.process_files`."
    ),
)

demux_parser.add_argument(
    "--model_name",
    type=str,
    default="WDX12_rna002_v0_4_3",
    help="Name of the model to use for classification. Default is `WDX12_rna002_v0_4_3`.",
)


demux_parser.add_argument(
    "--classif_block_size",
    type=int,
    default=500,
    help=(
        "Block size for pairwise distance matrix calculations during classification. "
        "Affects classification speed. Ideal settings depend on data size and available memory. "
        "Default is 500."
    ),
)


def parse_args() -> Tuple[str, Config]:
    args = parser.parse_args()

    if args.command == "fpts" or args.command == "resegment":
        if not args.config and not args.chemistry:
            parser.error("Either --config or --chemistry must be provided.")

    read_ids = []
    prev_results = []
    if args.command == "resegment":
        endswiths = [".csv"]
        basenameprefix = str("detected_boundaries")
        prev_results = input_to_filelist(
            args.prev_results, endswiths=endswiths, basenameprefix=basenameprefix
        )
        read_ids = np.hstack(
            [np.asarray(pd.read_csv(f).read_id) for f in prev_results]
        ).tolist()

    elif args.read_id_csv is not None:
        read_ids = pd.read_csv(
            args.read_id_csv,
        )[args.read_id_csv_colname].values

    preprocessed = False
    if args.command == "demux" and args.preprocessed:
        preprocessed = True
        endswiths = [".npz"]
        basenameprefix = str("barcode_fpts")
    else:
        endswiths = [".pod5"]
        basenameprefix = ""

    files = input_to_filelist(
        args.input, endswiths=endswiths, basenameprefix=basenameprefix
    )

    if len(files) == 0:
        print("No valid input files found.")
        print("Provided path: {}".format(args.input))
        exit(1)

    if args.output is None:
        args.output = os.getcwd()

    # create run dir
    if args.create_subdir:
        run_dir_name = (
            "warpdemux_"
            + __version__.replace(".", "_")
            + "_"
            + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        )
        run_dir = os.path.join(args.output, run_dir_name)
    else:
        run_dir = args.output
    os.makedirs(run_dir, exist_ok=True)

    input_config = InputConfig(
        files=files,
        read_ids=read_ids,
        preprocessed=preprocessed,
    )

    if args.command == "fpts" or args.command == "demux":
        debug_flag = args.DEBUG
        task_config = DetectTaskConfig(
            llr_return_trace=debug_flag,
        )
    elif args.command == "resegment":
        debug_flag = False
        detected_boundaries_dict = {
            str(k): v
            for prev_res in prev_results
            for k, v in pd.read_csv(prev_res)
            .set_index("read_id")
            .to_dict(orient="index")
            .items()
        }  # double loop to make type checker happy
        task_config = ResegmentTaskConfig(
            detect_result_dict=detected_boundaries_dict,
        )
    else:
        raise ValueError(f"Invalid command: {args.command}")

    batch_config = BatchConfig(
        num_proc=args.num_proc,
        batch_size=args.batch_size,
        minibatch_size=args.minibatch_size,
    )

    output_config = OutputConfig(
        output_dir=run_dir,
        save_llr_trace=args.command != "resegment" and args.DEBUG,
        save_dwell_time=args.save_dwell_times,
        save_boundaries=args.command != "resegment",
    )

    if args.command == "demux":
        spc = get_model_spc_config(args.model_name)
        cc = ClassificationConfig(
            model_name=args.model_name,
            block_size=args.classif_block_size,
            num_proc=args.num_proc,
        )
    else:
        if args.config:
            # TODO: document that the provided config file should also include the adapted config sections!
            spc = load_nested_config_from_file(args.config, SigProcConfig)
        else:
            spc = get_chemistry_specific_config(
                args.chemistry, load_adapted_config_first=True
            )
        cc = None

    spc.update_sig_preload_size()

    config = Config(
        input=input_config,
        task=task_config,
        batch=batch_config,
        output=output_config,
        sig_proc=spc,
        classif=cc,
    )

    print(config.sig_proc.sig_norm)
    print(config.sig_proc.sig_extract)
    print(config.sig_proc.segmentation)

    return args.command, config
