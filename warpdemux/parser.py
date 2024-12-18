"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import argparse
import ast
import datetime
import json
import os
import shutil
import sys
import uuid
from argparse import Action
from typing import Any, Iterable, List, NamedTuple, Optional, Type, Union

import pandas as pd
import toml
from adapted.io_utils import input_to_filelist

from warpdemux.config.classification import ClassificationConfig
from warpdemux.config.config import Config
from warpdemux.config.file_proc import (
    BatchConfig,
    InputConfig,
    OutputConfig,
    TaskConfig,
)
from warpdemux.config.utils import get_model_spc_config


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


class ArgumentParams(NamedTuple):
    name_or_flags: "tuple[str, ...]"
    action: Optional[Union[str, Type[Action]]] = None
    nargs: Optional[Union[int, str]] = None
    const: Any = None
    default: Any = None
    type: Optional[Any] = None
    choices: Optional[Iterable[Any]] = None
    required: bool = False
    help: Optional[str] = None
    metavar: Optional[Union[str, "tuple[str, ...]"]] = None
    dest: Optional[str] = None

    def add_to_parser(
        self, parser: argparse.ArgumentParser, required: Optional[bool] = None
    ):
        params = self._asdict()
        del params["name_or_flags"]
        if required is not None:
            params["required"] = required
        parser.add_argument(*self.name_or_flags, **params)


io_arguments = {
    "input": ArgumentParams(
        name_or_flags=("--input", "-i"),
        type=str,
        nargs="+",
        required=True,
        help=("Input file(s) or directory(s)."),
    ),
    "output": ArgumentParams(
        name_or_flags=("--output", "-o"),
        type=str,
        default=None,
        help=(
            "Path to where the run output folder should be created. "
            "Default is the current working directory."
        ),
    ),
    "save_fpts": ArgumentParams(
        name_or_flags=("--save_fpts",),
        type=str2bool,
        default=False,
        help=(
            "Whether to save the barcode fingerprints as .npz files. "
            "Default is False."
        ),
    ),
    "save_dwell_times": ArgumentParams(
        name_or_flags=("--save_dwell_times",),
        type=str2bool,
        default=False,
        help=(
            "Whether to save the dwell times per segment along with the fingerprints. "
            "Ignored if --save_fpts is False. Default is False."
        ),
    ),
    "save_boundaries": ArgumentParams(
        name_or_flags=("--save_boundaries",),
        type=str2bool,
        default=False,
        help=("Whether to save the boundaries as .csv files. Default is False."),
    ),
    "read_id_csv": ArgumentParams(
        name_or_flags=("--read_id_csv", "-r"),
        type=str,
        default=None,
        help=(
            "Path to a csv file containing read IDs to be processed. "
            "Should contain a 'read_id' column."
        ),
    ),
    "read_id_csv_colname": ArgumentParams(
        name_or_flags=("--read_id_csv_colname", "-c"),
        type=str,
        default="read_id",
        help=(
            "Column name in 'read_id_csv' containing the read IDs to be processed. "
            "Defaults to 'read_id'."
        ),
    ),
}

config_arguments = {
    "model_name": ArgumentParams(
        name_or_flags=("--model_name", "-m"),
        type=str,
        default="WDX4_rna004_v0_4_4",
        help=(
            "Name of the model to use for classification. "
            "The model name is directly linked to the preprocessing config that will be used, "
            "see details in models/model_files/config.toml."
        ),
    ),
    "export": ArgumentParams(
        name_or_flags=("--export", "-e"),
        type=str,
        default=None,
        help=(
            "Export custom configuration in format 'section.param=value1,section.param=value2'. "
            "Example: '--export cnn_boundaries.fallback_to_llr_short_reads=true,cnn_boundaries.polya_cand_k=15'. "
            "Alternatively, you can provide a path to a toml file containing the custom configuration. "
            "Example: '--export /path/to/config.toml'. "
            "Exported values will override the default values in the model config. None specified values will be left unchanged. "
            "Only use this argument if you know what you are doing. "
            "Default is None."
        ),
    ),
}

batch_arguments = {
    "ncores": ArgumentParams(
        name_or_flags=("--ncores", "-j"),
        type=int,
        help=("Number of cores to use for parallel processing."),
    ),
    "batch_size_output": ArgumentParams(
        name_or_flags=("--batch_size_output", "-b"),
        type=int,
        default=4000,
        help=("Number of reads per output file. Default is 4000."),
    ),
    "minibatch_size": ArgumentParams(
        name_or_flags=("--minibatch_size", "-s"),
        type=int,
        default=1000,
        help=("Number of reads per minibatch. Default is 1000."),
    ),
}

io_parser = argparse.ArgumentParser(
    add_help=False,
)
config_parser = argparse.ArgumentParser(
    add_help=False,
)
batch_parser = argparse.ArgumentParser(
    add_help=False,
)


for io_param in io_arguments:
    io_arguments[io_param].add_to_parser(io_parser)

for config_param in config_arguments:
    config_arguments[config_param].add_to_parser(config_parser)

for batch_param in batch_arguments:
    batch_arguments[batch_param].add_to_parser(batch_parser)

parser = argparse.ArgumentParser(
    description="WarpDemuX: Adapter barcode classification for nanopore direct RNA sequencing.",
)

subparsers = parser.add_subparsers(title="workflows", dest="command")

demux_parser = subparsers.add_parser(
    "demux",
    help="Demultiplex raw signal or preprocessed barcode fingerprints.",
    parents=[io_parser, config_parser, batch_parser],
)

prep_parser = subparsers.add_parser(
    "prep",
    help="Prepare data for WarpDemuX. This ignores most model-specific parameters.",
    parents=[io_parser, config_parser, batch_parser],
)

continue_parser = subparsers.add_parser(
    "continue",
    help="Continue from a previous (incomplete) run.",
    parents=[batch_parser],
)

continue_parser.add_argument(
    "continue_from",
    type=str,
    help="Path to a previous WarpDemuX output directory to continue processing from.",
)

parser_retry = subparsers.add_parser(
    "retry",
    help="Retry processing failed reads.",
    parents=[batch_parser],
)

parser_retry.add_argument(
    "retry_from",
    type=str,
    help="Path to a folder containing the results of a previous run.",
)


def parse_export_string(export_str):
    if export_str is None:
        return None

    if os.path.exists(export_str):
        print(f"Loading custom configuration from {export_str}", file=sys.stderr)
        toml_data = toml.load(export_str)
        export_vals = {}
        for section, params in toml_data.items():
            for attr, value in params.items():
                export_vals[(section, attr)] = value

        return export_vals
    else:
        export_vals = {}
        try:
            for pair in export_str.split(","):
                if "=" not in pair:
                    raise ValueError(f"Missing '=' in parameter pair: {pair}")

                key, value = pair.split("=", maxsplit=1)
                key = key.strip()
                if "." not in key:
                    raise ValueError(f"Missing '.' in parameter key: {key}")

                section, attr = key.split(".", maxsplit=1)
                if not section or not attr:
                    raise ValueError(f"Empty variable or attribute in key: {key}")

                # interpret value as bool, float or int if possible
                try:
                    # Handle boolean strings case-insensitively
                    if value.strip().lower() in ("true", "false"):
                        value = value.strip().lower() == "true"
                    else:
                        value = ast.literal_eval(value)
                except (ValueError, SyntaxError):
                    value = value.strip()

                export_vals[(section.strip(), attr.strip())] = value
            return export_vals
        except ValueError as e:
            raise argparse.ArgumentTypeError(
                f"Invalid export string format: {str(e)}. "
                "Must be in format 'section.attr=value1,section.attr2=value2'"
            )


def parse_args(in_args: Optional[List[str]] = None) -> Config:

    args = in_args or sys.argv[1:]
    args = parser.parse_args(args)

    if args.command == "continue":
        try:
            # load the parser arguments from the command.json file
            with open(os.path.join(args.continue_from, "command.json"), "r") as f:
                command_dict = json.load(f)
        except FileNotFoundError:
            parser.error(
                "No command.json file found in the continue_from directory. "
                "Please provide a valid continue-from directory."
            )

        # create a backup of the command.json file
        shutil.copy(
            os.path.join(args.continue_from, "command.json"),
            os.path.join(args.continue_from, "command_previous.json"),
        )

        if "retry_from" in command_dict:
            del command_dict["retry_from"]

        run_dir = args.continue_from
        args.__dict__.update(command_dict)  # update all params, including `command`
    elif args.command == "retry":
        try:
            # load the parser arguments from the command.json file
            with open(os.path.join(args.retry_from, "command.json"), "r") as f:
                command_dict = json.load(f)
        except FileNotFoundError:
            parser.error(
                "No command.json file found in the retry_from directory. "
                "Please provide a valid continue-from directory."
            )

        # create a backup of the command.json file
        shutil.copy(
            os.path.join(args.retry_from, "command.json"),
            os.path.join(args.retry_from, "command_previous.json"),
        )

        # if continue_from is set from a previous run, remove it from args
        if "continue_from" in command_dict:
            del command_dict["continue_from"]

        run_dir = args.retry_from
        args.__dict__.update(command_dict)  # update all params, including `command`

    else:
        args.output = args.output or os.getcwd()

        run_dir = os.path.join(
            args.output,
            "warpdemux_"
            + args.model_name
            + "_"
            + datetime.datetime.now().strftime("%Y%m%d_%H%M")
            + "_"
            + str(uuid.uuid4())[:8],
        )

    read_ids = []

    if args.read_id_csv is not None:
        read_ids = pd.read_csv(
            args.read_id_csv,
        )[args.read_id_csv_colname].values

    endswiths = [".pod5"]
    basenameprefix = ""

    files = input_to_filelist(
        args.input, endswiths=endswiths, basenameprefix=basenameprefix
    )

    if len(files) == 0:
        print("No valid input files found.")
        print("Provided path: {}".format(args.input))
        exit(1)

    input_config = InputConfig(
        files=files,
        read_ids=read_ids,
        continue_from=args.continue_from if "continue_from" in args else "",
        retry_from=args.retry_from if "retry_from" in args else "",
    )

    num_proc = (args.ncores or os.cpu_count()) or -1

    batch_config = BatchConfig(
        num_proc=num_proc,
        batch_size_output=args.batch_size_output,
        minibatch_size=args.minibatch_size,
    )

    if args.command == "prep":
        output_config = OutputConfig(
            output_dir=run_dir,
            save_dwell_time=args.save_dwell_times,
            save_fpts=True,
            save_boundaries=args.save_boundaries,
            save_predictions=False,
            output_subdir_fail=(
                "failed_reads_retry" if "retry_from" in args else "failed_reads"
            ),
        )
        task_config = TaskConfig(
            predict=False,
        )
    else:
        output_config = OutputConfig(
            output_dir=run_dir,
            save_dwell_time=args.save_dwell_times,
            save_fpts=args.save_fpts,
            save_boundaries=args.save_boundaries,
            output_subdir_fail=(
                "failed_reads_retry" if "retry_from" in args else "failed_reads"
            ),
        )
        task_config = TaskConfig(
            predict=True,
        )

    spc = get_model_spc_config(args.model_name)

    if args.export is not None:
        update_dict = parse_export_string(args.export)
        assert update_dict is not None
        for (section, attr), value in update_dict.items():
            section_dict = getattr(
                spc, section
            ).copy()  # Create a copy of the section dictionary
            section_dict[attr] = value  # Update the value
            setattr(spc, section, section_dict)  # Set the updated dictionary

    spc.update_primary_method()
    spc.update_sig_preload_size()

    cc = ClassificationConfig(
        model_name=args.model_name,
    )

    config = Config(
        input=input_config,
        batch=batch_config,
        output=output_config,
        sig_proc=spc,
        classif=cc,
        task=task_config,
    )

    os.makedirs(run_dir, exist_ok=True)

    # Create command.json file
    command_dict = vars(args)
    command_json_path = os.path.join(run_dir, "command.json")
    with open(command_json_path, "w") as f:
        json.dump(command_dict, f, indent=2)

    return config
