"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import importlib.resources as pkg_resources

import toml
from adapted.config.base import load_nested_config_from_file

from warpdemux.models import model_files
from warpdemux.config import config_files
from warpdemux.config.sig_proc import SigProcConfig


# TODO: update model naming to WDX[xx]_rna00[x]_130bps@v[x].[x].[x]
def get_model_spc_config(model_name: str) -> SigProcConfig:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    return get_config(model_config["spc"])


def get_model_spc_live_config(model_name: str) -> SigProcConfig:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    return get_config(model_config["spc_live"])


def get_model_num_bcs(model_name: str) -> int:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    return model_config["num_bcs"]


def get_config(config_name: str) -> SigProcConfig:
    with pkg_resources.path(config_files, f"{config_name}.toml") as config_path:
        return load_nested_config_from_file(config_path, SigProcConfig)
