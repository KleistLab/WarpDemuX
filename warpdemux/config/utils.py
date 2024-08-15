"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import importlib.resources as pkg_resources
from typing import Optional, Union, MutableMapping, Any

import toml
from adapted.config.base import load_nested_config_from_file, nested_config_from_dict
from adapted.config.sig_proc import (
    config_name_to_dict as adapted_config_name_to_dict,
    chemistry_specific_config_name as adapted_chemistry_specific_config_name,
)


from warpdemux.models import model_files
from warpdemux.config import config_files
from warpdemux.config.sig_proc import SigProcConfig

from warpdemux import __version__


# TODO: update model naming to WDX[xx]_rna00[x]_xxbps@v[x].[x].[x]
def get_model_spc_config(model_name: str) -> SigProcConfig:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    sqk = model_config["SQK"]  # sequencing kit, RNA002 or RNA004
    adapted_config_name = chemistry_specific_config_name(
        sqk
    )  # latest version based on submodule hash
    adapted_config_dict = adapted_config_name_to_dict(adapted_config_name)
    warpdemux_config_dict = config_name_to_dict(model_config["spc"])
    adapted_config_dict.update(
        warpdemux_config_dict
    )  # merge the two dictionaries, warpdemux_config_dict overwrites adapted_config_dict

    return nested_config_from_dict(adapted_config_dict, SigProcConfig)


def get_model_spc_live_config(model_name: str) -> SigProcConfig:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    return get_config(model_config["spc_live"])


def get_model_num_bcs(model_name: str) -> int:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    return model_config["num_bcs"]


def chemistry_specific_config_name(
    chemistry: str, version: Optional[str] = None
) -> str:
    if version is None:
        version = __version__
    speed = config_files.speeds[chemistry.lower()]

    return f"{chemistry.lower()}_{speed}@v{version}"


def get_chemistry_specific_config(
    chemistry: str,
    version: Optional[str] = None,
    load_adapted_config_first: bool = True,
) -> SigProcConfig:
    if chemistry.lower() not in ["rna002", "rna004"]:
        raise ValueError(f"Unknown chemistry: {chemistry}")
    if version is None:
        version = __version__

    if load_adapted_config_first and chemistry is None:
        raise ValueError(
            "chemistry must be provided if load_adapted_config_first is True"
        )
    if load_adapted_config_first and chemistry is not None:
        adapted_config_name = adapted_chemistry_specific_config_name(chemistry)
        adapted_config_dict = adapted_config_name_to_dict(adapted_config_name)
        warpdemux_config_dict = config_name_to_dict(
            chemistry_specific_config_name(chemistry, version)
        )
        adapted_config_dict.update(warpdemux_config_dict)
        return nested_config_from_dict(adapted_config_dict, SigProcConfig)
    else:
        return get_config(chemistry_specific_config_name(chemistry, version))


def get_config(
    config_name: str,
    load_adapted_config_first: bool = True,
    chemistry: Optional[str] = None,
) -> SigProcConfig:
    if load_adapted_config_first and chemistry is None:
        raise ValueError(
            "chemistry must be provided if load_adapted_config_first is True"
        )
    if load_adapted_config_first and chemistry is not None:
        adapted_config_name = adapted_chemistry_specific_config_name(chemistry)
        adapted_config_dict = adapted_config_name_to_dict(adapted_config_name)
        warpdemux_config_dict = config_name_to_dict(config_name)
        adapted_config_dict.update(warpdemux_config_dict)
        return nested_config_from_dict(adapted_config_dict, SigProcConfig)
    else:
        with pkg_resources.path(config_files, f"{config_name}.toml") as config_path:
            return load_nested_config_from_file(config_path, SigProcConfig)


def config_name_to_dict(config_name: str) -> Union[dict, MutableMapping[str, Any]]:
    with pkg_resources.path(config_files, f"{config_name}.toml") as config_path:
        return toml.load(config_path)
