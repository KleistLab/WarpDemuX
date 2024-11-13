"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import importlib.resources as pkg_resources
import logging
from typing import Any, MutableMapping, Optional, Union

import toml
from adapted.config.base import load_nested_config_from_file, nested_config_from_dict
from adapted.config.sig_proc import (
    chemistry_specific_config_name as adapted_chemistry_specific_config_name,
)
from adapted.config.sig_proc import config_name_to_dict as adapted_config_name_to_dict

from warpdemux import __version__
from warpdemux.config import config_files
from warpdemux.config.sig_proc import SigProcConfig
from warpdemux.models import model_files


# TODO: update model naming to WDX[xx]_rna00[x]_xxbps@v[x].[x].[x]
def get_model_spc_config(model_name: str) -> SigProcConfig:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    sqk = model_config["SQK"]  # sequencing kit, RNA002 or RNA004
    adapted_config_name = adapted_chemistry_specific_config_name(sqk)
    adapted_config_dict = adapted_config_name_to_dict(adapted_config_name)
    warpdemux_config_dict = config_name_to_dict(model_config["spc"])

    common_keys = set(adapted_config_dict.keys()) & set(warpdemux_config_dict.keys())
    # config dicts are nested, so we need to update nested dictionaries, rather than just merging the top level
    for key in common_keys:
        adapted_config_dict[key].update(
            warpdemux_config_dict[key]
        )  # merge the two dictionaries, warpdemux_config_dict overwrites adapted_config_dict
    keys_only_in_warpdemux = set(warpdemux_config_dict.keys()) - set(
        adapted_config_dict.keys()
    )
    for key in keys_only_in_warpdemux:
        adapted_config_dict[key] = warpdemux_config_dict[key]
    return nested_config_from_dict(adapted_config_dict, SigProcConfig)


def get_model_spc_live_config(model_name: str) -> SigProcConfig:
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)[model_name]
    return get_config(
        model_config["spc_live"],
        load_adapted_config_first=True,
        chemistry="rna002" if "rna002" in model_name else "rna004",
    )


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
        msg = f"Unknown chemistry: {chemistry}"
        logging.error(msg)
        raise ValueError(msg)
    if version is None:
        version = __version__

    if load_adapted_config_first and chemistry is None:
        msg = "chemistry must be provided if load_adapted_config_first is True"
        logging.error(msg)
        raise ValueError(msg)

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
        msg = "chemistry must be provided if load_adapted_config_first is True"
        logging.error(msg)
        raise ValueError(msg)

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
