import importlib.resources as pkg_resources

import toml

from warpdemux.models import model_files


def available_models():
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)
    return list(model_config.keys())
