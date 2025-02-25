import importlib.resources as pkg_resources

import numpy as np
import toml
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_sample_weight

from warpdemux.models import model_files


def available_models():
    with pkg_resources.path(model_files, "config.toml") as config_path:
        model_config = toml.load(config_path)
    return list(model_config.keys())


def confidence_margin(npa: np.ndarray) -> np.ndarray:
    sorted = np.sort(npa, axis=1)[:, ::-1]  # return sort in reverse, i.e. descending
    d = sorted[:, 0] - sorted[:, 1]
    return d


class WeightedStandardScaler(StandardScaler):
    def fit(self, X, y=None):
        sample_weights = compute_sample_weight("balanced", y)

        # Compute weighted mean and std
        self.mean_ = np.average(X, axis=0, weights=sample_weights)
        self.scale_ = np.sqrt(
            np.average((X - self.mean_) ** 2, axis=0, weights=sample_weights)
        )
        return self
