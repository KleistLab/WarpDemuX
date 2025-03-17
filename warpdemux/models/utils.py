import importlib.resources as pkg_resources
from typing import Optional, Tuple

import numpy as np
import pandas as pd
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

def predictions_to_df(y_pred: np.ndarray, y_prob: np.ndarray, conf: np.ndarray, label_mapper: dict) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "predicted_barcode": y_pred,
            "confidence_score": conf.round(3),
            **{f"p{label_mapper[i]:02d}": y_prob[:, i].round(4) for i in range(y_prob.shape[1])},
        }
    )

def process_probs(y_prob: np.ndarray, label_mapper: dict, thresholds: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Process probabilities to predictions and confidence scores.

    Args:
        y_prob: Probabilities of each class
        label_mapper: Label mapper
        thresholds: Thresholds for each class, confidence scores below these are filtered out
    """
    pred_idx = np.argmax(y_prob, axis=1)
    pred = np.array([label_mapper[i] for i in pred_idx])
    
    conf = confidence_margin(y_prob)
    
    if thresholds is not None:
        mask = conf < thresholds[pred_idx]
        pred[mask] = -1
    return pred, conf
