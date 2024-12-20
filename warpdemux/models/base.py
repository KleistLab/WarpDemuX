import logging
from abc import ABC, abstractmethod
from typing import Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from warpdemux.models.utils import confidence_margin


class BaseDTWModel(ABC):
    model: Optional[Union[Pipeline, SVC]] = None
    n_classes: Optional[int] = None
    _X: Optional[np.ndarray] = None
    label_mapper: Optional[Dict[int, int]] = None
    cal_dict: Optional[Dict[str, np.ndarray]] = None
    noise_class: bool = False

    @property
    def is_trained(self):
        return (
            self.model is not None
            # and self.model.fit_status_ == 0
            and self._X is not None
        )

    @property
    def num_bcs(self):
        if self.model is None:
            msg = "Model not trained yet."
            logging.error(msg)
            raise ValueError(msg)
        return self.model.classes_.size

    @abstractmethod
    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        nproc: int = -1,
        block_size: Optional[int] = None,
    ):
        """Train the model on the given data.

        Args:
            X: Training data
            y: Training labels (ordered ints, with highest label as outlier class)
            nproc: Number of processes (-1 for all cores, 1 to disable parallelization)
            block_size: Size of blocks for parallel processing
        """
        pass

    @abstractmethod
    def predict(
        self,
        X: np.ndarray,
        nproc: int = -1,
        block_size: Optional[int] = None,
        pbar: bool = False,
        pbar_kwargs: dict = {},
        return_df: bool = False,
        tam: Optional[str] = None,
    ) -> Union[Tuple[np.ndarray, np.ndarray], pd.DataFrame]:
        """Make predictions on new data.

        Args:
            X: Data to predict
            nproc: Number of processes
            block_size: Size of blocks for parallel processing
            pbar: Whether to show progress bar
            pbar_kwargs: Additional progress bar arguments
            return_df: Whether to return DataFrame instead of arrays
            tam: Target accuracy mode
        """
        pass

    def prob_to_pred(self, y_prob: np.ndarray, tam: Optional[str] = None) -> np.ndarray:
        if self.label_mapper is None:
            msg = "Label mapper not set."
            logging.error(msg)
            raise ValueError(msg)
        pred = np.array([self.label_mapper[i] for i in y_prob.argmax(axis=1)])
        if tam is not None:
            mask = self.get_tam_mask(y_prob, tam)
            pred[mask] = -1
        return pred

    def get_tam_mask(self, prob: np.ndarray, tam: str) -> np.ndarray:
        if self.cal_dict is None:
            msg = "Calibration dictionary not set."
            logging.error(msg)
            raise ValueError(msg)

        if tam not in self.cal_dict:
            msg = f"Target accuracy mode {tam} not in calibration dictionary."
            logging.error(msg)
            raise ValueError(msg)

        thresholds = self.cal_dict[tam]

        # TODO: handle noise class
        pred_idx = np.argmax(prob, axis=1)
        conf = confidence_margin(prob)
        mask = conf < thresholds[pred_idx]
        return mask

    def predictions_to_df(self, y_pred: np.ndarray, y_prob: np.ndarray) -> pd.DataFrame:
        if self.label_mapper is None:
            msg = "Label mapper not set."
            logging.error(msg)
            raise ValueError(msg)

        conf = confidence_margin(y_prob)

        return pd.DataFrame(
            {
                "predicted_barcode": y_pred,
                "confidence_score": conf.round(3),
                **{
                    f"p{self.label_mapper[i]:02d}": y_prob[:, i].round(4)
                    for i in range(y_prob.shape[1])
                },
            }
        )
