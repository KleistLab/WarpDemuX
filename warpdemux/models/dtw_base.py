import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC

from warpdemux.models.utils import process_probs as _process_probs


class BaseDTWModel(ABC):
    def __init__(self, model: Any, _X: np.ndarray, n_classes: int, label_mapper: Dict[int, int], noise_class: bool, thresholds: np.ndarray, block_size: int, window: int, penalty: float):
        self.block_size: int = block_size
        self.window: int = window
        self.penalty: float = penalty

        self.model: Union[Pipeline, SVC] = model
        self._X: np.ndarray = _X
        self.n_classes: int = n_classes
        self.label_mapper: Dict[int, int] = label_mapper
        self.noise_class: bool = noise_class

        self.thresholds: np.ndarray = thresholds


    @property
    def is_trained(self):
        return (
            self.model is not None
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
    def predict(
        self,
        X: np.ndarray,
        nproc: int = -1,
        block_size: Optional[int] = None,
        pbar: bool = False,
        pbar_kwargs: dict = {},
        return_df: bool = False,
    ) -> Union[Tuple[np.ndarray, np.ndarray], pd.DataFrame]:
        """Make predictions on new data.

        Args:
            X: Data to predict
            nproc: Number of processes
            block_size: Size of blocks for parallel processing
            pbar: Whether to show progress bar
            pbar_kwargs: Additional progress bar arguments
            return_df: Whether to return DataFrame instead of arrays
        """
        pass


    def get_thresholds(self) -> np.ndarray:
        return self.thresholds
    
    def process_probs(self, y_prob: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if self.label_mapper is None:
            msg = "Label mapper not set."
            logging.error(msg)
            raise ValueError(msg)
        
        thresholds = self.get_thresholds()
        pred, conf = _process_probs(y_prob, self.label_mapper, thresholds)
        return pred, conf
