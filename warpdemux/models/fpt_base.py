import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd

from warpdemux.models.utils import process_probs as _process_probs


class BaseFptModel(ABC):
    def __init__(self, model: Any, n_classes: int, label_mapper: Dict[int, int], noise_class: bool, thresholds: np.ndarray):
        
        """Fingerprinting-based model.

        Args:
            model: Pre-trained model
            n_classes: Number of classes (excluding noise class)
            label_mapper: Mapping from model output indices to barcode indices
            noise_class: Whether model includes a noise class
            thresholds: calibration thresholds for individual classes
        """
        
        self.model: Any = model
        self.n_classes: int = n_classes
        self.label_mapper: Dict[int, int] = label_mapper
        self.noise_class: bool = noise_class
        self.thresholds: np.ndarray = thresholds
        
    @property
    def is_trained(self):
        return (
            self.model is not None
        )
        

    @property
    def num_bcs(self) -> int:
        if self.n_classes is not None:
            return self.n_classes
        elif self.label_mapper is not None:
            return len(self.label_mapper) - self.noise_class
        else:
            msg = "Label mapper or n_classes not set."
            logging.error(msg)
            raise ValueError(msg)
        

    @abstractmethod
    def predict(
        self,
        X: np.ndarray,
        return_df: bool = False,
        **kwargs
    ) -> Union[Tuple[np.ndarray, np.ndarray], pd.DataFrame]:
        """Make predictions on new data.

        Args:
            X: Data to predict
            return_df: Whether to return DataFrame instead of arrays
            
            **kwargs: Additional arguments are accepted but ignored

        Returns:
            If return_df is False:
                Tuple of (predictions, confidences)
            If return_df is True:
                DataFrame with columns 'pred' and 'conf'
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