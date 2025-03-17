import logging
from typing import Dict, Optional, Tuple, Union

import catboost as cb
import numpy as np
import pandas as pd

from warpdemux.models.fpt_base import BaseFptModel
from warpdemux.models.utils import predictions_to_df


class Fpt_Boost(BaseFptModel):
    
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
        if not self.is_trained:
            msg = "Model not trained."
            logging.error(msg)
            raise ValueError(msg)
        
        if not self.label_mapper:
            msg = "Label mapper not set."
            logging.error(msg)
            raise ValueError(msg)

        # Get raw probabilities from model
        y_prob = self.model.predict_proba(X, thread_count=1)
        y_pred, conf = self.process_probs(y_prob)
        
        if return_df:
            return predictions_to_df(y_pred, y_prob, conf, self.label_mapper)
        return y_pred, conf
