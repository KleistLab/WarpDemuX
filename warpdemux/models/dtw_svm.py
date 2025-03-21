"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import logging
from typing import Literal, Optional, Tuple, Union, overload

import numpy as np
import pandas as pd
from sklearn import svm

from warpdemux.models.dtw_base import BaseDTWModel
from warpdemux.models.utils import predictions_to_df
from warpdemux.parallel_distances import distance_matrix_to


def pdist_kernel(pdist: np.ndarray, gamma: float = 1, pwr_dist: int = 1) -> np.ndarray:
    return np.exp(-gamma * np.power(pdist, pwr_dist))


class DTW_SVM(BaseDTWModel):
    def __init__(self, gamma: float = 1, pwr_dist: int = 1, C: float = 1, **kwargs):
        super().__init__(**kwargs)
        self.gamma: float = gamma
        self.pwr_dist: int = pwr_dist
        self.C: float = C

    @overload
    def predict(
        self,
        X: np.ndarray,
        nproc: int = -1,
        block_size: Optional[int] = None,
        pbar: bool = False,
        pbar_kwargs: dict = {},
        return_df: Literal[True] = True,
    ) -> pd.DataFrame: ...

    @overload
    def predict(
        self,
        X: np.ndarray,
        nproc: int = -1,
        block_size: Optional[int] = None,
        pbar: bool = False,
        pbar_kwargs: dict = {},
        return_df: Literal[False] = False,
    ) -> Tuple[np.ndarray, np.ndarray]: ...

    def predict(
        self,
        X: np.ndarray,
        nproc: int = -1,
        block_size: Optional[int] = None,
        pbar: bool = False,
        pbar_kwargs: dict = {},
        return_df: bool = False,
    ) -> Union[Tuple[np.ndarray, np.ndarray], pd.DataFrame]:
        """set nproc to 1 to disable parallelization, nproc=None to use all available cores"""

        if not self.is_trained:
            msg = "Model not trained yet."
            logging.error(msg)
            raise ValueError(msg)

        if X.ndim == 1:
            X = X.reshape(1, -1)

        if X.shape[1] != self._X.shape[1]:
            raise ValueError(
                "X must have the same number of columns as the training data "
                f" ({self._X.shape[1]})."
            )

        dtw_pdist = distance_matrix_to(
            X,
            self._X,
            window=self.window,
            penalty=self.penalty,
            block_size=self.block_size if block_size is None else block_size,
            n_jobs=nproc,
            pbar=pbar,
            pbar_kwargs=pbar_kwargs,
        )

        K = pdist_kernel(dtw_pdist, gamma=self.gamma, pwr_dist=self.pwr_dist)

        y_prob = self.model.predict_proba(K)
        y_pred, conf = self.process_probs(y_prob)

        if return_df:
            return predictions_to_df(y_pred, y_prob, conf, self.label_mapper)

        return y_pred, y_prob
