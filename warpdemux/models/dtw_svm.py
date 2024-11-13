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

from warpdemux.parallel_distances import distance_matrix_to


def pdist_kernel(pdist: np.ndarray, gamma: float = 1, pwr_dist: int = 1) -> np.ndarray:
    return np.exp(-np.power(pdist, pwr_dist) / gamma)


def confidence_margin(npa: np.ndarray) -> np.ndarray:
    sorted = np.sort(npa, axis=1)[:, ::-1]  # return sort in reverse, i.e. descending
    d = sorted[:, 0] - sorted[:, 1]
    return d


class DTW_SVM_Model:
    def __init__(
        self,
        C: float = 10.0,
        gamma: float = 1.0,
        pwr_dist: int = 1,
        block_size: int = 1000,
        window: int = 15,
        penalty: float = 0.1,
        noise_class: bool = True,
    ):
        self.C: float = C
        self.gamma: float = gamma
        self.model: Optional[svm.SVC] = None
        self.pwr_dist: int = pwr_dist
        self.block_size: int = block_size
        self.window: int = window
        self.penalty: float = penalty
        self.noise_class: bool = noise_class

        self.n_classes: Optional[int] = None
        self.X_train: Optional[np.ndarray] = None

        self.label_mapper: Optional[dict[int, int]] = None

    @property
    def is_trained(self):
        return (
            self.model is not None
            and self.model.fit_status_ == 0
            and self.X_train is not None
        )

    @property
    def num_bcs(self):
        if self.model is None:
            msg = "Model not trained yet."
            logging.error(msg)
            raise ValueError(msg)

        return self.model.classes_.size

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        nproc: int = -1,
        block_size: Optional[int] = None,
    ):
        """set nproc to 1 to disable parallelization, nproc=None to use all available cores
        D_stats should be dict of np.ndarray with keys "D_med" and "D_mad"

        When also fitting an outlier class, this class should have the highest label value.
        """

        # TODO: document expected input format for y: ordered ints, with the highest class label as the outlier class
        self.label_mapper = {i: u for i, u in enumerate(np.unique(y))}

        if self.noise_class:
            self.label_mapper[np.unique(y).size - 1] = -1

        dtw_pdist = distance_matrix_to(
            X,
            X,
            window=self.window,
            penalty=self.penalty,
            block_size=self.block_size if block_size is None else block_size,
            n_jobs=nproc,
        )

        K = pdist_kernel(dtw_pdist, gamma=self.gamma, pwr_dist=self.pwr_dist)

        self.model = svm.SVC(kernel="precomputed", probability=True, C=self.C)
        self.model.fit(K, y)

        # reduce model size
        support_indices = self.model.support_
        self.model.n_features_in_ = self.model.support_.size
        self.model.shape_fit_ = (self.model.support_.size, self.model.support_.size)
        self.model.support_ = np.arange(self.model.support_.size, dtype=np.int32)

        self.X_train = X[support_indices]

    # TODO: add confidence threshold with 0 class for unclassified reads
    def prob_to_pred(self, y_prob: np.ndarray) -> np.ndarray:
        if self.label_mapper is None:
            msg = "Label mapper not set."
            logging.error(msg)
            raise ValueError(msg)

        return np.array([self.label_mapper[i] for i in y_prob.argmax(axis=1)])

    def predictions_to_df(self, y_pred: np.ndarray, y_prob: np.ndarray) -> pd.DataFrame:
        if self.label_mapper is None:
            msg = "Label mapper not set."
            logging.error(msg)
            raise ValueError(msg)
        return pd.DataFrame(
            {
                "predicted_barcode": y_pred,  # -1 is the "no barcode (noise)" class
                "confidence_score": confidence_margin(y_prob).round(3),
                **{
                    f"p{self.label_mapper[i]:02d}": y_prob[:, i].round(4)
                    for i in range(y_prob.shape[1])
                },
            }
        )

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

        assert self.X_train is not None  # This informs the type checker
        assert self.model is not None  # This informs the type checker

        if X.ndim == 1:
            X = X.reshape(1, -1)

        if X.shape[1] != self.X_train.shape[1]:
            raise ValueError(
                "X must have the same number of columns as the training data "
                f" ({self.X_train.shape[1]})."
            )

        dtw_pdist = distance_matrix_to(
            X,
            self.X_train,
            window=self.window,
            penalty=self.penalty,
            block_size=self.block_size if block_size is None else block_size,
            n_jobs=nproc,
            pbar=pbar,
            pbar_kwargs=pbar_kwargs,
        )

        K = pdist_kernel(dtw_pdist, gamma=self.gamma, pwr_dist=self.pwr_dist)

        y_prob = self.model.predict_proba(K)
        y_pred = self.prob_to_pred(y_prob)

        if return_df:
            return self.predictions_to_df(y_pred, y_prob)

        return y_pred, y_prob
