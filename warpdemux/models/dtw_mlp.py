"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import logging
from typing import Literal, Optional, Tuple, Union, overload

import numpy as np
import pandas as pd

from warpdemux.models.dtw_base import BaseDTWModel
from warpdemux.models.utils import predictions_to_df
from warpdemux.parallel_distances import distance_matrix_to


class DTW_MLP(BaseDTWModel):

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
        """set nproc to 1 to disable parallelization, nproc=None to use all available cores

        """

        if not self.is_trained:
            msg = "Model not trained yet."
            logging.error(msg)
            raise ValueError(msg)

        assert self._X is not None  # This informs the type checker
        assert self.model is not None  # This informs the type checker

        if X.ndim == 1:
            X = X.reshape(1, -1)

        if X.shape[1] != self._X.shape[1]:
            raise ValueError(
                "X must have the same shape in axis 1 as the consensus sequences "
                f" ({self._X.shape})."
            )

        D = distance_matrix_to(
            X,
            self._X,
            window=self.window,
            penalty=self.penalty,
            block_size=self.block_size if block_size is None else block_size,
            n_jobs=nproc,
            pbar=pbar,
            pbar_kwargs=pbar_kwargs,
        )

        y_prob = self.model.predict_proba(D)
        y_pred, conf = self.process_probs(y_prob)

        if return_df:
            if self.label_mapper is None:
                raise ValueError("Label mapper is not set.")
            return predictions_to_df(y_pred, y_prob, conf, self.label_mapper)

        return y_pred, y_prob

    def num_bcs(self) -> int:
        if self.n_classes is not None:
            return self.n_classes
        if self.label_mapper is not None:
            return len(self.label_mapper) - self.noise_class
        raise ValueError("No number of barcodes available.")