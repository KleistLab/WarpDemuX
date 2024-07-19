from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional, Tuple

import numpy as np
from dtaidistance import dtw
from tqdm import tqdm


def compute_block_distance(
    block_indices: Tuple[np.ndarray, np.ndarray],
    X,
    window=None,
    penalty=None,
    **kwargs,
):
    # Extract the blocks from indices
    i, j = block_indices

    block_dtw_dist = dtw.distance_matrix(
        np.vstack([X[i], X[j]]),
        block=((0, i.size), (i.size, i.size + j.size)),
        parallel=False,
        use_c=True,
        only_triu=True,
        window=window,
        penalty=penalty,
        **kwargs,
    )[: i.size, i.size :].astype(np.float32)

    return (i, j, block_dtw_dist)


def distance_matrix_to(
    X,
    Y,
    window: Optional[int] = None,
    penalty: Optional[float] = None,
    block_size: Optional[int] = None,
    n_jobs: int = -1,
    pbar: bool = False,
    pbar_kwargs: dict = {},
):
    if n_jobs == 1:
        return dtw.distance_matrix(
            np.vstack([X, Y]),
            parallel=False,
            use_c=True,
            only_triu=True,
            window=window,
            penalty=penalty,
            block=((0, X.shape[0]), (X.shape[0], X.shape[0] + Y.shape[0])),
        )[: X.shape[0], X.shape[0] :].astype(np.float32)

    else:
        if block_size is None:
            raise ValueError("block_size must be specified when using parallel.")

        return parallel_distance_matrix_to(
            X,
            Y,
            block_size=block_size,
            n_jobs=n_jobs,
            window=window,
            penalty=penalty,
            pbar=pbar,
            pbar_kwargs=pbar_kwargs,
        )


def parallel_distance_matrix_to(
    X,
    Y,
    block_size: int = 1000,
    n_jobs: int = 6,
    window: Optional[int] = None,
    penalty: Optional[float] = None,
    pbar: bool = False,
    pbar_kwargs: dict = {},
    **kwargs,
):
    """
    Computes the pairwise distance matrix in parallel between two datasets, X and Y.

    This function concatenates X and Y vertically and then computes the pairwise
    distance matrix between all rows of X and all rows of Y using blocks.

    Parameters:
    -----------
    X : array-like
        First dataset, where each row is an instance. Shape (n_samples_X, n_features).

    Y : array-like
        Second dataset, where each row is an instance. Shape (n_samples_Y, n_features).

    block_size : int, optional (default=1000)
        Size of the block for pairwise computation. Larger block sizes can
        leverage efficient matrix operations but might consume more memory.

    n_jobs : int, optional (default=6)
        Number of parallel jobs.

    Returns:
    --------
    numpy.ndarray
        The pairwise distance matrix where the rows correspond to instances in X
        and columns correspond to instances in Y. Shape (n_samples_X, n_samples_Y).
    """

    return parallel_distance_matrix(
        np.vstack([X, Y]),
        block_size=block_size,
        n_jobs=n_jobs,
        subset=((0, X.shape[0]), (X.shape[0], X.shape[0] + Y.shape[0])),
        window=window,
        penalty=penalty,
        pbar=pbar,
        pbar_kwargs=pbar_kwargs,
        **kwargs,
    )


def parallel_distance_matrix(
    X,
    block_size: int = 1000,
    n_jobs: int = 6,
    subset: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None,
    window: Optional[int] = None,
    penalty: Optional[float] = None,
    pbar: bool = False,
    pbar_kwargs: dict = {},
    **kwargs,
):
    """
    Computes the pairwise distance matrix in parallel using block-based computation.

    :param X: Input data matrix.
    :param block_size: Size of the block for pairwise computation.
    :param n_jobs: Number of parallel jobs. -1 means use all available CPUs.
    :param subset: Tuple containing two range-tuples to define the subset.
                   E.g., ((x1, x2), (x3, x4)) computes the pairwise distances
                   between rows x1:x2 and rows x3:x4.
    :return: Pairwise distance matrix.
    """
    if subset:
        (r1_start, r1_end), (r2_start, r2_end) = subset
    else:
        r1_start, r1_end, r2_start, r2_end = 0, X.shape[0], 0, X.shape[0]

    blocks = [
        (
            np.arange(i, min(i + block_size, r1_end)),
            np.arange(j, min(j + block_size, r2_end)),
        )
        for i in np.arange(r1_start, r1_end, block_size)
        for j in np.arange(r2_start, r2_end, block_size)
    ]

    # Compute pairwise distances in parallel
    result_matrix = np.zeros((r1_end - r1_start, r2_end - r2_start), dtype=np.float32)

    with ProcessPoolExecutor(max_workers=n_jobs if n_jobs != -1 else None) as executor:
        futures = [
            executor.submit(compute_block_distance, b, X, window, penalty, **kwargs)
            for b in blocks
        ]

        if pbar:
            for future in tqdm(
                as_completed(futures),
                total=len(blocks),
                desc="Computing Kernel Matrix",
                **pbar_kwargs if pbar_kwargs else {},
            ):
                i, j, distances = future.result()
                result_matrix[np.ix_(i - r1_start, j - r2_start)] = distances
        else:
            for future in as_completed(futures):
                i, j, distances = future.result()
                result_matrix[np.ix_(i - r1_start, j - r2_start)] = distances

    return result_matrix
