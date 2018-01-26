# K-nearest neighbor smoothing for UMI-filtered scRNA-Seq data
# (Python 3 implementation, depends on scikit-learn and numpy.)

# Authors:
#   Florian Wagner <florian.wagner@nyu.edu>
#   Yun Yan <yun.yan@nyumc.org>
# Copyright (c) 2017 New York University

import time
import sys
from math import log, ceil

from sklearn.metrics.pairwise import pairwise_distances
import numpy as np

import click


def _freeman_tukey_transform(X):
    """Returns the Freeman-Tukey transformed data."""
    return np.sqrt(X) + np.sqrt(X+1)


def _calculate_pairwise_distances(X, num_jobs=1):
    """Calculates all pairwise distances for X.

    Performs median-normalization and variance stabilization, before
    calculating distances using the Euclidean metric."""
    # median-normalize
    num_transcripts = np.sum(X, axis=0)
    T = (np.median(num_transcripts) / num_transcripts) * X
    # stabilize variance
    F = _freeman_tukey_transform(T)
    # calculate distances
    D = pairwise_distances(F.T, n_jobs=num_jobs, metric='euclidean')
    return D


def knn_smoothing(X, k, num_jobs=1):
    """K-nearest neighbor smoothing.

    Parameters
    ----------
    X : np.ndarray (should of type float/np.float64)
        The UMI count matrix.
    k : int
        The number of nearest neighbors to use for smoothing.
    num_jobs: int
        The number of threads to use. See scikit-learn's
        documentation of the `pairwise_distances` function.
    """
    assert k < X.shape[1], 'Specified k should be smaller than #cells.'
    assert isinstance(X, np.ndarray), 'Input should be a numpy.ndarray.'

    num_powers = ceil(log(k+1)/log(2))
    S = X.copy()

    for p in range(1, num_powers+1):
        k_step = min(pow(2,p)-1, k)
        print('Step %d/%d: Smooth using k=%d'
              % (p, num_powers, k_step)); sys.stdout.flush()

        # determine cell-cell distances based on smoothed matrix
        t0 = time.time()
        D = _calculate_pairwise_distances(S, num_jobs=num_jobs)
        t1 = time.time()
        print('Calculating the pair-wise distances took %.1f s.'
              % (t1-t0)); sys.stdout.flush()

        # sort the distances and generate new smoothed matrix
        t0 = time.time()
        A = np.argsort(D, axis=1, kind='mergesort')
        for j in range(X.shape[1]):
            ind = A[j, :(k_step+1)]
            S[:, j] = np.sum(X[:, ind], axis=1)

        t1 = time.time()
        print('Calculating the smoothed expression matrix took %.1f s.'
              %(t1-t0)); sys.stdout.flush()

    return S


@click.command()
@click.option('--k', default=0, help='Number of K.')
@click.option('--fpath', help='Input UMI-count matrix.')
@click.option('--saveto', help='Output smoothed UMI-count matrix.')
@click.option('--sep', help='File sep when reading fpath.', default='\t')
def main(k, fpath, saveto, sep):
    import pandas as pd

    expr = pd.read_csv(fpath, index_col=0, sep=sep).astype(np.float64)
    print('Performing KNN Smoothing with k={}'.format(k)); sys.stdout.flush()
    genes = expr.index
    cells = expr.columns
    expr = knn_smoothing(X=expr.values, k=k)
    expr = pd.DataFrame(expr, index=genes, columns=cells)
    expr.to_csv(saveto, sep=sep)


if __name__ == '__main__':
    main()

