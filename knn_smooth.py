# K-nearest neighbor smoothing for high-throughput scRNA-Seq data
# (Python 3 implementation, depends on scikit-learn.
#  The command-line version also depends on click and pandas.)

# Authors:
#   Florian Wagner <florian.wagner@nyu.edu>
#   Yun Yan <yun.yan@nyumc.org>
# Copyright (c) 2017, 2018 New York University

import time
import sys
from math import log, ceil
import hashlib

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.decomposition import PCA
import numpy as np


def _median_normalize(X):
    """Performs median-normalization.

    Parameters
    ----------
    X : numpy.ndarray
        A p-by-n expression matrix containing UMI counts for p genes and n
        cells.

    Returns
    -------
    numpy.ndarray
        A p-by-n expression matrix containing the normalized UMI counts.

    Notes
    -----
    We first determine the median total UMI count per cell, and then scale
    each expression profile so that its total UMI count equals that number.
    This normalization method was originally described as "Model I" in
    Grün et al., Nature Methods 2014).
    """
    num_transcripts = np.sum(X, axis=0)
    X_norm = (np.median(num_transcripts) / num_transcripts) * X
    return X_norm


def _freeman_tukey_transform(X):
    """Applies the Freeman-Tukey transformation, y = sqrt(x) + sqrt(x+1).
    
    Parameters
    ----------
    X : numpy.ndarray
        A p-by-n expression matrix containing UMI counts for p genes and n
        cells (usually after median-normalization).

    Returns
    -------
    numpy.ndarray
        A p-by-n expression matrix containing the Freeman-Tukey-transformed
        UMI counts.

    Notes
    -----
    The Freeman-Tukey transformation serves to stabilize the variance of
    Poisson-distributed random variables. For X ~ Pois(l) with l >= 1, Freeman
    and Tukey (1953) show that Var(X) = 1 (+- 6%).
    """
    return np.sqrt(X) + np.sqrt(X+1)


def _calculate_pc_scores(matrix, d, seed=0):
    """Projects the cells onto their first d principal components.

    Input
    -----
    X: `numpy.ndarray`
        A p-by-n expression matrix containing the UMI counts for p genes and n
        cells.

    Returns
    -------
    `numpy.ndarray`
        A d-by-n matrix containing the coordinates of n cells in d-dimensional
        principal component space.

    Notes
    -----
    We perform median-normalization and Freeman-Tukey-transformation to the UMI
    counts, before performing PCA. Median-normalization serves to counteract
    efficiency noise (Grün et al., 2014), whereas Freeman-Tukey transformation
    stabilizes the technical variance of the data. While PCA does not require
    homoskedastic data, variance-stabilization ensures that the increased
    technical variance of highly expressed genes does not result in the first
    PCs being biased towards highly expressed genes.

    We specify svd_solver='randomized', which invokes the randomized algorithm
    by Halko et al. (2009) to efficiently calculate the first d principal
    components. (We assume that d << min(p, n-1).)
    """
    # median-normalize
    tmatrix = _median_normalize(matrix)
    # Freeman-Tukey transform
    tmatrix = _freeman_tukey_transform(tmatrix)
    pca = PCA(n_components=d, svd_solver='randomized', random_state=seed)
    t0 = time.time()
    tmatrix = pca.fit_transform(tmatrix.T).T
    t1 = time.time()
    var_explained = np.cumsum(pca.explained_variance_ratio_)[-1]
    print('\tPCA took %.1f s.' % (t1-t0)); sys.stdout.flush()
    print('\tThe fraction of variance explained by the top %d PCs is %.1f %%.'
          % (d, 100*var_explained))

    return tmatrix


def _calculate_pairwise_distances(X, num_jobs=1):
    """Calculates the distances between all cells in X.
    
    Input: numpy.ndarray
        A d-by-n matrix containing the coordinates of n cells in d-dimensional
        space.
    Output: numpy.ndarray
        A n-by-n matrix containing all pairwise distances between the cells.

    Notes
    -----
    This uses the Euclidean metric.
    """
    D = pairwise_distances(X.T, n_jobs=num_jobs, metric='euclidean')
    return D


def knn_smoothing(X, k, d=10, dither=0.03, seed=0):
    """K-nearest neighbor smoothing for UMI-filtered single-cell RNA-Seq data.
    
    This function implements an improved version of the kNN-smoothing 2
    algorithm by Wagner et al.
    (https://www.biorxiv.org/content/early/2018/04/09/217737).

    Parameters
    ----------
    X : numpy.ndarray
        A p-by-n expression matrix containing UMI counts for p genes and n
        cells. Must contain floating point values, i.e. dtype=np.float64.
    k : int
        The number of neighbors to use for smoothing.
    d : int, optional
        The number of principal components to use for identifying neighbors.
        Default: 10.
    dither : float, optional
        Amount of dither to apply to the partially smoothed and PCA-transformed
        data in each step. Specified as the fraction of the range of the
        cell scores for each PC. Default: 0.03.
    seed : int, optional
        The seed for initializing the pseudo-random number generator used by
        the randomized PCA algorithm. This usually does not need to be changed.
        Default: 0.

    Returns
    -------
    numpy.ndarray
        A p-by-n expression matrix containing the smoothed expression values.
        The matrix is not normalized. Therefore, even though efficiency noise
        is usually dampened by the smoothing, median-normalization of the
        smoothed matrix is recommended.
    
    Raises
    ------
    ValueError
        If X does not contain floating point values.
        If k is invalid (k < 1, or k >= n).
        If d is invalid (d < 1 or d > # principal components).
    """
    
    np.random.seed(seed)

    if not (X.dtype == np.float64 or X.dtype == np.float32):
        raise ValueError('X must contain floating point values! '
                         'Try X = np.float64(X).')

    p, n = X.shape
    num_pcs = min(p, n-1)  # the number of principal components

    if k < 1 or k > n:
        raise ValueError('k must be between 1 and and %d.' % n)
    if d < 1 or d > num_pcs:
        raise ValueError('d must be between 1 and %d.' % num_pcs)

    print('Performing kNN-smoothing v2.1 with k=%d, d=%d, and dither=%.3f...'
          % (k, d, dither))
    sys.stdout.flush()

    t0_total = time.time()

    if k == 1:
        num_steps = 0
    else:
        num_steps = ceil(log(k)/log(2))
    
    S = X.copy()
    
    for t in range(1, num_steps+1):
        k_step = min(pow(2, t), k)
        print('Step %d/%d: Smooth using k=%d' % (t, num_steps, k_step))
        sys.stdout.flush()
        
        Y = _calculate_pc_scores(S, d, seed=seed)
        if dither > 0:
            for l in range(d):
                ptp = np.ptp(Y[l, :])
                dy = (np.random.rand(Y.shape[1])-0.5)*ptp*dither
                Y[l, :] = Y[l, :] + dy
            

        # determine cell-cell distances using smoothed matrix
        t0 = time.time()
        D = _calculate_pairwise_distances(Y)
        t1 = time.time()
        print('\tCalculating pair-wise distance matrix took %.1f s.' % (t1-t0))
        sys.stdout.flush()
        
        t0 = time.time()
        A = np.argsort(D, axis=1, kind='mergesort')
        for j in range(X.shape[1]):
            ind = A[j, :k_step]
            S[:, j] = np.sum(X[:, ind], axis=1)

        t1 = time.time()
        print('\tCalculating the smoothed expression matrix took %.1f s.'
              % (t1-t0))
        sys.stdout.flush()

    t1_total = time.time()
    print('kNN-smoothing finished in %.1f s.' % (t1_total-t0_total))
    sys.stdout.flush()

    return S


if __name__ == '__main__':
    import click
    import pandas as pd

    @click.command()
    @click.option('-k', type=int,
                  help='The number of neighbors to use for smoothing.')
    @click.option('-d', default=10, show_default=True,
                  help='The number of principal components used to identify '
                       'neighbors.')
    @click.option('--dither', default=0.03, show_default=True,
                  help='The amount of dither to apply to the partially '
                       'smoothed and PCA-transformed data in each step. '
                       'Specified as the faction of range of the scores of '
                       'each PC.')
    @click.option('-f', '--fpath', help='The input UMI-count matrix.')
    @click.option('-o', '--saveto', help='The output matrix.')
    @click.option('-s', '--seed', default=0, show_default=True,
                  help='Seed for pseudo-random number generator.')
    @click.option('--sep', default='\t', show_default=False,
                  help='Separator used in input file. The output file will '
                       'use this separator as well.  [default: \\t]')
    @click.option('--test', is_flag=True,
                  help='Test if results for test data are correct.')
    def main(k, d, dither, fpath, saveto, seed, sep, test):

        print('Loading the data...', end=' '); sys.stdout.flush()
        t0 = time.time()
        matrix = pd.read_csv(fpath, index_col=0, sep=sep).\
                astype(np.float64)
        t1 = time.time()
        print('done. (Took %.1f s.)' % (t1-t0)); sys.stdout.flush()
        p, n = matrix.shape
        print('The expression matrix contains %d genes and %d cells.' % (p, n))
        sys.stdout.flush()
        print()

        S = knn_smoothing(matrix.values, k, d=d, dither=dither, seed=seed)
        print()

        print('Writing results to "%s"...' % saveto, end=' ')
        sys.stdout.flush()
        t0 = time.time()
        matrix = pd.DataFrame(S, index=matrix.index, columns=matrix.columns)
        matrix.to_csv(saveto, sep=sep)
        t1 = time.time()
        print('done. (Took %.1f s.)' % (t1-t0))

        if test:
            with open(saveto, 'rb') as fh:
                h = str(hashlib.md5(fh.read()).hexdigest())
                if h == 'c8ee70f41b141b781041075e280661ff':
                    print('Test successful!!!')
                else:
                    raise ValueError('Output not correct!')

    main()
