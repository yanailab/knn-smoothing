# K-nearest neighbor smoothing for high-throughput scRNA-Seq data
# (Python 3 implementation, depends on scikit-learn.
#  The command-line version also depends on click and pandas.)

# Authors:
#   Florian Wagner <florian.wagner@nyu.edu>
#   Yun Yan <yun.yan@nyumc.org>
# Copyright (c) 2017, 2018 New York University

import click
import pandas as pd
import hashlib
from .knn_smooth import knn_smoothing

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
    print('Loading the data...', end=' ');
    sys.stdout.flush()
    t0 = time.time()
    matrix = pd.read_csv(fpath, index_col=0, sep=sep). \
        astype(np.float64)
    t1 = time.time()
    print('done. (Took %.1f s.)' % (t1 - t0));
    sys.stdout.flush()
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
    print('done. (Took %.1f s.)' % (t1 - t0))

    if test:
        with open(saveto, 'rb') as fh:
            h = str(hashlib.md5(fh.read()).hexdigest())
            if h == 'c8ee70f41b141b781041075e280661ff':
                print('Test successful!!!')
            else:
                raise ValueError('Output not correct!')

if __name__ == '__main__':
    main()
