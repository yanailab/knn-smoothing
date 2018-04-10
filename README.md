## kNN-smoothing for high-throughput single-cell RNA-Seq data

This repository contains reference Python, R, and Matlab implementations of the kNN-smoothing and kNN-smoothing 2 algorithms ([Wagner et al., 2017](https://www.biorxiv.org/content/early/2018/04/09/217737)) for smoothing UMI-filtered single-cell RNA-Seq data.

### New! Version 2 of the algorithm released! (4/9/2018)

Version 2 is a major improvement over our original algorithm, and performs much better whenever the data contains distinct cell populations with very similar expression profiles. Version 2 completely replaces the original version. It takes two parameters (`k` and `d`). `k` is the number of neighbors to use for smoothing (same as in the original version), and `d` is the number of principal components used for determining the nearest neighbors in each smoothing step. For most applications, the default value of `d=10` works well. See our [preprint](https://www.biorxiv.org/content/early/2018/04/09/217737) for a discussion of how to choose `k` and `d`.

### Running kNN-smoothing from the command-line

Follow these instructions to run the Python implementation of kNN-smoothing from the command-line. This is the recommend method to run kNN-smoothing if you don't usually do your data analysis in Python, or if you prefer to work on the command-line.

1. Install dependencies

   Make sure you have Python 3 and the Python packages `scikit-learn`, `pandas`, and `click` installed. The easiest way to install Python 3 as well as these packages is to download and install [Anaconda](https://github.com/yanailab/CEL-Seq-pipeline/blob/133912cd4ceb20af0c67627ab883dfce8b9668df/sample_sheet_example.txt) (select the "Python 3.6 version").

2. Download the GitHub repository

   [Download](https://github.com/yanailab/knn-smoothing/archive/master.zip) this GitHub repository, and extract the contents into a folder.

3. Test running the script

   To run the script, change into the folder where you extracted the files, and run (on Linux/Mac):
    
   ``` bash
   python3 knn_smooth.py --help
   ```

   You should see the following output:

    ```
	Usage: knn_smooth.py [OPTIONS]

	Options:
	  -k INTEGER          The number of neighbors to use for smoothing.
	  -d INTEGER          The number of principal components used to identify
		              neighbors. Set to 0 in order to invoke old version of
		              kNN-smoothing (not recommended).  [default: 10]
	  -f, --fpath TEXT    The input UMI-count matrix.
	  -o, --saveto TEXT   The output matrix.
	  -s, --seed INTEGER  Seed for pseudo-random number generator.  [default: 0]
	  --sep TEXT          Separator used in input file. The output file will use
		              this separator as well.  [default: \t]
	  --help              Show this message and exit.
    ```

4. Make sure your expression matrix file is formatted correctly

   By default, the script expects your expression matrix to be stored as a tab-separated plain-text file, with gene labels contained in the first column, and cell labels contained in the first row (the top-left "cell" in the matrix can either be empty or contain the first cell label).

   If your file uses a separator other than the tab character, you must specify it by passing the `--sep` argument to the script. For example, if you're using comma-separated values (csv), pass `--sep ,`.  This will also affect the separator used in the output file.

5. Run smoothing!

   Let's say your (tab-separated) expression matrix file is called `expression.tsv`, and you saved it in the same directory as the "knn_smooth.py" script. Then, to run smoothing with `k=15`, you would use:

   ``` bash
   python3 knn_smooth.py --k 15 --fpath expression.tsv --saveto expression_smoothed.tsv
   ```

   This will produce a smoothed matrix called `expression_smoothed.tsv`.


### Example

  Run kNN-smoothing 2 from the command-line:
  ``` bash
  $ python3 knn_smooth.py -k 127 -d 5 -f smoothing_test_data1.tsv -o smoothing_test_results.tsv
  ```

  Output:
  ```
	Loading the data... done. (Took 3.1 s.)
	The expression matrix contains 19208 genes and 2000 cells.

	Performing kNN-smoothing 2 with k=127 and d=5...
	Step 1/7: Smooth using k=1
		PCA took 1.4 s.
		The fraction of variance explained by the top 5 PCs is 4.8 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 0.7 s.
	Step 2/7: Smooth using k=3
		PCA took 1.4 s.
		The fraction of variance explained by the top 5 PCs is 6.3 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 0.8 s.
	Step 3/7: Smooth using k=7
		PCA took 1.4 s.
		The fraction of variance explained by the top 5 PCs is 8.4 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 1.0 s.
	Step 4/7: Smooth using k=15
		PCA took 1.4 s.
		The fraction of variance explained by the top 5 PCs is 11.7 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 1.6 s.
	Step 5/7: Smooth using k=31
		PCA took 1.4 s.
		The fraction of variance explained by the top 5 PCs is 16.9 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 3.0 s.
	Step 6/7: Smooth using k=63
		PCA took 1.5 s.
		The fraction of variance explained by the top 5 PCs is 25.4 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 5.7 s.
	Step 7/7: Smooth using k=127
		PCA took 1.4 s.
		The fraction of variance explained by the top 5 PCs is 38.7 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 10.6 s.
	kNN-smoothing finished in 37.9 s.

	Writing results to "smoothing_test_results.tsv"... done. (Took 24.0 s.)
  ```
