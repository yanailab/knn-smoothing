## kNN-smoothing for single-cell RNA-Seq data

This repository contains reference Python, R, and Matlab implementations of the kNN-smoothing and kNN-smoothing 2 algorithms ([Wagner et al., 2017](https://www.biorxiv.org/content/early/2018/01/24/217737)) for UMI-filtered single-cell RNA-Seq data.

### New! Version 2 of the algorithm released! (4/9/2018)



### Overview of the different implementations

Even though all implementations should produce the same results, the Python implementation (`knn_smooth.py`) currently runs much faster than the R implementation. This is despite our efforts to optimize the performance of the R implementation (`knn_smooth.R`).

If you have a large dataset and/or want to apply smoothing with a large `k`, we therefore recommmend that you either call the Python function `knn_smoothing()` from `knn_smooth.py`, or that you run the Python implementation of kNN-smoothing from the command-line (see below).

### Running kNN-smoothing from the command-line

Follow these instructions to run the Python implementation of kNN-smoothing from the command-line. Since the Python implementation is currently much faster than the R implementation (see above), this is the recommend method to run kNN-smoothing if you don't usually do your data analysis in Python, or if you prefer to work on the command-line.

1. Install dependencies

   Make sure you have Python 3 and the Python packages `pandas`,  `scikit-learn` and `click` installed. The easiest way to install Python 3 as well as these packages is to download and install [Anaconda](https://github.com/yanailab/CEL-Seq-pipeline/blob/133912cd4ceb20af0c67627ab883dfce8b9668df/sample_sheet_example.txt) (select the "Python 3.6 version").

2. Download the GitHub repository

   [Download](https://github.com/yanailab/knn-smoothing/archive/master.zip) this GitHub repository, and extract the contents into a folder.

3. Test running the script

   To run the script, change into the folder where you extracted the files, and run (on Linux/Mac):
    
   ``` bash
   python3 knn_smooth.py --help
   ```

   You should see the following output:

            Usage: knn_smooth.py [OPTIONS]
            
            Options:
            --k INTEGER    Number of K.
            --d INTEGER
            --fpath TEXT   Input UMI-count matrix.
            --saveto TEXT  Output smoothed UMI-count matrix.
            --sep TEXT     File sep when reading fpath.
            --help         Show this message and exit.


4. Make sure your expression matrix file is formatted correctly

   By default, the script expects your expression matrix to be stored as a tab-delimited plain-text file, with gene labels contained in the first column, and cell labels contained in the first row (the top-left "cell" in the matrix can either be empty or contain the first cell label).

   You can change the separator by passing the `--sep` argument to the script. For example, if you're using comma-separated values (csv), pass `--sep ,`. This will also affect the separator used in the output file.

5. Run smoothing!

   Let's say your (tab-delimited) expression matrix file is called `expression.tsv`, and you saved it in the same directory as the "knn_smooth.py" script. Then, to run smoothing with `k=15`, you would use:

   ``` bash
   python3 knn_smooth.py --k 15 --fpath expression.tsv --saveto expression_smoothed.tsv
   ```

   This will produce a smoothed matrix called `expression_smoothed.tsv`.


### Sample output

  ``` bash
  $ python3 knn_smooth.py -k 127 -d 5 -f smoothing_test_data1.tsv -o smoothing_test_results.tsv
  ```

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
