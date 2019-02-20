# PEGASUS-WINGS
Ward clustering using Inter Node distances of Gene Score arrays to define signficance (WINGS)

Authors: Melissa R. McGuirl, Samuel Pattillo Smith, Bjorn Sandstede, Sohini Ramachandran

# Description: 
This repository contains MATLAB code for WINGS, an algorithim which clusters phenotypes using gene score as the input feature vectors with a thresholded Ward hierarchical clustering algorithm. The identified clusters are ranked by significance. 
 
### Usage:

```
matlab 
WINGS(gene_scores, outdir, scale, plotOn, saveOn, cluser_alg, cluster_dist, sz_thresh)

```

# Inputs:
1) gene_scores is a .txt or .csv file containing a matrix of PEGASUS gene score feature vectors of size M + 1 X N + 1,
where M = number of gene, N = number of phenotypes with the first row containing all the gene names/labels and
the first column containing all of the phenotype names/labels.

2) outdir is the directory where the user wants to save the output images and files (optional, default = './').

3) scale is the scale on which the scores will be clustered. Options: "raw" (non-transformed) or "log10"
(-log10transformed). NOTE: WE ASSUME THE INPUTTED GENE SCORE MATRIX IS NOT TRANSFORMED. (optional, default = 'log10').

4) plotOn is a plot indicator variable. If plotOn = 1, this function will plot and save the dendrogram and sorted
branch length figure (optional, default = 1).

5) saveOn is a save indicator variable. If saveOn = 1, this function will save text files with (1) the estimated
significant clusters with corresponding branch length threshold and (2) all clusters ranked by branch length.
(optional, default = 1).

6) cluster_alg is the clustering method that describes how to measure the distance between clusters in the hierarchical clustering algorithm. Options include "average", "centroid", "complete", "median", "single" "ward", and
"weighted." See "doc linkage" (method) for more details. (optional, default = "ward").

7) cluster_dist is the distance used to measure differences between phenotypes. Some example options include "euclidean", "squaredeuclidean", "seuclidean", "cityblock", "cosine", "hammning", "jaccard", or any other metric that is
accepted by pdist. (optional, default = "euclidean").

8) sz_thresh is the size threshold used to define a cluster. All clusted of size greater than sz_thresh will not be recognized as "true clusters." (optional, default = ceil(N/3)).

# Outputs: 
Dendrogram and branch length figures (if plotOn = 1) and text files of the estimated
significant clusters with corresponding branch length threshold and (2) all clusters ranked by branch length.
(if saveOn = 1).

# Sample Run: 
```
cd src
WINGS('../examples/sim_example_small.csv', '../outputs', 'raw', 1, 1, 'ward', 'euclidean', 4)
```

# Data and Outputs: 
Sample inputs are contained in the 'examples' folder and sample outputs and figures are contained in the 'outputs' folder. 'sim_example.csv' and 'sim_example_small.csv' are toy simulations, whereas 'UKBiobank_small_raw.txt', 'UKBiobank_expanded_raw.txt', and 'UKBiobank_expanded_log10.txt' correspond to the empirical PEGASUS scores for the UKBiobank data (as described in the paper). The sample results in 'outputs' correspond to the outputs from running the sample run above. 

# Requirements:
MATLAB with the Statistics and Machine Learning Toolbox

# References:
This work is based upon the paper: [add citation] 

# Contact:
Corresponding Author for Code: Melissa R. McGuirl, melissa_mcguirl@brown.edu
