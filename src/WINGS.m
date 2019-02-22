%% function: cluster_phenos(gene_scores, outdir, scale, plotOn, saveOn, cluster_alg, cluster_dist, sz_thresh)
%%
%% Description: This is a function that takes as input a matrix of PEGASUS gene score feature vectors for a collection of
%% phenotypes and runs a threshold hierarchical clustering algorithm to identify and rank significant phenotype clusters
%%
%% Inputs:
%% 1) gene_scores is a .txt or .csv file containing a matrix of PEGASUS gene score feature vectors of size M + 1 X N + 1, 
%% where M = number of genes, N = number of phenotypes with the first column containing all the gene names/labels and the first 
%% row containing all of the phenotype names/labels.
%%
%% 2) outdir is the directory where the user wants to save the output images and files (optional, default = './').
%%
%% 3) scale is the scale on which the scores will be clustered. Options: "raw" (non-transformed) or "log10"
%% (-log10transformed). NOTE: WE ASSUME THE INPUTTED GENE SCORE MATRIX IS NOT TRANSFORMED. (optional, default = 'log10').
%%
%% 4) plotOn is a plot indicator variable. If plotOn = 1, this function will plot and save the dendrogram and sorted
%% branch length figure (optional, default = 1).
%%
%% 5) saveOn is a save indicator variable. If saveOn = 1, this function will save text files with (1) the estimated
%% significant clusters with corresponding branch length threshold and (2) all clusters ranked by branch length.
%% (optional, default = 1).
%%
%% 6) cluster_alg is the clustering method that describes how to measure the distance between clusters in the hierarchical
%% clustering algorithm. Options include "average", "centroid", "complete", "median", "single" "ward", and
%% "weighted." See "doc linkage" (method) for more details. (optional, default = "ward").
%%
%% 7) cluster_dist is the distance used to measure differences between phenotypes. Some example options include
%% "euclidean", "squaredeuclidean", "seuclidean", "cityblock", "cosine", "hammning", "jaccard", or any other metric that is
%% accepted by pdist. (optional, default = "euclidean").
%%
%% 8) sz_thresh is the size threshold used to define a cluster. All clusted of size greater than sz_thresh will not be
%% recognized as "true clusters." (optional, default = ceil(N/3)).
%%
%% Outputs: Dendrogram and branch length figures (if plotOn = 1) and text files of the estimated
%% significant clusters with corresponding branch length threshold and (2) all clusters ranked by branch length.
%% (if saveOn = 1).
%%
%% Sample Run: cluster_phenos(path_to_gene_scores, './', 'log10', 1, 1, 'ward', 'euclidean', 9)
%%
%% 
%% Corresponding Author: Melissa R. McGuirl, melissa_mcguirl@brown.edu

function [] = cluster_phenos(gene_scores, outdir, scale, plotOn, saveOn, cluster_alg, cluster_dist, sz_thresh)

% first check if optional inputs given, otherwise set to default value.
if ~exist('gene_scores', 'var')
    fprintf('Error: minimum one input variable required. Sample run: cluster_phenos(path_to_gene_score_file)')
    return
end

if ~exist('outdir', 'var')
    scale = './';
end

if ~exist('scale', 'var')
    scale = 'log10';
end

if ~exist('plotOn', 'var')
    plotOn = 1;
end

if ~exist('saveOn', 'var')
    saveOn = 1;
end

if ~exist('cluster_alg', 'var')
    cluster_alg = 'ward';
end

if ~exist('cluster_dist', 'var')
    cluster_dist = 'euclidean';
end

% load in data and extract gene score vector and labels
data = importdata(gene_scores);
% pegasus_scores = data(2:end, 2:end);

pegasus_scores = data.data;
pheno_labels = data.textdata(1,2:end);

% pheno_labels = cell(1,size(pegasus_scores,2));
% for k = 1:size(pegasus_scores,2)
%     pheno_labels{k} = char(string((k-1)));
% end

%gene_labels = data.textdata(2:end,1);  % uncomment for gene labels
N_phenos = size(pheno_labels,2);
%N_genes = size(gene_labels,1);  % uncomment for number of genes

% get cluster size threshold if not specified by User
if ~exist('sz_thresh', 'var')
    sz_thresh = ceil(N_phenos/3);
end

% cluster gene scores 
if strcmp(scale,'raw')
    Z = linkage(pegasus_scores',cluster_alg, cluster_dist);
elseif strcmp(scale,'log10')
    Z = linkage(-log10(pegasus_scores)',cluster_alg, cluster_dist);
else
    fprintf('Error: scale must be "raw" or "log10"')
    return
end

% define cluster key for branch lenghts (pointers to parents, exclude singleton clusters)
Z_key = N_phenos+1:2*N_phenos -1 ;

% define empty variable for branch lengths
branch_lengths = zeros(N_phenos - 1, 1);

% go through each internal branch and find the corresponding branch length
% using the linkage output (note: root length = 0).
for k = 1:N_phenos-2
    if ~isempty(find(Z(:,1) == Z_key(k))) %pointer in first column
        branch_lengths(k) = abs(Z(find(Z(:,1) == Z_key(k)),3) - Z(k,3));
    else %pointer in second column
        branch_lengths(k) = abs(Z(find(Z(:,2) == Z_key(k)),3) - Z(k,3));
    end
end

% define cell array for all clusters 
all_clusters = {};
all_cluster_labels = {};

% sort branch lengths
[sorted_branches, sort_ind] = sort(branch_lengths, 'descend');

% go through and find complete cluster of phenotypes
% associated with each branch 

for k = 1:length(branch_lengths)
    members =  Z(sort_ind(k),1:2);
    while max(members) > N_phenos
        % replace each internal branch with children until all members are
        % individual phenoypes
        children = Z(find(Z_key == max(members)),1:2);
        members = [members(members ~= max(members)) children];
    end
    
    all_clusters{end + 1} = members;
    temp_labels = join(pheno_labels(members));
    all_cluster_labels{end + 1} =temp_labels{1};
 
end

% now only extract clusters of size < sz_thredh
small_clusters = [];
for k = 1:N_phenos - 1
    if size(split(all_cluster_labels{k}),1) <  sz_thresh
        small_clusters(end  +1) = k;
    end
end

% extract branch lengths corresponding to "true" clusters
sig_branches = sorted_branches(small_clusters);

% find the branch length threshold
branchLength_gap = get_BL_thresh(sig_branches);

% identify all branches whose branch length is above the branch length threshold and
% find the corresponding clusters
big_branches = find(branch_lengths >= branchLength_gap);
clusters = {};
if ~isempty(big_branches)
    for k = 1:length(big_branches)
        members =  Z(big_branches(k),1:2);
        while max(members) > N_phenos
            %        replace each internal branch with children until all members are
            %        individual phenoypes
            children = Z(find(Z_key == max(members)),1:2);
            members = [members(members ~= max(members)) children];
        end
        if size(members,2) > 1 && size(members,2) < sz_thresh
            clusters{end + 1} = members;
        end
    end
    
else
    fprintf('Warning: no clusters identified')

end
N = length(clusters);

%get labels for each cluster
cluster_labels = cell(N,1);
for k = 1:N
    temp= join(pheno_labels(clusters{k}),',');
    temp = split(temp, ',');
    cluster_labels{k} = temp;
end

% plot dendrogram and ranks clusters with branch lengths if specified
if plotOn == 1
     plot_outputs(Z, pheno_labels, N, clusters, small_clusters, sorted_branches, branchLength_gap, all_cluster_labels, outdir)
end

% save outputs to text files if specified
if saveOn == 1
    save_outputs(outdir, branchLength_gap, cluster_labels, all_cluster_labels, sorted_branches)
end

end

