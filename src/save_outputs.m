%% function: save_outputs(outdir, branchLength_gap, cluster_labels, all_cluster_labels, sorted_branches)
%%
%% Description: This is a function that takes as input the outputs from cluster_phenos and saves text files with (1) the estimated
%% significant clusters with corresponding branch length threshold and (2) all clusters ranked by branch length.
%% (optional, default = 1).
%% 
%% Corresponding Author: Melissa R. McGuirl, melissa_mcguirl@brown.edu

function [] = save_outputs(outdir, branchLength_gap, cluster_labels, all_cluster_labels, sorted_branches)
 
   fileID = fopen([outdir '/estimated_clusters.txt'], 'w');
    fprintf(fileID,'%12s',num2str(branchLength_gap));
    fprintf(fileID,'\n');
    for k = 1:size(cluster_labels,1)
        temp = join(cluster_labels{k});
        fprintf(fileID,'%12s', temp{:});
        fprintf(fileID,'\n');
    end
    fclose(fileID);
    
    fileID = fopen([outdir '/all_clusters.txt'], 'w');
    for k = 1:size(all_cluster_labels,2)
        temp = join(all_cluster_labels{k});
        fprintf(fileID,'%12s',num2str(sorted_branches(k)));
        fprintf(fileID, '\t');
        fprintf(fileID,'%12s', temp);
        fprintf(fileID,'\n');
    end
    fclose(fileID);
    
end

