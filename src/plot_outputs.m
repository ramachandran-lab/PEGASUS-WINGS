%% function: plot_outputs(Z, pheno_labels, N, clusters, small_clusters, sorted_branches, branchLength_gap, all_cluster_labels, outdir)
%%
%% Description: This is a function that takes as input the outputs from cluster_phenos and plots the corresponding dendrogram and sorted 
%% clusters by branch length. 
%% 
%% Corresponding Author: Melissa R. McGuirl, melissa_mcguirl@brown.edu

function [] =  plot_outputs(Z, pheno_labels, N, clusters, small_clusters, sorted_branches, branchLength_gap, all_cluster_labels, outdir)

f = figure(1);
d = dendrogram(Z,0,'Labels', pheno_labels);
set(d,'LineWidth',2.5)

cmap = hsv(N);
set( d, 'Color', 'k');
for k = 1:N
    cluster_color = cmap(k,:);
    for j = 1:length(clusters{k})
        if  ~isempty(find(Z(:,1) == clusters{k}(j)))
            set( d(find(Z(:,1) == clusters{k}(j))), 'Color', cluster_color );
        else
            set( d(find(Z(:,2) == clusters{k}(j))), 'Color', cluster_color );
        end
    end
end
xtickangle(90)
title('Dendrogram')
set(gca,'fontSize', 10)
saveas(f,[outdir '/WINGS_Output_dendro.jpg'])

g = figure(2);
hold on
scatter([1:1:length(small_clusters)], sorted_branches(small_clusters), 150,'filled');
hold on
plot([1:1:length(small_clusters)],branchLength_gap*ones(1,length(small_clusters)) , 'r--', 'LineWidth', 4)
ylabel('Branch Length')
xlabel('Identified Clusters')
xticks([1:1:length(small_clusters)])
xticklabels( all_cluster_labels(small_clusters));
xtickangle(90)
legend('Branch Lengths', 'Branch Length Cut Off')
title('Sorted Branch Lengths')
suptitle('Hierarchical Clustering Output')
set(gca,'fontSize', 10)
saveas(g,[outdir '/WINGS_Output_BL.jpg'])

end

