%% function: branchLength_gap = get_BL_thresh(sig_branches)
%%
%% Description: This is a function that takes as input a collection of branch lengths and identifies a distinguishing branch 
%% length threshold. 
%%
%% Input: collection of branch lengths
%%
%% Output: branch length threshold
%% 
%% Corresponding Author: Melissa R. McGuirl, melissa_mcguirl@brown.edu

function branchLength_gap = get_BL_thresh(sig_branches)


% find the branch length threshold
diffs = abs(diff(sig_branches));
TF = isoutlier(diffs, 'median');
big_branches_temp = (sig_branches > median(sig_branches));
TF = logical(TF.*big_branches_temp(1:end-1));
branchLength_gap = min(sig_branches(TF));

% check for extreme cases when branch length threshold does not exist and 
% redefine it based only on branch length gap outliers or one standard
% deviation away from the mean.
if isempty(branchLength_gap)
    fprintf('Warning: branch length threshold is insignificant\n')
    branchLength_gap = 10^15;
%     TF = isoutlier(diffs, 'median');
%     branchLength_gap = min(sig_branches(TF));
end

% if isempty(branchLength_gap)
%     fprintf('Warning: No branch length outliers \n')
%     branchLength_gap = mean(sig_branches) + std(sig_branches);
% end


end

