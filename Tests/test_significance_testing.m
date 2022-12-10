% Goal: score the FC in many windows of resting-state data and show that
% the data is significant relative to the null distribution about 10% of the
% time (either above the 95th percentile or below the 5th percentile)

root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/';
cd(root_path)
addpath(fullfile(root_path, 'Scripts'))

load(fullfile(root_path, 'Null', 'DK80_BOLD_HCP_100subs.mat'))
load(fullfile(root_path, 'Null', 'all_Cp_null_DK62_061222_nosign.mat'))
load(fullfile(root_path, 'Null', 'all_Cw_null_DK62_061222_nosign.mat'))

num_iter = 500;
num_rois = 62;
[all_Cp, all_Cw] = null_dist(num_iter,num_rois);
sig_ps = zeros(num_networks,num_networks);
sig_ws = zeros(num_networks,num_networks);
for r = 1:num_networks
    for c = 1:num_networks
        % Produce a histogram of values for each cell of the scoring
        % matrices in the surrogate study set (Cp_vals and Cw_vals are 
        % 1x500 vectors). This histogram is the null distribution for that 
        % cell.
        sig_p = 0; % should be equal to no more than 10% of 500 (50) in the end
        sig_w = 0; % should be equal to no more than 10% of 500 (50) in the end
        null_Cp_vals = squeeze(all_Cp_null(r,c,:));
        null_Cw_vals = squeeze(all_Cw_null(r,c,:));
        top95_p = prctile(null_Cp_vals,95);
        top95_w = prctile(null_Cw_vals,95);        
        bottom5_p = prctile(null_Cp_vals,5);
        bottom5_w = prctile(null_Cw_vals,5);
        for i = 1:n_iter
            if all_Cp(r,c,i) > top95_p | all_Cp(r,c,i) < bottom5_p
                sig_p = sig_p + 1;
            end
            if all_Cw(r,c,i) > top95_w | all_Cw(r,c,i) < bottom5_w
                sig_w = sig_w + 1;
            end 
        end
        sig_ps(r,c) = sig_p;
        sig_ws(r,c) = sig_w;
    end
end   

assert(mean(mean(sig_ps)) > 40);
assert(mean(mean(sig_ps)) < 50);
assert(mean(mean(sig_ws)) > 40);
assert(mean(mean(sig_ws)) < 50);
