% Used for significance testing (significance_testing.m), this function 
% determines the pairs of ROIs that are functionally connected in 
% resting-state HCP data. For each surrogate study, we randomly select two
% windows of HCP data and compute FC (the Fisher z-transformed Pearson
% correlation) between each window for N HCP subjects, where N is the
% average number of subjects in a psychedelic FC study. We then determine
% which ROIs are significantly correlated across windows. 

% Inputs:
% - T_HCP: Number of timepoints/TRs in the HCP dataset.
% - TR: Average number of TRs in the psychedelic dataset.
% - ts: HCP timeseries (100x100x1189 tensor).
% - num_rois: Number of ROIs in the DK parcellation of the HCP dataset.
%   num_rois must be either 62 or 80; if num_rois is 62, then we measure FC
%   in only the cortical ROIs, but if it's 80 then we measure FC across all
%   ROIs, including subcortical ones. 
% - N_HCP: Number of subjects in the HCP dataset (100).
% - N: Average number of subjects in the psychedelic datatset.

% Outputs:
% - surrogate_study: Matrix of pairs of ROIs that are functionally 
%   connected, will be updated by compute_surrogate_study. Note that 
%   significance_testing.m iterates compute_surrogate_study several times, 
%   once for each study in the psychedelic dataset. (nx3 matrix, where n is 
%   the number of functionally connected pairs, columns 1 and 2 are the 
%   indices of the ROIs, and column 3 is the sign of the connectivity.)

function surrogate_study = compute_surrogate_study(T_HCP, TR, ts, num_rois, N_HCP, N)         

% Randomly select two windows of HCP data of length TR
rands = randperm(T_HCP - TR, 2);
start_window1 = rands(1);
end_window1 = start_window1 + TR;
start_window2 = rands(2);
end_window2 = start_window2 + TR;

% Randomly select N subjects out of the 100 HCP subjects
sampled_subs = randperm(N_HCP, N);
% rs is a tensor which contains the num_rois x num_rois correlation matrices
% for all subjects.
rs1 = zeros(N,num_rois,num_rois);
rs2 = zeros(N,num_rois,num_rois);
surrogate_study = [];
dowhile = true;

while dowhile
    
    for j = 1:N
        if num_rois == 62
            % select only the cortical ROIs
            window1 = [ts{sampled_subs(j)}.signal_filt(1:31,start_window1:end_window1); ts{sampled_subs(j)}.signal_filt(50:80,start_window1:end_window1)];
            window2 = [ts{sampled_subs(j)}.signal_filt(1:31,start_window2:end_window2); ts{sampled_subs(j)}.signal_filt(50:80,start_window2:end_window2)];
        elseif num_rois == 80
            % select both cortical and subcortical ROIs
            window1 = ts{sampled_subs(j)}.signal_filt(:,start_window1:end_window1);
            window2 = ts{sampled_subs(j)}.signal_filt(:,start_window2:end_window2);
        else
            error('num_rois must be either 62 or 80');
        end
        c1 = [atanh(corr(window1'))]; % atanh = Fisher z-transformation
        c2 = [atanh(corr(window2'))];
        clear window1 window2
        c1(isinf(c1)) = 1;
        c2(isinf(c2)) = 1;
        rs1(j,:,:) = c1; % store each subject's correlation matrix in the rs tensor 
        rs2(j,:,:) = c2;
        clear c1 c2
    end
    
    % Determine which ROIs are significantly connected in the resting-state
    % data. We use a Wilcoxon signed-rank test to determine significance.
    sigFC = zeros(num_rois);
    for r = 1:num_rois
        for c = 1:num_rois
            if r > c % (correlation matrix is symmetric)
                p = signrank(rs1(:,r,c),rs2(:,r,c));
                if p <= 0.05
                    sigFC(r,c) = 1;
                end
            end
        end
    end
    
    clear rs1 rs2

    if any(any(sigFC)) % if there are any nonzero elements in sigFC, continue; otherwise, repeat the while loop
        [r,c] = ind2sub(size(sigFC),find(sigFC==1)); % indices of elements of sigFC that are equal to 1 (ROIs that are signficantly connected)
        surrogate_study = [surrogate_study; [r,c]]; % append these ROIs to the surrogate_study vector
        dowhile = false;
    end
        
end