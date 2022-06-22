% Used for significance testing (significance_testing.m), this function 
% determines the pairs of ROIs that are functionally connected in 
% resting-state HCP data.

% Inputs:
% - T_HCP: Number of timepoints/TRs in the HCP dataset.
% - TR: Average number of TRs in the psychedelic dataset.
% - ts: HCP timeseries (100x100x1189 tensor).
% - num_rois: Number of ROIs in the parcellation of the HCP dataset, which
%   is 100.
% - N_HCP: Number of subjects in the HCP dataset.
% - N: Average number of subjects in the psychedelic datatset.
% - surrogate_study: Matrix of pairs of ROIs that are functionally 
%   connected, will be updated by compute_surrogate_study. Note that 
%   significance_testing.m iterates compute_surrogate_study several times, 
%   once for each study in the psychedelic dataset. (nx3 matrix, where n is 
%   the number of functionally connected pairs, columns 1 and 2 are the 
%   indices of the ROIs, and column 3 is the sign of the connectivity.)

% Outputs:
% - surrogate_study: See above.

function surrogate_study = compute_surrogate_study(T_HCP, TR, ts, num_rois, N_HCP, N, surrogate_study)         
  
start_window1 = randi(T_HCP - TR);
end_window1 = start_window1 + TR;
start_window2 = randi(T_HCP - TR);
end_window2 = start_window2 + TR;
sampled_subs = randi(N_HCP, [1,N]);
num_subs = length(sampled_subs);
% rs is a cell array where each element is a matrix
% consisting of the (i,j)th element of all matrices in rs
rs1 = zeros(num_subs,1,100,100);
rs2 = zeros(num_subs,1,100,100);
dowhile = true;

while dowhile
    
    for j = 1:num_subs
        window1 = ts{sampled_subs(j)}(:,start_window1:end_window1);
        window2 = ts{sampled_subs(j)}(:,start_window2:end_window2);
        c1 = [atanh(corr(window1'))];
        c2 = [atanh(corr(window2'))];
        clear window1 window2
        c1(isinf(c1)) = 1;
        c2(isinf(c2)) = 1;
        for r = 1:num_rois
            for c = 1:num_rois
                rs1(j,1,r,c) = c1(r,c);
                rs2(j,1,r,c) = c2(r,c);
            end
        end
        clear c1 c2
    end

    for r = 1:num_rois
        for c = 1:num_rois
            if r >= c
                fc2 = rs2(:,:,r,c);
                fc1 = rs1(:,:,r,c);
                [h,p] = ttest(fc2, fc1);
                if p <= 0.01
                    surrogate_study = [surrogate_study; [r, c, sign(mean([fc2; fc1]))]];
                    dowhile = false;
                end
                clear fc2 fc1 
            end
        end
    end
    clear rs1 rs2
    
end