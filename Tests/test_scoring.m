root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/CodeReview/meta-analysis';
cd(root_path)
addpath(fullfile(root_path, 'Scripts'))

C_p = zeros(7);
C_w = zeros(7);
s = -1;
num_rois = 2;
num_secondaries = 4;
op1 = [0.1 0 0 0.4 0 0.5 0]; % partial overlap vector of ROI 1, in which each element corresponds to the proportion of a Yeo network that overlaps with the ROI
op2 = [0.3 0.5 0 0.1 0 0.1 0]; % partial overlap vector of ROI 2, in which each element corresponds to the proportion of a Yeo network that overlaps with the ROI
ow1 = 6; % winner-takes-all (WTA) overlap of ROI 1, corresponding to the Yeo network that overlaps most with ROI 1 
ow2 = 2;% winner-takes-all (WTA) overlap of ROI 2, corresponding to the Yeo network that overlaps most with ROI 2
[C_p, C_w, ~, ~] = scoring('other', s, op1, op2, ow1, ow2, C_p, C_w, C_p, C_w, num_rois, num_secondaries);
assert(all(all(C_p == [-0.02 -0.03 0 -0.0450 0 -0.05 0; -0.03 0 0 -0.045 0 -0.05 0; 0 0 0 0 0 0 0; -0.045 -0.045 0 -0.025 0 -0.055 0; 0 0 0 0 0 0 0; -0.05 -0.05 0 -0.055 0 -0.03 0; 0 0 0 0 0 0 0])));
assert(all(all(C_w == [0 0 0 0 0 0 0; 0 0 0 0 0 -0.1 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 -0.1 0 0 0 0 0; 0 0 0 0 0 0 0]))); 

% Note: these are the same overlap vectors that were created in
% test_calculate_overlap.m between ROIs situated in various Yeo networks.
op_vislim = [0.25 0 0 0 0.75 0 0];
ow_vislim = 5;
op_visFPN = [0.5 0 0 0 0 0.5 0];
ow_visFPN = 1;
[C_p, C_w, ~, ~] = scoring('other', s, op_vislim, op_visFPN, ow_vislim, ow_visFPN, C_p, C_w, C_p, C_w, num_rois, num_secondaries);
assert(all(all(round(C_p,4) == round([-0.0575 -0.03 0 -0.0450 -0.0625 -0.0875 0; -0.03 0 0 -0.045 0 -0.05 0; 0 0 0 0 0 0 0; -0.045 -0.045 0 -0.025 0 -0.055 0; -0.0625 0 0 0 0 -0.0625 0; -0.0875 -0.05 0 -0.055 -0.0625 -0.03 0; 0 0 0 0 0 0 0],4))));
% above, we round to the ten-thousandths place because of floating-point
% arithmetic in MATLAB, which sometimes causes numbers to be very slightly
% different from their true values