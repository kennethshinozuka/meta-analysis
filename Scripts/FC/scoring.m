% Function that updates the scoring matrix given two functionally connected
% ROIs in a study and their overlaps with the Yeo networks. The scoring
% algorithm is as follows:

% Partial:
% C(y1,y2) = C(y2,y1) = (s(ROI1||ROI2) * (1/n)(o_y1 + o_y2)/2)
% where y1, y2 are the Yeo networks that either ROI #1 or #2 overlaps with,
% s(ROI1||ROI2) is the sign of the connectivity between ROIs #1 and #2,
% n is the number of ROIs in the study, o_y1 is the relative overlap
% between ROI #1 and y1, and o_y2 is the relative overlap between ROI #2
% and y2. Note that any changes made to C(y1,y2) are also made to C(y2,y1),
% such that the scoring matrix is symmetric.

% Winner-takes-all (WTA):
% C(y1,y2) = C(y2,y1) = s(ROI1||ROI2) * 1/n
% (see above for definition of variables)

% Inputs:
% - datatype: A study's method for measuring functional connectivity. This
%   is either 'other' or 'seed-to-voxel'.
% - signs: The sign of the connectivity between the ROIs. This is either a
%   scalar (1 or -1) when the datatype is 'other' or a nx1 vector if the 
%   datatype is 'seed-to-voxel', where n is the number of Yeo networks that 
%   contain at least one voxel in the spatial map of seed-to-voxel 
%   connectivity.
% - o1p: The overlaps between ROI #1 and the Yeo networks according to the
%   partial paradigm. This is a 1x7 vector where each element is the
%   proportion of voxels in the ROI that overlap with a particular Yeo 
%   network.
% - o2p: The overlaps between ROI #2 and the Yeo networks according to the
%   partial paradigm (1x7 vector).
% - o1w: The overlaps between ROI #1 and the Yeo networks according to the
%   WTA paradigm. This is a scalar representing the 
%   single Yeo network that contains the most voxels in the ROI.
% - o2w: The overlaps between ROI #2 and the Yeo networks according to the
%   WTA paradigm (scalar).
% - C_p: The partial scoring matrix that is to be updated through the
%   algorithm. The initial state of C_p, before any studies have been
%   scored, is a 7x7 matrix of zeroes.
% - C_w: The WTA scoring matrix that is to be updated through the
%   algorithm. The initial state of C_p, before any studies have been
%   scored, is a 7x7 matrix of zeroes.
% - num_rois: The number of pairs of functionally connected ROIs in the 
%   study.

% Outputs:
% - C_p: The updated partial scoring matrix (7x7 matrix).
% - C_w: The updated WTA scoring matrix (7x7 matrix).

function [C_p, C_w] = scoring(datatype, signs, o1p, o2p, o1w, o2w, C_p, C_w, num_rois)

    % Update the partial scoring matrix
    
    % If neither ROI #1 nor ROI #2 overlap with any of the Yeo networks,
    % then return the inputted partial scoring matrix
    if (1 && all(o1p == 0)) | (1 && all(o2p == 0))
        C_p;
    else
        nonzero_o1p = o1p(find(o1p~=0)); % amount of nonzero overlap between ROI #1 and each Yeo network
        inds_nonzero_o1p = find(o1p~=0); % Yeo networks that have nonzero overlap with ROI #1
        nonzero_o2p = o2p(find(o2p~=0));
        inds_nonzero_o2p = find(o2p~=0);
        arr_nonzero_o1p = [nonzero_o1p; inds_nonzero_o1p];
        arr_nonzero_o2p = [nonzero_o2p; inds_nonzero_o2p];
        % Cycle through each pair of Yeo networks that ROI #1 or ROI #2
        % overlaps with.
        for ind1 = 1:size(arr_nonzero_o1p,2)
            for ind2 = 1:size(arr_nonzero_o2p,2)
                y1 = arr_nonzero_o1p(2,ind1); % Yeo network that ROI #1 overlaps with
                y2 = arr_nonzero_o2p(2,ind2); % Yeo network that ROI #2 overlaps with
                if strcmp(datatype, 'seed_to_voxel')
                    % In seed-to-voxel connectivity, ROI #2 is always the 
                    % map of voxels that is connected to the seed (ROI #1).
                    s2p = find(signs(1,:) == y2); % sign of the connectivity between ROI #1 and the region of ROI #2 that overlaps with the Yeo network at hand
                    score_p = signs(2,s2p) * (arr_nonzero_o1p(1,ind1) + arr_nonzero_o2p(1,ind2))/2/num_rois;
                else
                    score_p = signs * (arr_nonzero_o1p(1,ind1) + arr_nonzero_o2p(1,ind2))/2/num_rois;
                end
                C_p(y1,y2) = C_p(y1,y2) + score_p;
                C_p(y2,y1) = C_p(y2,y1) + score_p;
            end
        end
    end
    
    % Update the WTA scoring matrix
    
    % If neither ROI #1 nor ROI #2 overlap with any of the Yeo networks,
    % then return the inputted partial scoring matrix
    if o1w == 0 | o2w == 0
        C_w;
    else
        % In seed-to-voxel connectivity, ROI #2 is always the 
        % map of voxels that is connected to the seed (ROI #1).
        if strcmp(datatype, 'seed_to_voxel') 
            s2w = find(signs(1,:) == o2w); % sign of the connectivity between ROI #1 and the region of ROI #2 that overlaps with the Yeo network at hand
            score_w = signs(2,s2w) * 1/num_rois;
        else
            score_w = signs * 1/num_rois;
        end
        C_w(o1w,o2w) = C_w(o1w,o2w) + score_w;
        C_w(o2w,o1w) = C_w(o2w,o1w) + score_w;
    end

end



