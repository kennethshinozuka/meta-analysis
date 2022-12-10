root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/';
cd(root_path)
addpath(fullfile(root_path, 'Scripts'))

  % Check that all spatial maps of the individual Yeo networks overlap 100% with
% each of the Yeo networks, according to calculate_overlap

yeo_atlas_path = 'Atlases/yeo/MNI/Individual Networks';
cd(yeo_atlas_path)

[~, op_vis, ow_vis] = calculate_overlap('other', 'yeo_visual.nii.gz', 'yeo');
assert(all(op_vis == [1 0 0 0 0 0 0]));
assert(ow_vis == 1);

[~, op_SMN, ow_SMN] = calculate_overlap('other', 'yeo_SMN.nii.gz', 'yeo');
assert(all(op_SMN == [0 1 0 0 0 0 0]));
assert(ow_SMN == 2);

[~, op_VAN, ow_VAN] = calculate_overlap('other', 'yeo_VAN.nii.gz', 'yeo');
assert(all(op_VAN == [0 0 1 0 0 0 0]));
assert(ow_VAN == 3);

[~, op_DAN, ow_DAN] = calculate_overlap('other', 'yeo_DAN.nii.gz', 'yeo');
assert(all(op_DAN == [0 0 0 1 0 0 0]));
assert(ow_DAN == 4);

[~, op_lim, ow_lim] = calculate_overlap('other', 'yeo_limbic.nii.gz', 'yeo');
assert(all(op_lim == [0 0 0 0 1 0 0]));
assert(ow_lim == 5);

[~, op_FPN, ow_FPN] = calculate_overlap('other', 'yeo_FPN.nii.gz', 'yeo');
assert(all(op_FPN == [0 0 0 0 0 1 0]));
assert(ow_FPN == 6);

[~, op_DMN, ow_DMN] = calculate_overlap('other', 'yeo_DMN.nii.gz', 'yeo');
assert(all(op_DMN == [0 0 0 0 0 0 1]));
assert(ow_DMN == 7);

% Create spatial maps from coordinates of individual voxels that have known
% alignments with different Yeo networks. Note: You will need FSL installed 
% on your device.

mkdir(fullfile(root_path, 'TestSpheres'))
outDir = fullfile(root_path, 'TestSpheres');
addpath(fullfile(root_path, 'Toolboxes', 'marsbar'))
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
cd(outDir)

c_vis = [19 -62 -7];
sphereROI_vis = maroi_sphere(struct('centre', [c_vis(1) c_vis(2) c_vis(3)], 'radius', 2)); % 2 is the radius of the sphere in mm^3 (note: Yeo atlas resolution is 2mm^3)
outName_vis = fullfile(outDir, sprintf('roi_vis'));
save_as_image(sphereROI_vis, [outName_vis '.nii']);

c_lim = [0 39 -21];
sphereROI_lim = maroi_sphere(struct('centre', [c_lim(1) c_lim(2) c_lim(3)], 'radius', 4)); 
outName_lim = fullfile(outDir, sprintf('roi_lim'));
save_as_image(sphereROI_lim, [outName_lim '.nii']);

c_FPN = [47 40 17];
sphereROI_FPN = maroi_sphere(struct('centre', [c_FPN(1) c_FPN(2) c_FPN(3)], 'radius', 2));
outName_FPN = fullfile(outDir, sprintf('roi_FPN'));
save_as_image(sphereROI_FPN, [outName_FPN '.nii']);

cmd = ['fslmaths ' outName_vis ' -add ' outName_FPN ' roi_vis_FPN'];
system(cmd);
cmd = ['fslmaths ' outName_vis ' -add ' outName_lim ' roi_vis_lim'];
system(cmd);

[~, op_vislim, ow_vislim] = calculate_overlap('other', 'roi_vis_lim.nii.gz', 'yeo');
assert(all(op_vislim == [0.25 0 0 0 0.75 0 0]));
assert(ow_vislim == 5);

[~, op_visFPN, ow_visFPN] = calculate_overlap('other', 'roi_vis_FPN.nii.gz', 'yeo');
assert(all(op_visFPN == [0.5 0 0 0 0 0.5 0]));
% note ow_visFPN = 1;
