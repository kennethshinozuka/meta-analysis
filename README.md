# Psychedelic fMRI meta-analysis: Code review

Scripts Folder:

This folder contains a portion of the scripts that I wrote, specifically for two segments of the meta-analysis: functional connectivity and the expression of receptors that psychedelics bind to (e.g. 5-HT2A). 

FC (Functional Connectivity):

1. FC_meta_analysis.m computes the scoring matrix for the psychedelic FC dataset, which represents the overall connectivity between Yeo networks. It calls functions implemented in calculate_overlap.m, which calculates the overlaps between ROIs and Yeo networks, and scoring.m, which contains the scoring algorithm.
- Note: You need to edit the root filepath in line 45 of FC_meta_analysis.m and line 30 of calculate_overlap.m. 

2. significance_testing.m performs significance testing on the resulting scoring matrix. It calls the function compute_surrogate_study.m, which calculates FC in a null, resting-state dataset. significance_testing.m then compares the scoring matrix for the psychedelic dataset to a histogram of null FC values.
- Note 1: As written in the comments, there is a for loop in significance_testing.m (specifically, in the 'Generate surrogate study set' section) which you may have to run multiple times in order to complete all 500 iterations.
- Note 2: You need to edit the root filepath in line 20 of significance_testing.m.

3. parcellation_consistency.m checks whether the scoring matrices corresponding to different parcellations of the same resting-state data are consistent with one another. We do this because the psychedelic FC dataset utilises several different parcellations, ranging from 27 to 268 ROIs. The code calls the Bash script extract_ROI_maps.sh.
- Note 1: You will need to run extract_ROI_maps.sh on Terminal as it uses FSLmaths. There is a way to run Bash scripts that execute FSL commands on MATLAB, but for some reason I kept getting the error "fslmaths: command not found".
- Note 2: You need to edit the root filepath in line 11.

Receptors (i.e. expression profiles of receptors that psychedelics bind to):

1. genes_to_MNI.m computes the mean expression level of the genes encoding each receptor, as quantified by the Allen Human Brain Atlas (AHBA), in each Yeo network. It calls the function calculate_overlap_genes.m and several functions provided by the Allen Brain Atlas API, which are in the AHBA subdirectory.


Data Folder:

FC:

1. FC_indices.xlsx contains the names/indices of the pairs of ROIs that are functionally connected in each study. It is used in FC_meta_analysis.m. There are also FC indices for each of the individual psychedelic drugs.

2. ...coordinates.xlsx contains the MNI coordinates of ROIs in studies that define ROIs as 10mm-radius spheres around coordinates, or as networks of such spheres.

3. ...seed.nii contains spatial maps of ROIs in studies that define ROIs based on an anatomically-derived mask (e.g. from an atlas) or a data-driven network (e.g. via ICA). 

4. Diagram FC inputs.png contains a diagram that I made of the various kinds of inputs that go into the FC meta-analysis.

5. Parcellation Consistency > ...FC_matrix.csv contains the ROI-to-ROI FC values for each parcellation of resting-state data from the NYU Test-Retest dataset. parcellation_consistency.m uses this data. Parcellation Consistency > ...atlas.nii contains each parcellation tested by parcellation_consistency.m.

6. Null > schaefer100_coordinates.xlsx contains the coordinates of the 100 ROIs in the Schaefer100 parcellation, which are used to score the resting-state HCP data in significance_testing.m. Schaefer100_BOLD_HCP.mat contains the timeseries of each ROI in the resting-state HCP data.

Receptors:

(1) gene_probes.xlsx contains the IDs of the AHBA probes that were used to measure the expression of each gene that encodes a receptor that psychedelics bind to.


