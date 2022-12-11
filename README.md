In this repository, we share the code for the functional connectivity (FC) segment of "The Hierarchy of the Psychedelic Experience: A Systematic Review and Meta-Analysis." The code for the other segments of the meta-analysis (phenomenology, GingerALE, and pharmacology) is quite straightforward so we did not include them here, but we are happy to provide them upon request. (Please email kenneth.shinozuka@psych.ox.ac.uk if you are interested.)

In order to run any of the code in the `Scripts` and `Tests` folders, there are **two prerequisites**:
1. You MUST update the `root_path` at the top of each file, which contains the path to the directory where you have saved this repository.
2. To run `compute_surrogate_study.m` and `null_dist.m`, you MUST download the HCP timeseries from 100 unrelated subjects, run `ts = subject(1:100)` then `save(ts, fullfile(root_path, 'Null', 'DK80_BOLD_HCP_100subs.mat')`. You can get the timeseries via this link: http://www.kringelbach.org/ndte/dk80/hcp1003_REST1_LR_dbs80.mat

The `Scripts` folder contains six files/directories:
1. `calculate_overlap.m`, a function that calculates the overlaps between Yeo networks and regions of interest. 
2. `scoring.m`, a function that determines the score of the overall FC between Yeo networks (or Yeo + Tian networks, if subcortical connectivity is included). 
3. `compute_surrogate_study.m`, a function that is called by `null_dist.m` to help create the null distribution.
4. `null_dist.m`, a function that generates the null distribution.
5. `FC_meta_analysis.m`, which calls all four of the above functions. Note that **you cannot run this code** without the FC data from the literature, which we cannot share because it is not available publicly.
6. `natsortfiles`, a folder of functions for sorting files numerically. This is called by `null_dist.m`.

The `Tests` folder contains three files:
1. `test_calculate_overlap.m`, which tests `calculate_overlap.m`.
2. `test_scoring.m`, which tests `scoring.m`.
3. `test_significance_testing.m`, which tests both `compute_surrogate_study.m` and `null_dist.m`.
To run all of the tests at once, simply `cd root_path` then execute `runtests('tests/').

There are three other folders: 
1. `Atlases`, which contains the Desikan-Killiany-62 (DK-62) and Desikan-Killiany-80 (DK-80) parcellations for the null distribution, as well as the Yeo networks and Yeo + Tian networks (subcortical)
2. `Null`, which contains the output of running `null_dist.m` with both the DK-62 and DK-80 parcellations
3. `Toolboxes,` which contains the MarsBaR toolbox for producing spheres around regions of interest. This is used by `calculate_overlap.m` and the corresponding test code.
