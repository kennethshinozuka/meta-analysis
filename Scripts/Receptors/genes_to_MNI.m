% This code computes the expression of the gene of each receptor, as
% determined by the Allen Human Brain Atlas (AHBA), in each Yeo network.
% We are implementing the algorithm by Burt et al. (2018): 
% https://www.nature.com/articles/s41593-018-0195-0

root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/Code Review/';
addpath(fullfile(root_path, 'Scripts', 'Receptors', 'AHBA'))

% IDs of the six postmortem human brains in AHBA 
specimen_names{1} = 'H0351.2001';
specimen_names{2} = 'H0351.2002';
specimen_names{3} = 'H0351.1009';
specimen_names{4} = 'H0351.1012';
specimen_names{5} = 'H0351.1015';
specimen_names{6} = 'H0351.1016';

% genes encoding receptors that at least one psychedelic (LSD, psilocybin,
% or ayahuasca) binds to
p.receptor = {'5-HT1A', '5-HT1B', '5-HT1D', '5-HT1E', '5-HT2A', '5-HT2B', '5-HT2C', '5-HT5', '5-HT6', '5-HT7', 'D1', 'D2', 'D3', 'D4', 'D5', 'ADRA1A', 'ADRA1B', 'ADRA2A', 'ADRA2B', 'ADRA2C', 'ADRB1', 'ADRB2', 'H1', 'MAO-A'};
% IDs of probes used to measure the expression of each gene
p.probe = {};
[num_probes, txt_probes, raw_probes] = xlsread(fullfile(root_path, 'Data', 'Receptors', 'gene_probes.xlsx'));
for c = 1:size(txt_probes,2)
    a = txt_probes(2:end,c);
    p.probe{c} = a(~cellfun('isempty', a));
end

% some functions below provided by: http://api.brain-map.org/examples/brodmann/index.html
gene_explevels = [];
for r = 1:size(p.receptor,2)
    subj_explevels = [];
    for s = 1:size(specimen_names,2)
        specimen = download_specimen(specimen_names{s});
        scores_probes = []; 
        explevels_array = []; % expression levels across all probes for a single subject
        samples_array = {}; % sample data across all probes for a single subject, which includes the MRI coordinates of the sample
        for pr = 1:size(p.probe{r},1)
            probe_name = p.probe{r}{pr};
            probe = download_probe(probe_name);
            [samples, explevels] = download_expression(probe.id, specimen.donor_id); % expression levels of the gene at each sample measured by the probe
            samples_array{pr} = samples; 
            explevels_array = [explevels_array explevels];
        end
        % If there are two probes for a subject, use the one with the
        % higher variance.
        if size(p.probe{r},1) == 2
            vars = [var(explevels_array(:,1)) var(explevels_array(:,2))]; 
            best_probe_ind = find(vars==max(vars));
        % If there are more than two probes, use the one whose expression
        % levels have the highest correlation with those of other probes.
        else
            corrs_explevels = corrcoef(explevels_array);
            similarity = sum(corrs_explevels,1);
            best_probe_ind = find(similarity==max(similarity));
        end
        best_probe = p.probe{r}{best_probe_ind};
        best_explevels = explevels_array(:,best_probe_ind);
        best_samples = samples_array{best_probe_ind};
        best_coords = transform_samples(best_samples,specimen.alignment3d)'; % MNI coordinates of the samples measured by the best probe
        mean_yeo_explevels = calculate_overlap_genes(best_coords, best_explevels); % mean expression level of the gene in each Yeo network for one subject
        subj_explevels = [subj_explevels zscore(mean_yeo_explevels)]; % vector of mean expression levels of the gene in each Yeo network across all subjects
    end
    gene_explevels = [gene_explevels zscore(mean(subj_explevels, 2))]; % mean expression level of the gene in each Yeo network, averaged across subjects
end
