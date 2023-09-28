%% This loads the outputs from applyCARLARealCCEPsLoop to summarize the optimal CCEP cutoff size for stim sites
%   in each tissue type across the 4 subjects
%
%   2023/09/27
%
%    If this code is used in a publication, please cite the manuscript:
%    "CARLA: Adjusted common average referencing for cortico-cortical evoked potential data"
%    by H Huang, G Ojeda Valencia, NM Gregg, GM Osman, MN Montoya,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    CARLA manuscript package.
%    Copyright (C) 2023 Harvey Huang
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%

% subjects and runs corresponding to each subject to load
subs = {'sub-1', 'sub-2', 'sub-3', 'sub-4'};
runs = {[1, 2, 3, 5], [12, 13, 14], [2, 3, 4], [1, 2]};

for ss = 1:length(subs)

    sub = subs{ss};
    fprintf('Sub: %s\n', sub);

    %% Obtain stim sites across all runs, number of trials
    
    elecsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_electrodes.tsv', sub));
    elecs = ieeg_readtableRmHyphens(elecsPath);
    SOZ = elecs.name(contains(elecs.seizure_zone, 'SOZ'));
    
    sites = cell(0, 2); % name, number of trials
    for rr = 1:length(runs{ss})
        run = runs{ss}(rr);
        fprintf('%s run %d\n', sub, run);
        
        eventsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_events.tsv', sub, run));
        annotPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_eventsAnnot.tsv', sub, run));
    
        events = ieeg_readtableRmHyphens(eventsPath, 'electrical_stimulation_site', 1);
        eventsStatusOrig = events.status;
        events(~strcmpi(events.status, 'good'), :) = [];
        
        sitesRun = groupby(events.electrical_stimulation_site); % all sites in current run
        
        % read manual annotations on events if it exists
        try
            annots = readtable(annotPath, 'FileType', 'text', 'Delimiter', '\t');
            annots = strip(annots.status_description);
        catch
            warning('Unable to load event annotation for events');
            clear annots
        end
        
        if exist('annots', 'var')
            annots(~strcmpi(eventsStatusOrig, 'good')) = [];
            assert(length(annots) == height(events), 'Error: length of event annotations does not match events table height');
        end   
        
        % Remove stim sites in SOZ
        isSoz = any(ismember(split(sitesRun(:, 1), '-', 2), SOZ), 2);
        sitesRun(isSoz, :) = [];
        
        % remove "all" annotated trials and calculate final number of trials
        nTrials = cellfun(@length, sitesRun(:, 2));
        if exist('annots', 'var')
            for ii = 1:size(sitesRun, 1)
                annotsCurr = annots(sitesRun{ii, 2});
                nTrsRm = sum(strcmpi(annotsCurr, 'all'));
                nTrials(ii) = nTrials(ii) - nTrsRm;
            end
        end
        
        sites = [sites; sitesRun(:, 1), num2cell(nTrials)];
         
    end
    
    % ensure no stim site has less than 8 trials
    sites = cell2table(sites, 'VariableNames', {'site', 'trials'});
    sites(sites.trials < 8, :) = [];
    
    %% Identify Destrieux labels for each stim site
    
    % subject-specific Freesurfer derivatives folder
    FSdir = fullfile(dataPath, 'derivatives', 'freesurfer', sub);
    xyzs = ieeg_getPairXyzs(sites.site, elecs);
    labs = ieeg_getLabelXyzDestrieux(xyzs, FSdir); % within radius of 3
    
    labGroup = repmat({'n/a'}, length(labs), 1);
    labGroup(contains(labs, {'Thalamus', 'Hippocampus', 'Amygdala'})) = {'subgray'};
    labGroup(contains(labs, 'White_Matter')) = {'WM'};
    labGroup(startsWith(labs, {'lh', 'rh'})) = {'cortgray'};
    
    fprintf('Subcortical Gray: %d\nCortical Gray: %d\nWhite Matter: %d\nn/a: %d\n', ...
            sum(strcmp(labGroup, 'subgray')), sum(strcmp(labGroup, 'cortgray')), sum(strcmp(labGroup, 'WM')), sum(strcmp(labGroup, 'n/a')));
    
    sites.labs = labs;
    sites.labGroup = labGroup;
    
    %% load the nCAR files and plot distributions
    
    rng('default');
    
    outdir = fullfile('output', ''realCCEPsLoop', sub);
    
    carSz = zeros(height(sites), 2);
    for ii = 1:height(sites)
        M = readmatrix(fullfile(outdir, sprintf('%s_%s_nCAR.txt', sub, sites.site{ii})));
        carSz(ii, :) = M;
    end
    
    sites.carSz = carSz(:, 1);
    sites.numChs = carSz(:, 2);
    sites.pct = 100*carSz(:, 1)./carSz(:, 2);
    
    % plot overall distribution
    figure; histogram(sites.pct, 10);
    
    % plot distributions by label group and hemisphere
    pctByGroup = cell(6, 1);
    
    % Right hemi
    pctByGroup{1} = sites.pct(strcmp(labGroup, 'cortgray') & startsWith(lower(labs), 'r'));
    pctByGroup{2} = sites.pct(strcmp(labGroup, 'subgray') & startsWith(lower(labs), 'r'));
    pctByGroup{3} = sites.pct(strcmp(labGroup, 'WM') & startsWith(lower(labs), 'r'));
    
    % Left hemi
    pctByGroup{4} = sites.pct(strcmp(labGroup, 'cortgray') & startsWith(lower(labs), 'l'));
    pctByGroup{5} = sites.pct(strcmp(labGroup, 'subgray') & startsWith(lower(labs), 'l'));
    pctByGroup{6} = sites.pct(strcmp(labGroup, 'WM') & startsWith(lower(labs), 'l'));
    
    
    figure('Position', [200, 200, 250, 300]); hold on
    bar(1:3, cellfun(@mean, {[pctByGroup{1}; pctByGroup{4}], [pctByGroup{2}; pctByGroup{5}], [pctByGroup{3}; pctByGroup{6}]}), ...
        'FaceColor', [0.8, 0.8, 0.8]); % mean for stim sites in both hemispheres
    hsR = jitterplot(pctByGroup(1:3), brighten([252, 141, 98]/255, -0.5));
    set(hsR, 'Marker', '.', 'MarkerSize', 16);
    hsL = jitterplot(pctByGroup(4:6), brighten([102, 194, 165]/255, -0.5));
    set(hsL, 'Marker', '.', 'MarkerSize', 16);
    ylim([0, 100]);
    set(gca, 'XTick', 1:3, 'XTickLabels', {'cortgray', 'subgray', 'WM'});
    
    [p, tbl, stats] = anova1(sites.pct, sites.labGroup, 'off');
    fprintf('ANOVA p = %0.02e\n', p);
    
    saveas(gcf, fullfile(outdir, 'summaryPctCAR'), 'png');
    saveas(gcf, fullfile(outdir, 'summaryPctCAR'), 'svg');
    
    writetable(sites, fullfile(outdir, 'summaryTbl.tsv'), 'Delimiter', '\t', 'FileType', 'text');

end
