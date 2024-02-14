%% This loops through stim sites for each subject and saves CARLA results, with fewer other outputs along the way
% Then run fig7_summarizeCARLARealCCEPs to summarize CARLA results for Figure 7
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
%% Paths, load data for each CCEP run

outdir = fullfile('output', 'realCCEPsLoop');

% subjects and runs corresponding to each subject to load
subs = {'sub-1', 'sub-2', 'sub-3', 'sub-4'};
runs = {[1, 2, 3, 5], [12, 13, 14], [2, 3, 4], [1, 2]};

cmSens = [1, 165/255, 0]; % orange color used for detecting optimum

for ss = 1:length(subs)

    sub = subs{ss};
    mkdir(outdir, sub);

    % iterate across runs in subject
    for rr = 1:length(runs{ss})
        
        run = runs{ss}(rr);
    
        % known issue: Matlab crashing when readMef3 loads mef data from relative path. Temporary fix = fullfile with pwd. If you specified an absolute "dataPath", remove "pwd".
        mefPath = fullfile(pwd, dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_ieeg.mefd', sub, run));
        channelsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_channels.tsv', sub, run));
        eventsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_events.tsv', sub, run));
        annotPath = fullfile(dataPath, 'derivatives', 'event_annotations', sub, sprintf('%s_ses-ieeg01_task-ccep_run-%02d_eventsAnnot.tsv', sub, run));
        elecsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_electrodes.tsv', sub));
    
        % Load CCEP data from each subject
        mef = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
        mef.loadMefAll;
        mef.highpass;
        mef.loadMefTrials([-1, 1]);
        
        tt = mef.tt;
        srate = mef.srate;
        channels = mef.channels;
        events = mef.evts;
        data = mef.data;
    
        sites = groupby(events.electrical_stimulation_site);
    
        % read manual annotations on events
        clear annots
        try % subject 4 has no annots because no interictal activity, hence the try
            annots = readtable(annotPath, 'FileType', 'text', 'Delimiter', '\t');
            annots = annots.status_description;
        catch
            warning('Unable to load event annotation for events');
        end
        
        % check that annots match events and remove bad events from annots
        if exist('annots', 'var')
            eventsTemp = ieeg_readtableRmHyphens(eventsPath, 'electrical_stimulation_site', 1); % to remove bad events from annotations
            assert(length(annots) == height(eventsTemp), 'Error: length of event annotations does not match events table height');
            annots(~strcmpi(eventsTemp.status, 'good')) = [];
        end        
            
        elecs = ieeg_readtableRmHyphens(elecsPath);
        SOZ = elecs.name(contains(elecs.seizure_zone, 'SOZ'));
        disp('SOZs:');
        disp(SOZ');
    
        %% Iterate across sites in current run
    
        for kk = 1:size(sites, 1)
    
            site = sites{kk};
            idxesSite = sites{kk, 2};
    
            if any(ismember(split(site, '-'), SOZ))
                fprintf('Skipping stim site %s in SOZ\n', site);
                continue
            end
            fprintf('Current site: %s\n', site);
    
            % extract good channels at this stim site
            goodChs = strcmp(channels.type, 'SEEG') & strcmp(channels.status, 'good') & ~ismember(upper(channels.name), getNeighborChs(split(site,'-'), 2)); % exclude sites 2 away from current
            dataSite = data(goodChs, :, idxesSite);
            chsSite = channels.name(goodChs);
    
            try
                % remove bad trials and set channels to nan in singular trials, according to annotations
                annotsSite = annots(idxesSite); % get annotations corresponding to this stim site
                dataSite = rmBadTrialsAnnots(dataSite, chsSite, annotsSite);
            end
            
            % minimum floor on number of trials for any stim site
            if size(dataSite, 3) < 8
                fprintf('Skipping stim site %s with %d trials\n', site, size(dataSite, 3));
                continue
            end
    
            %% Remove channels that have any nan trials, apply CARLA
    
            rng('default');
    
            nanChs = any(isnan(dataSite), [2, 3]); % any channels with nan in any trial
            dataSiteNoNan = dataSite(~nanChs, :, :);
            chsSiteNoNan = chsSite(~nanChs);
            
            % calculate average interchannel R-squared at each possible adjusted CAR size
            rsqs = crossChCorr(tt, dataSiteNoNan, srate);
            save(fullfile(outdir, sub, sprintf('%s_%s_interChRsq.mat', sub, site)), 'rsqs');
    
            % Perform CAR
            [Vcar, CAR, stats] = CARLA(tt, dataSiteNoNan, srate, true);
            nCAR = length(stats.chsUsed);
    
            % ZminMean plot
            f = figure('Position', [200, 200, 600, 300]); hold on
            errorbar(mean(stats.zMinMean, 2), std(stats.zMinMean, 0, 2), 'k.-', 'MarkerSize', 10, 'CapSize', 1); % SD of bootstrapped means
            plot(nCAR, mean(stats.zMinMean(nCAR, :), 2), '*', 'Color', cmSens);
            xlim([0, length(chsSiteNoNan)+1]); ylim([-inf, 0]);
            yline(0, 'Color', 'k');
            saveas(gcf, fullfile(outdir, sub, sprintf('%s_%s_zMinMean', sub, site)), 'svg');
            saveas(gcf, fullfile(outdir, sub, sprintf('%s_%s_zMinMean', sub, site)), 'png');
            close(f);
            
            % save the optimum CAR size and number of total channels to text
            fid = fopen(fullfile(outdir, sub, sprintf('%s_%s_nCAR.txt', sub, site)), 'w');
            fprintf(fid, '%d\n%d', nCAR, length(chsSiteNoNan));
            fclose(fid);
            
    
        end
        
    end
    
end