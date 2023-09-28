%% This script generates outputs for CARLA applied to single stim sites in real data (Fig 5B-F, Fig 6C-F)
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
%% Load all data for selected subject

% 'sub' and 'run' variables should have been defined in main script

elecsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_electrodes.tsv', sub));
elecs = ieeg_readtableRmHyphens(elecsPath);

mefPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_ieeg.mefd', sub, run));
channelsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_channels.tsv', sub, run));
eventsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_events.tsv', sub, run));
annotPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_eventsAnnot.tsv', sub, run));

% Load CCEP data
mef = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
mef.loadMefAll;
mef.highpass;
mef.loadMefTrials([-1, 1]);

mef.plotOutputs(1, [], 200);

tt = mef.tt;
srate = mef.srate;
channels = mef.channels;
events = mef.evts;
data = mef.data;

sites = groupby(events.electrical_stimulation_site);

outdir = fullfile('output', 'realCCEPs');
mkdir(outdir);

try
    annots = readtable(annotPath, 'FileType', 'text', 'Delimiter', '\t');
    annots = annots.status_description;
    assert(length(annots) == height(events), 'Error: length of event annotations does not match events table height');
catch
    warning('Unable to load event annotation for events');
end

% show the seizure onset zones which are excluded from stimulation
SOZ = elecs.name(contains(elecs.seizure_zone, 'SOZ'));
disp('SOZs:');
disp(SOZ');

clear *Path*

%% Extract data from site of interest, NaN out individual bad trials

% 'site' variable should have been defined in main script

idxesSite = sites{strcmp(site, sites(:, 1)), 2}; % indices corresponding to this stim site
nTrials = length(idxesSite);

goodChs = strcmp(channels.type, 'SEEG') & strcmp(channels.status, 'good') & ~ismember(upper(channels.name), getNeighborChs(split(site,'-'), 2)); % exclude sites 2 away from current

dataSite = data(goodChs, :, idxesSite);
chsSite = channels.name(goodChs);

% remove bad trials and set channels to nan in individual trials, according to annotations
annotsSite = annots(idxesSite); % get annotations corresponding to this stim site
dataSite = rmBadTrialsAnnots(dataSite, chsSite, annotsSite);

figure('Position', [200, 200, 1200, 1200]);
subplot(1, 2, 1); % first half of channels
ieeg_plotTrials(tt, mean(dataSite(1:floor(length(chsSite)/2), :, :), 3, 'omitnan')', 200, chsSite(1:floor(length(chsSite)/2)));
xlim([-0.5, 1]);
subplot(1, 2, 2); % first half of channels
ieeg_plotTrials(tt, mean(dataSite(ceil(length(chsSite)/2):end, :, :), 3, 'omitnan')', 200, chsSite(ceil(length(chsSite)/2):end));
xlim([-0.5, 1]);

%% Apply CARLA and plot variance, zminmean plots

rng('default');

nanChs = any(isnan(dataSite), [2, 3]); % any channels with nan in any trial

dataSiteNoNan = dataSite(~nanChs, :, :);
chsSiteNoNan = chsSite(~nanChs);

[Vcar, CAR, stats] = CARLA(tt, dataSiteNoNan, srate, true);
nCAR = length(stats.chsUsed);

cmSens = [1, 165/255, 0]; % orange color for max

% ZminMean plot
figure('Position', [200, 200, 600, 300]); hold on
errorbar(mean(stats.zMinMean, 2), std(stats.zMinMean, 0, 2), 'k.-', 'MarkerSize', 10, 'CapSize', 1); % SD of bootstrapped means
plot(nCAR, mean(stats.zMinMean(nCAR, :), 2), '*', 'Color', cmSens);
xlim([0, length(chsSiteNoNan)+1]); ylim([-1.5, 0]);
yline(0, 'Color', 'k');
saveas(gcf, fullfile(outdir, sprintf('%s_%s_zMinMean', sub, site)), 'svg');
saveas(gcf, fullfile(outdir, sprintf('%s_%s_zMinMean', sub, site)), 'png');

% Variance plot
vars = stats.vars(stats.order); % sort in order
figure('Position', [200, 200, 600, 300]); hold on
plot(vars, 'k.-', 'MarkerSize', 10);
xline(nCAR + 0.5, 'Color', cmSens);
xlim([0, length(chsSiteNoNan)+1]);
saveas(gcf, fullfile(outdir, sprintf('%s_%s_covar', sub, site)), 'png');
saveas(gcf, fullfile(outdir, sprintf('%s_%s_covar', sub, site)), 'svg');

%% Plot channels sorted by CARLA order and CAR

% Plot (pre-CAR) data with line noise removed, in order of variance
dataSite2Plot = zeros(size(dataSiteNoNan));
for tr = 1:size(dataSiteNoNan, 3)
    dataSite2Plot(:, :, tr) = ieeg_notch(dataSiteNoNan(:, :, tr)', srate, 60)';
end
dataSite2Plot = dataSite2Plot(stats.order, :, :);

yspace = 120;

% make excluded channels gray
cm = zeros(length(chsSiteNoNan), 3);
cm(nCAR+1:end, :) = repmat([0.5, 0.5, 0.5], length(chsSiteNoNan)-nCAR, 1);

figure('Position', [200, 200, 1200, 1200]);
nHalf = floor(length(chsSiteNoNan)/2); % number at which to break channels into 2 columns
subplot(1, 2, 1); % first half of channels
ieeg_plotTrials(tt, mean(dataSite2Plot(1:nHalf, :, :), 3)', yspace, 1:nHalf, cm(1:nHalf, :));
xlim([-0.1, 0.5]);
subplot(1, 2, 2); % second half of channels
ieeg_plotTrials(tt, mean(dataSite2Plot(nHalf+1:end, :, :), 3)', yspace, nHalf+1:length(chsSiteNoNan), cm(nHalf+1:end, :));
xlim([-0.1, 0.5]);
saveas(gcf, fullfile(outdir, sprintf('%s_%s_sortedChs', sub, site)), 'png');
print(gcf, fullfile(outdir, sprintf('%s_%s_sortedChs', sub, site)),'-depsc2', '-r300', '-painters');

% Plot the adjusted common average itself
figure('Position', [200, 200, 600, 200]);
hold on
plot(tt, CAR, 'Color', [0.5, 0.5, 0.5]);
plot(tt, mean(CAR, 2), 'k-', 'LineWidth', 1.5);
xlim([-0.1, 0.5]); ylim([-1200, 1200]);
saveas(gcf, fullfile(outdir, sprintf('%s_%s_CAR', sub, site)), 'png');
print(gcf, fullfile(outdir, sprintf('%s_%s_CAR', sub, site)), '-depsc2', '-r300', '-painters');

%% Plot all channels in one column (for Figure S1)

% downsample plotting signals to make less taxing on illustrator
dataSite2PlotDs = nan(size(dataSite2Plot, 1), length(tt)/4, size(dataSite2Plot, 3));
for ii = 1:size(dataSite2Plot, 3)
    dataSite2PlotDs(:, :, ii) = downsample(dataSite2Plot(:, :, ii)', 4)';
end
ttDs = tt(1:4:end); % downsample time vector

yspace = 80;

figure('Position', [200, 200, 400, 1200]);
ys = ieeg_plotTrials(ttDs, mean(dataSite2PlotDs, 3)', yspace, 1:length(chsSiteNoNan), cm);
xlim([-0.1, 0.5]); ylim([ys(end) - 7*yspace, ys(1) + 3*yspace]);
saveas(gcf, fullfile(outdir, sprintf('%s_%s_sortedChs_column', sub, site)), 'png');
print(gcf, fullfile(outdir, sprintf('%s_%s_sortedChs_column', sub, site)), '-depsc2', '-r300', '-painters');

%% Plot examples of channels before and after CARLA

% chs2Plot should be defined in main script

for ii = 1:length(chs2Plot)

    ch = chs2Plot(ii);

    fprintf('Ch %s\n', chsSiteNoNan{stats.order(ch)});
    
    sigRaw = squeeze(dataSiteNoNan(stats.order(ch), :, :)); % only high-pass, no notch
    sigNotch = squeeze(dataSite2Plot(ch, :, :)); % only notch, no re-reference (order already sorted)
    sigCar = squeeze(Vcar(stats.order(ch), :, :));
    
    % plot raw signal
    figure('Position', [200, 200, 450, 450]);
    subplot(3, 1, 1); hold on;
    plot(tt, sigRaw, 'Color', [0.5, 0.5, 0.5]);
    plot(tt, mean(sigRaw, 2), 'k-', 'LineWidth', 1.5);
    xlim([-0.1, 0.5]); ylim([-1200, 1200]);
    
    % plot notched signal
    subplot(3, 1, 2); hold on;
    plot(tt, sigNotch, 'Color', [0.5, 0.5, 0.5]);
    plot(tt, mean(sigNotch, 2), 'k-', 'LineWidth', 1.5);
    xlim([-0.1, 0.5]); ylim([-600, 600]);
    
    % plot CARLA signal
    subplot(3, 1, 3); hold on;
    plot(tt, sigCar, 'Color', [0.5, 0.5, 0.5]);
    plot(tt, mean(sigCar, 2), 'k-', 'LineWidth', 1.5);
    xlim([-0.1, 0.5]); ylim([-600, 600]);
    
    saveas(gcf, fullfile(outdir, sprintf('%s_%s_exCh%d', sub, site, ch)), 'png');
    print(gcf, fullfile(outdir, sprintf('%s_%s_exCh%d', sub, site, ch)), '-depsc2', '-r300', '-painters');

end

