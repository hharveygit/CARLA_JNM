%% This script uses output from applyCARLARealCCEPsLoop to calculate cross-channel rsqs for each adjusted CAR method
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

outdir = fullfile('output', 'realCCEPsLoop');

% subjects
subs = {'sub-1', 'sub-2', 'sub-3', 'sub-4'};

%% Load .mat files and CARLA results

% sub, stim site, z-transformed interch correlation matrix, CARLA result
sites = cell(0, 4);

for ss = 1:length(subs)
    
    sub = subs{ss};
    
    mats = glob(fullfile(outdir, sub, sprintf('%s_*-*_interChRsq.mat', sub)));
    
    for jj = 1:length(mats)
        
        % note which site is being added currently
        eles = split(mats{jj}, '_');
        load(mats{jj});
        M = readmatrix(fullfile(outdir, sub, sprintf('%s_%s_nCAR.txt', sub, eles{2})));
        
        % 3rd dim of zs == number of recording channels + 1 (corresponding to no reference)
        assert(size(rsqs, 3) == M(2)+1, 'Mismatch between number of channels in CARLA output and correlation output'); 
        
        sites = [sites; {sub}, eles(2), {rsqs}, {M}];
    end
end

assert(size(sites, 1) == 82, 'Error: expecting 82 total stim sites across the 4 subjects');

%% Summarize mean correlation for no re-referencing, 50% covariance, CARLA, and full rereferencing

func = (@(x) mean(x, 'all', 'omitnan'));

% each row = stim site. Order of columns = no rereferencing, 50% rereferencing, full rereferencing, CARLA
rsqsAll = nan(size(sites, 1), 5);
for ii = 1:size(sites, 1)
    rsqs = sites{ii, 3};
    nChs = sites{ii, 4}(2); % number of recording channels
    nCarla = sites{ii, 4}(1); % optimal carla size
    
    rsqsAll(ii, 1) = func(rsqs(:, :, 1)); % mean correlation without referencing
    rsqsAll(ii, 2) = func(rsqs(:, :, end)); % mean correlation at full CAR
    rsqsAll(ii, 3) = func(rsqs(:, :, floor(nChs*0.25) + 1)); % mean correlation at 25%
    rsqsAll(ii, 4) = func(rsqs(:, :, floor(nChs*0.5) + 1)); % mean correlation at 50%
    rsqsAll(ii, 5) = func(rsqs(:, :, nCarla + 1)); % mean correlation using CARLA optimum
end

figure('Position', [200, 200, 200, 300]);
boxplot(rsqsAll, 'Symbol', '.');
xlim([0, 6]); ylim([0, 1]);
saveas(gcf, fullfile(outdir, 'rsCompsAllsubs'), 'png');
saveas(gcf, fullfile(outdir, 'rsCompsAllsubs'), 'svg');

ps = nan(4, 1);
for ii = 1:4
    ps(ii) = signrank(rsqsAll(:, ii), rsqsAll(:, 5));
end
psBonf = ps*4;
disp(psBonf);

%% Plot histogram of pair-wise stim site differences between each reference condition and CARLA

rsqsDiff = rsqsAll - rsqsAll(:, 5);
figure('Position', [200, 200, 200, 300]);
boxplot(rsqsDiff(:, 1:4), 'Symbol', '.');
yline(0, 'Color', [0.2, 0.2, 0.2]);
xlim([0, 5]); ylim([-0.4, 0.8]);
saveas(gcf, fullfile(outdir, 'rsDiffAllsubs'), 'png');
saveas(gcf, fullfile(outdir, 'rsDiffAllsubs'), 'svg');

%% Plot matrices for sub1 stimsite 1

sub = 'sub-1';
site = 'RMO8-RMO9';

idx = find(strcmp(sites(:, 1), sub) & strcmp(sites(:, 2), site));
rmat = sites{idx, 3};
nChs = sites{idx, 4}(2);
nCAR = sites{idx, 4}(1);

% no rereferencing, standard CAR, naive low variance CAR, CARLA
labels = {'none', 'CAR', 'bottom25', 'bottom50', 'CARLA'};
rmatSelect = rmat(:, :, [1, nChs+1, floor(nChs*0.25)+1, floor(nChs*0.5)+1, nCAR+1]);

for ii = 1:5
    
    figure;
    imagesc(rmatSelect(:, :, ii), [0, 1]);
    axis square;
    colormap(hot);
    colorbar;
    
    saveas(gcf, fullfile(outdir, sprintf('rsMap_%s_%s-%s', labels{ii}, sub, site)), 'png');
    saveas(gcf, fullfile(outdir, sprintf('rsMap_%s_%s-%s', labels{ii}, sub, site)), 'svg');
    
end
