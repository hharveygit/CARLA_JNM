%% This script takes outputs from figS2a_simCCEPSNRs, and summarizes the accuracy for each level of responsiveness and SNR
% Generates outputs for Figure S2B
%
%   2024/02/12
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
%% Configure parameters and directories

parCm = brighten(parula(21), -0.5);
parCm(1, :) = [0, 0, 0]; % set first row to black for extra contrast
parCm(4, :) = [0.4940 0.1840 0.5560]; % purple

rootdir = 'simLoopSNR';
nChs = 50; % total number of channels
nReps = 30;
noiseCoefs = 0.4:0.1:1.0; % which noise coefficient levels to load
nResponsives = 0:5:40; % number of responsive channels, from 0 to 80% by 10% increments

%% Read and store accuracy values at each subset of responsiveness, save to file

% each field is in order of forloop in simCCEP_testSNRs: noise x responsiveness x nReps
accuracy = struct();
accuracy.TP = nan(length(noiseCoefs), length(nResponsives), nReps);
accuracy.TN = nan(length(noiseCoefs), length(nResponsives), nReps);
accuracy.FN = nan(length(noiseCoefs), length(nResponsives), nReps);
accuracy.FP = nan(length(noiseCoefs), length(nResponsives), nReps);

% similar structure for snrs but cell to accommodate differeing amounts
snrs = cell(length(noiseCoefs), length(nResponsives), nReps);

for nC = 1:length(noiseCoefs)
    for nR = 1:length(nResponsives)
        
        fprintf('.');
        
        indir = fullfile('output', rootdir, sprintf('noise%0.1f_nchs%d-%d', noiseCoefs(nC), nChs, nResponsives(nR)));
    
        for rep = 1:nReps
    
            T_acc = readtable(fullfile(indir, sprintf('accuracy_rep%d.txt', rep)), 'FileType', 'text', 'Delimiter', '\t');
            accuracy.TP(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TP'));
            accuracy.TN(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TN'));
            accuracy.FN(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FN'));
            accuracy.FP(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FP'));
    
            temp = load(fullfile(indir, sprintf('snr_rep%d.mat', rep)));
            snrs{nC, nR, rep} = temp.snr;
    
        end
        
    end
end
fprintf('\n');
save(fullfile('output', rootdir, 'accuracies_snrs_all.mat'), 'accuracy', 'snrs');

%% Plot FN, FP vs resp traces for each noise level

load(fullfile('output', rootdir, 'accuracies_snrs_all.mat'), 'accuracy', 'snrs');

cmSens = [1, 165/255, 0]; % orange = sensitive threshold

% Based on more sensitive, first-peak optimum
figure('Position', [200, 200, 600, 900]);
subplot(2, 1, 1); % FN = responsive channels included in the CAR
hold on
for nc = 1:length(noiseCoefs)
    FN = squeeze(accuracy.FN(nc, :, :));
    iqr = prctile(FN, [25, 75], 2);
    med = median(FN, 2);
    plot(nResponsives, med, '-o', 'Color', parCm(nc*3-2, :), 'LineWidth', 1);
    errorbar(nResponsives, med, med - iqr(:, 1), iqr(:, 2) - med, 'Color', parCm(nc*3-2, :), 'LineWidth', 1, 'CapSize', 8);
end
hold off
xlim([-1, 41]); ylim([-1, 50]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50, 'box', 'off');
ylabel('Misses');
title('Responsive Channels Missed');

subplot(2, 1, 2);
hold on
for nc = 1:length(noiseCoefs)
    FP = squeeze(accuracy.FP(nc, :, :));
    iqr = prctile(FP, [25, 75], 2);
    med = median(FP, 2);
    plot(nResponsives, med, '-o', 'Color', parCm(nc*3-2, :), 'LineWidth', 1);
    errorbar(nResponsives, med, med - iqr(:, 1), iqr(:, 2) - med, 'Color', parCm(nc*3-2, :), 'LineWidth', 1, 'CapSize', 8);
end
hold off
xlim([-1, 41]); ylim([-1, 50]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50, 'box', 'off');
xlabel('Number of responsive channels / 50');
ylabel('Misses');
title('Nonresponsive Channels Missed');
saveas(gcf, fullfile('output', rootdir, sprintf('nchs%d_FNFP', nChs)), 'png');
saveas(gcf, fullfile('output', rootdir, sprintf('nchs%d_FNFP', nChs)), 'svg');

% Calculate mean SNR for each noise level (geo bc lognormal and arith)
snrGeomean = nan(length(noiseCoefs), 1);
snrGeoSD = nan(length(noiseCoefs), 1);
snrMean = nan(length(noiseCoefs), 1);
snrSD = nan(length(noiseCoefs), 1);
for ii = 1:length(noiseCoefs)
    snrsCurrCell = reshape(snrs(ii, :, :), 1, []);
    snrsCurr = [];
    for jj = 1:length(snrsCurrCell)
        snrsCurr = [snrsCurr; snrsCurrCell{jj}(:)];
    end
    snrsCurr(isnan(snrsCurr)) = [];
    snrGeomean(ii) = geomean(snrsCurr);
    snrGeoSD(ii) = 10^std(log10(snrsCurr));
    snrMean(ii) = mean(snrsCurr);
    snrSD(ii) = std(snrsCurr);
end

