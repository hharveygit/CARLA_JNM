%% This script takes outputs from figS1a_simCCEPTotalChs, and summarizes the accuracy for each total channel size and level of responsiveness
% Generates outputs for Figure S1A
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

parCm = brighten(parula(15), -0.5); % black is the 50 ch condition, yellow is lowest number of channels (10)
parCm(1, :) = [0, 0, 0];
parCm(4, :) = [0.4940 0.1840 0.5560];

% Configure the name of the folders to load accuracy scores from.
rootdir = 'simLoopNChs';
rootdir50 = 'simLoopNoGlob'; % We do not regenerate the 50 Chs condition, just pull from previous outputs

nChsToTest = [50, 25, 20, 15, 10]; % total number of channels
fracResp = 0:0.2:0.8; % fraction of all channels responsive
nReps = 30;

%% Read and store accuracy values at each number of total channels, then save to file

% each field is in order of forloop in simCCEP_testNChs: nchs x responsiveness x nReps
accuracy = struct();
accuracy.TP = nan(length(nChsToTest), length(fracResp), nReps);
accuracy.TN = nan(length(nChsToTest), length(fracResp), nReps);
accuracy.FN = nan(length(nChsToTest), length(fracResp), nReps);
accuracy.FP = nan(length(nChsToTest), length(fracResp), nReps);
accuracy.TPglob = nan(length(nChsToTest), length(fracResp), nReps); % for the global maxes
accuracy.TNglob = nan(length(nChsToTest), length(fracResp), nReps);
accuracy.FNglob = nan(length(nChsToTest), length(fracResp), nReps);
accuracy.FPglob = nan(length(nChsToTest), length(fracResp), nReps);

% First store the 50 ch condition by pulling from previous output folder
for nR = 1:length(fracResp)
    fprintf('.');
    nResp = round(fracResp(nR)*50);
    indir = fullfile('output', rootdir50, sprintf('nchs50-%d', nResp));
    for rep = 1:nReps

        T_acc = readtable(fullfile(indir, sprintf('accuracy_50-%d_rep%d.txt', nResp, rep)), 'FileType', 'text', 'Delimiter', '\t');
        accuracy.TP(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TP'));
        accuracy.TN(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TN'));
        accuracy.FN(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FN'));
        accuracy.FP(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FP'));
        accuracy.TPglob(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TPglob'));
        accuracy.TNglob(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TNglob'));
        accuracy.FNglob(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FNglob'));
        accuracy.FPglob(1, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FPglob'));

    end
end

% now load for each total number of channels from rootdir
for nC = 2:length(nChsToTest)
    for nR = 1:length(fracResp)
        fprintf('.');
        nChs = nChsToTest(nC);
        indir = fullfile('output', rootdir, sprintf('nchs%d-%d', nChs, round(fracResp(nR)*nChs)));
        for rep = 1:nReps
            T_acc = readtable(fullfile(indir, sprintf('accuracy_rep%d.txt', rep)), 'FileType', 'text', 'Delimiter', '\t');
            accuracy.TP(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TP'));
            accuracy.TN(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TN'));
            accuracy.FN(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FN'));
            accuracy.FP(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FP'));
            accuracy.TPglob(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TPglob'));
            accuracy.TNglob(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'TNglob'));
            accuracy.FNglob(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FNglob'));
            accuracy.FPglob(nC, nR, rep) = T_acc.Var2(strcmp(T_acc.Var1, 'FPglob'));
        end
    end
end
fprintf('\n');
save(fullfile('output', rootdir, 'accuracies_all.mat'), 'accuracy');

%% Plot FN, FP vs resp traces for each noise level, for first-peak method, as bars

load(fullfile('output', rootdir, 'accuracies_all.mat'), 'accuracy');

% displacements for the bars
disps = [-0.04, -0.02, 0, 0.02, 0.04];

% Based on more sensitive, first-peak optimum
figure('Position', [200, 200, 700, 800]);
subplot(2, 1, 1); % FN = responsive channels included in the CAR
hold on
for nc = 1:size(accuracy.FN, 1) % starting at 2 excludes the 50
    FN = squeeze(accuracy.FN(nc, :, :));
    FN = FN ./ (fracResp'*nChsToTest(nc)); % normalize by number of responsive channels
    FN(isnan(FN)) = 0; % when there are 0 responsive channels, change these nans to 0
    iqr = prctile(FN, [25, 75], 2);
    med = median(FN, 2);
    bar(fracResp + disps(nc), med, 0.08, 'FaceColor', parCm(nc*3-2, :));
    errorbar(fracResp + disps(nc), med, med - iqr(:, 1), iqr(:, 2) - med, 'Color', parCm(nc*3-2, :), 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 8);
end
hold off
xlim([-0.05, 0.85]); ylim([-0.05, 1.05]);
set(gca, 'xtick', 0:0.2:0.8, 'xticklabels', 0:20:80, 'ytick', 0:0.2:1, 'yticklabels', 0:20:100, 'box', 'off');
ylabel('% Responsive Channels Missed');

subplot(2, 1, 2);
hold on
for nc = 1:size(accuracy.FN, 1)
    FP = squeeze(accuracy.FP(nc, :, :));
    FP = FP ./ ((1-fracResp)'*nChsToTest(nc)); % normalize by number of non-responsive channels
    iqr = prctile(FP, [25, 75], 2);
    med = median(FP, 2);
    bar(fracResp + disps(nc), med, 0.08, 'FaceColor', parCm(nc*3-2, :));
    errorbar(fracResp + disps(nc), med, med - iqr(:, 1), iqr(:, 2) - med, 'Color', parCm(nc*3-2, :), 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 8);
end
hold off
xlim([-0.05, 0.85]); ylim([-0.05, 1.05]);
set(gca, 'xtick', 0:0.2:0.8, 'xticklabels', 0:20:80, 'ytick', 0:0.2:1, 'yticklabels', 0:20:100, 'box', 'off');
xlabel('Fraction of Channels with Response');
ylabel('% Non-responsive Channels Missed');
saveas(gcf, fullfile('output', rootdir, 'FNFP_bars_acrossNChs'), 'png');
saveas(gcf, fullfile('output', rootdir, 'FNFP_bars_acrossNChs'), 'svg');
