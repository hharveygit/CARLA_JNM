%% This script takes outputs from simCCEPLoop, and evaluates accuracy at each level of responsiveness.
% Generates outputs for Figure 4C,D and 6B
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

%% Configure loading directories
%
nchs = 50; % total number of channels

% rename the output folder depending on whether Aglobal was 0 to avoid overwriting Figure 4 outputs with Figure 6 outputs
if ~Aglobal % no global signal
    rootdir = fullfile('output', 'simLoopNoGlob');
else % yes global signal
    rootdir = fullfile('output', 'simLoopWithGlob');
end
movefile(fullfile('output', 'simLoop'), rootdir);

%% Read and store accuracy values at each subset of responsiveness

% each field is reps x nresp (starting at 0). XXglob = at cutoff found using global maximum
accuracy = struct();
accuracy.TP = nan(30, nchs+1);
accuracy.TN = nan(30, nchs+1);
accuracy.FN = nan(30, nchs+1);
accuracy.FP = nan(30, nchs+1);
accuracy.TPglob = nan(30, nchs+1);
accuracy.TNglob = nan(30, nchs+1);
accuracy.FNglob = nan(30, nchs+1);
accuracy.FPglob = nan(30, nchs+1);

for nresp = 0:nchs
    
    fprintf('.');
    
    indir = fullfile(rootdir, sprintf('nchs%d-%d', nchs, nresp));

    for rep = 1:30

        T_acc = readtable(fullfile(indir, sprintf('accuracy_%d-%d_rep%d.txt', nchs, nresp, rep)), 'FileType', 'text', 'Delimiter', '\t');
        accuracy.TP(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'TP'));
        accuracy.TN(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'TN'));
        accuracy.FN(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'FN'));
        accuracy.FP(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'FP'));
        accuracy.TPglob(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'TPglob'));
        accuracy.TNglob(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'TNglob'));
        accuracy.FNglob(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'FNglob'));
        accuracy.FPglob(rep, nresp+1) = T_acc.Var2(strcmp(T_acc.Var1, 'FPglob'));

    end
    
end
fprintf('\n');

%% Identify examples that correspond to median values at each # responsive channels

FNmed = median(accuracy.FN);
FPmed = median(accuracy.FP);
FNglobMed = median(accuracy.FNglob);
FPglobMed = median(accuracy.FPglob);

if contains(rootdir, 'simLoopWithGlob')
    idxes = 26; % 50% for main figure with global signal
    FNmedSelect = FNmed(idxes);
    FPmedSelect = FPmed(idxes);
    FNglobmedSelect = FNglobMed(idxes);
    FPglobmedSelect = FPglobMed(idxes);
    fprintf('Median FN at 25 resp chs = %0.1f\n', FNmedSelect);
    fprintf('Median FP at 25 = %0.1f\n', FPmedSelect);
    fprintf('Median FN (global) at 25 = %0.1f\n', FNglobmedSelect);
    fprintf('Median FP (global) at 25 = %0.1f\n', FPglobmedSelect);

elseif contains(rootdir, 'simLoopNoGlob')
    idxes = [1, 11, 21, 31, 41, 46]; % every 20%, for main figure with no global signal
    FNmedSelect = FNmed(idxes);
    FPmedSelect = FPmed(idxes);
    FNglobmedSelect = FNglobMed(idxes);
    FPglobmedSelect = FPglobMed(idxes);
    fprintf('Median FN at 0, 10, 20, 30, 40, 45 resp chs = %0.1f, %0.1f, %0.1f, %0.1f, %0.1f, %0.1f\n', FNmedSelect(:));
    fprintf('Median FP at 0, 10, 20, 30, 40, 45 = %0.1f, %0.1f, %0.1f, %0.1f, %0.1f, %0.1f\n', FPmedSelect(:));
    fprintf('Median FN (global) at 0, 10, 20, 30, 40, 45 = %0.1f, %0.1f, %0.1f, %0.1f, %0.1f, %0.1f\n', FNglobmedSelect(:));
    fprintf('Median FP (global) at 0, 10, 20, 30, 40, 45 = %0.1f, %0.1f, %0.1f, %0.1f, %0.1f, %0.1f\n', FPglobmedSelect(:));

end

% which rep to use that best matches the median
repSelect = zeros(size(idxes));
for ii = 1:length(idxes)
    rep = find(accuracy.FN(:, idxes(ii)) == floor(FNmedSelect(ii)) & ...
               accuracy.FP(:, idxes(ii)) == floor(FPmedSelect(ii)) & ...
               accuracy.FNglob(:, idxes(ii)) == floor(FNglobmedSelect(ii)) & ...
               accuracy.FPglob(:, idxes(ii)) == floor(FPglobmedSelect(ii)), 1, 'first');
    if isempty(rep), rep = nan; end
    repSelect(ii) = rep;
end
fprintf('Rep to use at responsiveness: %d \n', repSelect(:));

%% Absolute values of FP, FN

cmSens = [1, 165/255, 0]; % orange = sensitive threshold

% Using global optimum
% up to 90% responsiveness (45/50 responsive channels)
figure('Position', [200, 200, 600, 600]);
subplot(2, 1, 1); % FN = responsive channels included in the CAR
boxplot(accuracy.FNglob(:, 1:46), 'Positions', 0:45, 'Colors', 'b', 'PlotStyle', 'compact');
xlim([-1, 46]); ylim([-1, 50]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50, 'box', 'off');
ylabel('Number of incorrectly assigned channels');
title('False Negative (Resp. Chs in CAR)');

subplot(2, 1, 2);
boxplot(accuracy.FPglob(:, 1:46), 'Positions', 0:45, 'Colors', 'b', 'PlotStyle', 'compact');
xlim([-1, 46]); ylim([-1, 50]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50, 'box', 'off');
xlabel('Number of responsive channels / 50');
ylabel('Number of incorrectly assigned channels');
title('False Positive (Non-Resp. Chs in CAR)');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_FNFPglob', nchs)), 'png');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_FNFPglob', nchs)), 'svg');

% For sensitive method
% up to 90% responsiveness (45/50 responsive channels)
figure('Position', [200, 200, 600, 600]);
subplot(2, 1, 1); % FN = responsive channels included in the CAR
boxplot(accuracy.FN(:, 1:46), 'Positions', 0:45, 'Colors', cmSens, 'PlotStyle', 'compact');
xlim([-1, 46]); ylim([-1, 50]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50, 'box', 'off');
ylabel('Number of incorrectly assigned channels');
title('False Negative (Resp. Chs in CAR)');

subplot(2, 1, 2);
boxplot(accuracy.FP(:, 1:46), 'Positions', 0:45, 'Colors', cmSens, 'PlotStyle', 'compact');
xlim([-1, 46]); ylim([-1, 50]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50, 'box', 'off');
xlabel('Number of responsive channels / 50');
ylabel('Number of incorrectly assigned channels');
title('False Positive (Non-Resp. Chs in CAR)');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_FNFP', nchs)), 'png');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_FNFP', nchs)), 'svg');

%% Sensitivity specificity plot

% sensitivity and specificity
sensglob = accuracy.TPglob ./ (accuracy.TPglob + accuracy.FNglob);
specglob = accuracy.TNglob ./ (accuracy.TNglob + accuracy.FPglob);
sens = accuracy.TP ./ (accuracy.TP + accuracy.FN); % more important. need to sensitively detect responsive channels
spec = accuracy.TN ./ (accuracy.TN + accuracy.FP);

% balanced accuracy
accglob = (sensglob + specglob)/2;
acc = (sens + spec)/2;


% Using global optimum
cmBlue = [0, 0, 139; 115, 147, 179]/255;
figure('Position', [200, 200, 600, 300]); hold on
errorbar(0:45, mean(sensglob(:, 1:46)), std(sensglob(:, 1:46)), '-o', 'MarkerSize', 5, 'MarkerFaceColor', cmBlue(1, :), 'Color', cmBlue(1, :));
errorbar(0:45, mean(specglob(:, 1:46)), std(specglob(:, 1:46)), '--s', 'MarkerSize', 6, 'Color', cmBlue(2, :));
xlim([-1, 46]); ylim([-0.2, 1.2]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50);
xlabel('Number of responsive channels / 50');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_sensSpecglob', nchs)), 'png');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_sensSpecglob', nchs)), 'svg');

% Using sensitive method
cmOrange = [139, 64, 0; 255, 172, 28]/255;
figure('Position', [200, 200, 600, 300]); hold on
errorbar(0:45, mean(sens(:, 1:46)), std(sens(:, 1:46)), '-o', 'MarkerSize', 5, 'MarkerFaceColor', cmOrange(1, :), 'Color', cmOrange(1, :));
errorbar(0:45, mean(spec(:, 1:46)), std(spec(:, 1:46)), '--s', 'MarkerSize', 6, 'Color', cmOrange(2, :));
xlim([-1, 46]); ylim([-0.2, 1.2]);
set(gca, 'xtick', 0:5:50, 'xticklabels', 0:5:50);
xlabel('Number of responsive channels / 50');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_sensSpec', nchs)), 'png');
saveas(gcf, fullfile(rootdir, sprintf('nchs%d_sensSpec', nchs)), 'svg');

