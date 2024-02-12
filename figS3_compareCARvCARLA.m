%% Compare CARLA vs CAR re-referenced signals on same set of channels. Run this script as an addendum after running applyCARLARealCCEPs
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

%% Sort the CARLA output and calculate the standard CAR-referenced signals

% sort the CARLA output channels by CARLA order
VcarlaSorted = Vcar(stats.order, :, :);

% calculate the standard common average and subtract from signals
CARstandard = mean(dataSiteNoNan, 1, 'omitnan');
VcarStandard = dataSiteNoNan - CARstandard;
VcarStandardSorted = VcarStandard(stats.order, :, :);

%% Zoom into the first 50 channels, compare side by side CAR vs CARLA

yspace = 50;
numChsToPlot = 50;

figure('Position', [200, 200, 1000, 1400]);

% the CARLA channels
subplot(1, 2, 1);
ieeg_plotTrials(tt, mean(VcarlaSorted(1:numChsToPlot, :, :), 3)', yspace, 1:numChsToPlot, [0, 0, 0]);
xlim([-0.1, 0.5]);

% the standard CAR channels
subplot(1, 2, 2);
ieeg_plotTrials(tt, mean(VcarStandardSorted(1:numChsToPlot, :, :), 3)', yspace, 1:numChsToPlot, [0, 0, 0]);
xlim([-0.1, 0.5]);

saveas(gcf, fullfile(outdir, sprintf('%s_%s_sortedChsCARvCARLA_1to%d', sub, site, numChsToPlot)), 'png');
print(gcf, fullfile(outdir, sprintf('%s_%s_sortedChsCARvCARLA_1to%d', sub, site, numChsToPlot)),'-depsc2', '-r300', '-painters');
