%% Figure 2: Methods illustration of how CARLA works, on single-trial simulated CCEP data
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
%% A) Create data, add stim artifacts and responses to the last 2 channels

srate = 4800;
tt = -0.5:1/srate:1-1/srate;

rng(50);

nchs = 6; % number of simulated channels
ntrs = 1; % number of trials
chs = arrayfun(@(x) sprintf('ch%d', x), 1:nchs, 'UniformOutput', false)';

% stores the data
V0 = zeros(length(tt), nchs);

Aart = 50 + rand(nchs, 1)*5; % slightly different artifact amplitudes for each channel
artifact = sin(2*pi*600*tt)';
artifact(tt < 0 | tt > 0.002) = 0;
V0 = V0 + artifact*Aart';

chsResp = [nchs-1, nchs];

sig = genRandSig(tt, length(chsResp), 150);

V1 = V0;
V1(:, chsResp) = V0(:, chsResp) + sig;

figure('Position', [200, 200, 400, 800]); ieeg_plotTrials(tt, V1, 80, [], [], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Channels');


%% B) Add common noise (line noise and reference brown noise) to each trial

rng('default');

% no line noise in this demonstration, pretend that this has been removed by filtering
% Create line noise: random noise phase for each harmonic (60, 120, 180)
%phNoise = rand(1, 3)*2*pi;
%LN = sin(2*pi*60*tt - phNoise(1))' + 0.3*sin(2*pi*120*tt - phNoise(2))' + 0.1*sin(2*pi*180*tt - phNoise(3))';

% Create some brown noise and low-pass as common signal
brownNoiseRef = cumsum(0.4*randn(2*length(tt), 1));
brownNoiseRef = ieeg_highpass(brownNoiseRef, srate, true);
brownNoiseRef = brownNoiseRef((0.5*length(tt)+1) : 1.5*length(tt));

V3 = V1 + brownNoiseRef;

figure('Position', [200, 200, 400, 800]); ieeg_plotTrials(tt, V3, 80, [], [], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Channels');


%% C) add random brown noise to each channel

rng(4);

noiseRand = cumsum(0.4*randn(2*length(tt), nchs)); % give double the number of time points so we can highpass it
noiseRand = ieeg_highpass(noiseRand, srate, true);

noiseRand = noiseRand((0.5*length(tt)+1) : 1.5*length(tt), :);

V4 = V3 + noiseRand;

figure('Position', [200, 200, 400, 800]); ieeg_plotTrials(tt, V4, 80, [], [], 'LineWidth', 1.5);
xlim([-0.1, 0.5])
xlabel('Time (s)'); ylabel('Channels');

%% Apply CARLA and show outputs

rng('default');

[Vcar, CAR, stats] = CARLA(tt, V4, srate);

outdir = fullfile('output', 'simulationTrial');
mkdir(outdir);

% Plot in increasing order of variance (Figure 2A)
V4sorted = V4(:, stats.order);
figure('Position', [200, 200, 200, 500]); ieeg_plotTrials(tt, V4sorted, 100, [], [], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Channels');
xlim([-0.1, 0.5]);
xline(0.01, 'Color', 'r');
xline(0.3, 'Color', 'r');
saveas(gcf, fullfile(outdir, 'trialExample'), 'svg');
saveas(gcf, fullfile(outdir, 'trialExample'), 'png');

% Plot zmin (Figure 2C)
figure('Position', [200, 200, 400, 300]); hold on
plot(stats.zMinMean, 'k-o', 'LineWidth', 1, 'MarkerFaceColor', 'k');
xlim([1, 6]); ylim([-1.5, 0]);
saveas(gcf, fullfile(outdir, 'trialExampleRmin'), 'svg');
saveas(gcf, fullfile(outdir, 'trialExampleRmin'), 'png');

% Plot variance in increasing order (Figure 2D)
vars = stats.vars(stats.order);
figure('Position', [200, 200, 400, 250]); hold on
plot(vars, 'k-o', 'LineWidth', 1, 'MarkerFaceColor', 'k');
xlim([1, 6]);
saveas(gcf, fullfile(outdir, 'trialExampleVar'), 'svg');
saveas(gcf, fullfile(outdir, 'trialExampleVar'), 'png');

%% Plot CAR-rereferenced signal at a few sizes

% how many channels to use for common average to illustrate re-referenced signals (Figure 2B)
ns = [2, 4, 5];

for ii = 1:length(ns)

    n = ns(ii);

    CARcurr = mean(V4sorted(:, 1:n), 2);
    V4Carn = V4sorted - mean(V4sorted(:, 1:n), 2);
    
    figure('Position', [200, 200, 150, 500]); ieeg_plotTrials(tt, V4Carn, 100, [], [], 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Channels');
    xlim([0.01, 0.3]);
    saveas(gcf, fullfile(outdir, sprintf('trialExampleCAR%dch', n)), 'svg');
    saveas(gcf, fullfile(outdir, sprintf('trialExampleCAR%dch', n)), 'png');

end

