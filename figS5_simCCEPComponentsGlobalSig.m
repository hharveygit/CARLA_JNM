%% Figure S2. Create simulated CCEP data and their components, like Figure 1, but this time containing a global signal
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

srate = 4800;
tt = -0.5:1/srate:1-1/srate;

nchs = 10; % number of simulated channels
ntrs = 12; % number of trials
chs = arrayfun(@(x) sprintf('ch%d', x), 1:nchs, 'UniformOutput', false)';

outdir = fullfile('output', 'simComponentsGlobal');
mkdir(outdir)

%% i) Create data, add stim artifacts

% stores the data
V0 = zeros(length(tt), nchs);

Aart = 50 + rand(nchs, 1)*6; % slightly different artifact amplitudes for each channel
artifact = sin(2*pi*600*tt)';
artifact(tt < 0 | tt > 0.002) = 0;
V0 = V0 + artifact*Aart';

figure('Position', [200, 200, 400, 800]); ieeg_plotTrials(tt, V0, 100);
xlabel('Time (s)'); ylabel('Channels');

%% A) Add unique input signals at a subset of channels

rng(14);

chsResp = 1:4;

A = 100;
% using previous version of signal generator for this illustration to preserve backwards compatibility
sig = genRandSigOrig(tt, length(chsResp), A);

V1 = V0;
V1(:, chsResp) = V0(:, chsResp) + sig;

figure('Position', [200, 200, 400, 800]); ieeg_plotTrials(tt, V1, 100, [], [], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Channels');

%% B) Add a global noise across all chs and trials, e.g. from reference.

rng(10);

Aglobal = 20; % same amplitude of noise on all channels

sigGlob = genRandSigOrig(tt, 1, Aglobal);

V2 = V1 + sigGlob;

figure('Position', [200, 200, 400, 800]); ieeg_plotTrials(tt, V2, 100, [], [], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Channels');

%% C) Add noise common to all channels (line noise and some brown noise from reference)

rng('default');

% add trials, so now put V in ch x time x trial format
V3 = repmat(V2', 1, 1, ntrs);

% random line noise phases for each trial, trials x harmonic (60, 120, 180)
phLN = rand(ntrs, 3)*2*pi;
LN = zeros(length(tt), ntrs);
for ii = 1:ntrs
    LN(:, ii) = 8*sin(2*pi*60*tt - phLN(ii, 1)) + 2*sin(2*pi*120*tt - phLN(ii, 2)) + 1*sin(2*pi*180*tt - phLN(ii, 3));
end

% brown noise shared across channels (from ref electrode
BN = cumsum(0.4*randn(2*length(tt), ntrs));
BN = ieeg_highpass(BN, srate, true);
BN = BN((0.5*length(tt)+1) : 1.5*length(tt), :);

noiseCommon = LN + BN;

V3 = V3 + shiftdim(noiseCommon, -1);

figure('Position', [200, 200, 600, 400]);
subplot(1, 2, 1); ieeg_plotTrials(tt, V3(:, :, 1)', 100, [], [], 'LineWidth', 1.5);
title('V at one trial');
xlabel('Time (s)'); ylabel('Channels');

subplot(1, 2, 2); ieeg_plotTrials(tt, mean(V3, 3)', 100, [], [], 'LineWidth', 1.5);
title('V across trials');
xlabel('Time (s)'); ylabel('Channels');


%% D) add random brown noise across all channels

rng('default');

noiseRand = cumsum(0.4*randn(nchs, 2*length(tt), ntrs), 2); % give double the number of time points so we can highpass it
for ii = 1:nchs
    noiseRand(ii, :, :) = ieeg_highpass(squeeze(noiseRand(ii, :, :)), srate, true);
end
noiseRand = noiseRand(:, (0.5*length(tt)+1) : 1.5*length(tt), :);

V4 = V3 + noiseRand;

figure('Position', [200, 200, 600, 400]);
subplot(1, 2, 1); ieeg_plotTrials(tt, V4(:, :, 1)', 100, [], [], 'LineWidth', 1.5);
title('V at one trial');
xlabel('Time (s)'); ylabel('Channels');

subplot(1, 2, 2); ieeg_plotTrials(tt, mean(V4, 3)', 100, [], [], 'LineWidth', 1.5);
title('V across trials');
xlabel('Time (s)'); ylabel('Channels');

saveas(gcf, fullfile(outdir, 'V4'), 'svg');
saveas(gcf, fullfile(outdir, 'V4'), 'png');


%% Save construction components for channels of interest

ylims = [-120, 120];

% Plot select channels for figure (3 and 8)
chs2Plot = [3, 8];

for ii = 1:length(chs2Plot)

    ch = chs2Plot(ii);

    % plot full simulation of all trials at data
    figure('Position', [200, 200, 600, 300]); hold on
    plot(tt, squeeze(V4(ch, :, :)), 'Color', [0.5, 0.5, 0.5]);
    plot(tt, mean(squeeze(V4(ch, :, :)), 2), 'k', 'LineWidth', 1.5);
    hold off
    ylim(ylims);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    saveas(gcf, fullfile(outdir, sprintf('V4_ch%d', ch)), 'png');
    saveas(gcf, fullfile(outdir, sprintf('V4_ch%d', ch)), 'svg');
    
    % plot the true signal
    figure('Position', [200, 200, 200, 300]);
    plot(tt, V1(:, ch), 'k', 'LineWidth', 2);
    ylim(ylims); xlim([0, 0.5]);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    saveas(gcf, fullfile(outdir, sprintf('V1_ch%d', ch)), 'png');
    saveas(gcf, fullfile(outdir, sprintf('V1_ch%d', ch)), 'svg');
    
    % plot the global noise
    figure('Position', [200, 200, 200, 300]);
    plot(tt, sigGlob, 'k', 'LineWidth', 2);
    ylim(ylims); xlim([0, 0.5]);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    saveas(gcf, fullfile(outdir, 'globalSig'), 'png');
    saveas(gcf, fullfile(outdir, 'globalNoise'), 'svg');
    
    % plot examples of the common noise
    figure('Position', [200, 200, 200, 300]);
    ieeg_plotTrials(tt, noiseCommon(:, 1:3), 50, [], [0.5, 0.5, 0.5]); % plot at one trial
    xlim([0, 0.5]);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    saveas(gcf, fullfile(outdir, 'lineNoise_3trs'), 'png');
    saveas(gcf, fullfile(outdir, 'lineNoise_3trs'), 'svg');
    
    % plot the random noise
    figure('Position', [200, 200, 200, 300]); hold on
    plot(tt, squeeze(noiseRand(ch, :, :)), 'Color', [0.5, 0.5, 0.5]);
    hold off
    xlim([0, 0.5]); ylim(ylims);
    xlabel('Time (s)'); ylabel('Voltage (\muV)');
    saveas(gcf, fullfile(outdir, sprintf('randNoise_ch%d', ch)), 'png');
    saveas(gcf, fullfile(outdir, sprintf('randNoise_ch%d', ch)), 'svg');

end
