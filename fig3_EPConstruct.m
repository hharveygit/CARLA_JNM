%% This script plots an example of enveloped sinusoids that add to form a simulated evoked potential (Figure 3)
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
%% Configure signal parameters

rng(4);

srate = 4800;
tt = -0.5:1/srate:1-1/srate;

A = 100;

tau1 = 0.01 + rand*0.02; % time constant to dampen the faster sinusoid (s1)
tau2 = 0.06 + rand*0.08;
tauOn1 = 0.005; % time constant of subtracted exponential decay for both
tauOn2 = 0.025; % time constant of subtracted exponential decay for both

f1 = 8 + rand*4; % freq (hz) of first sinusoid (mean period = 100 ms)
f2 = 1 + rand*2;  % freq of second sinusoid, (mean period = 500 ms)

ph1 = rand*2*pi; % phase of first sinusoid
ph2 = rand*2*pi; % phase of second sinusoid

outdir = fullfile('output', 'EPConstruct');
mkdir(outdir);

%% Construct and save full signal (3C)

sig = A*( (exp(-tt/tau1) - exp(-tt/tauOn1)).*sin(2*pi*f1*tt - ph1) + (exp(-tt/tau2) - exp(-tt/tauOn2)).*sin(2*pi*f2*tt - ph2) );
sig(tt < 0) = 0; % make causal

figure('Position', [200, 200, 200, 100]);
plot(tt, sig, 'k-');
xlabel('Time (s)');
xlim([-0.1, 0.5]);
ylim([-80, 80]);

saveas(gcf, fullfile(outdir, 'signal'), 'png');
saveas(gcf, fullfile(outdir, 'signal'), 'svg');

%% Plot faster part of signal (3A)

% the direct/faster signal component
env1 = (exp(-tt/tau1) - exp(-tt/tauOn1));
sig1 = A*sin(2*pi*f1*tt - ph1);
env1(tt < 0) = 0; % make causal
sig1(tt < 0) = 0;

figure('Position', [200, 200, 200, 100]);
plot(tt, env1, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
xline(tau1, 'Color', 'r'); xline(tauOn1, 'Color', 'r');
xlim([-0.1, 0.5]);
ylim([-1, 1]);
saveas(gcf, fullfile(outdir, 'env1'), 'png');
saveas(gcf, fullfile(outdir, 'env1'), 'svg');

figure('Position', [200, 200, 200, 100]);
plot(tt, sig1, 'k-');
xlim([-0.1, 0.5]);
ylim([-120, 120]);
saveas(gcf, fullfile(outdir, 'sig1'), 'png');
saveas(gcf, fullfile(outdir, 'sig1'), 'svg');

figure('Position', [200, 200, 200, 100]);
plot(tt, env1.*sig1, 'k-');
xlim([-0.1, 0.5]);
ylim([-80, 80]);
saveas(gcf, fullfile(outdir, 'env1sig1'), 'png');
saveas(gcf, fullfile(outdir, 'env1sig1'), 'svg');

%% Plot slower part of signal (3B)

env2 = (exp(-tt/tau2) - exp(-tt/tauOn2));
sig2 = A*sin(2*pi*f2*tt - ph2);
env2(tt < 0) = 0; % make causal
sig2(tt < 0) = 0;

figure('Position', [200, 200, 200, 100]);
plot(tt, env2, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
xlim([-0.1, 0.5]);
ylim([-1, 1]);
xline(tau2, 'Color', 'r'); xline(tauOn2, 'Color', 'r');
saveas(gcf, fullfile(outdir, 'env2'), 'png');
saveas(gcf, fullfile(outdir, 'env2'), 'svg');

figure('Position', [200, 200, 200, 100]);
plot(tt, sig2, 'k-');
xlim([-0.1, 0.5]);
ylim([-120, 120]);
saveas(gcf, fullfile(outdir, 'sig2'), 'png');
saveas(gcf, fullfile(outdir, 'sig2'), 'svg');

figure('Position', [200, 200, 200, 100]);
plot(tt, env2.*sig2, 'k-');
xlim([-0.1, 0.5]);
ylim([-80, 80]);
saveas(gcf, fullfile(outdir, 'env2sig2'), 'png');
saveas(gcf, fullfile(outdir, 'env2sig2'), 'svg');

%% Save full signal (again), using the constructed components

figure('Position', [200, 200, 200, 100]);
plot(tt, env1.*sig1 + env2.*sig2, 'k-');
xlabel('Time (s)');
xlim([-0.1, 0.5]);
ylim([-80, 80]);
