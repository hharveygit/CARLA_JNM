%% This script is derived from simCCEPLoop, except an outer loop varies the amplitude of brown noise (mulitiplicative coefficient)
% We simulate 30 reps for each combination of SNR and responsiveness and test CARLA's accuracy
% Generate data with responses at 0, 10, 20, 30, 40, 50, 60, 70, 80 percent of all channels
% The last section of this script saves examples of how the same signal looks at varying levels of SNR (Figure S2A)
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
%% Control parameters of simulated data

srate = 4800;
tt = -0.5:1/srate:1-1/srate;

nchs = 50; % number of simulated channels
ntrs = 12; % number of trials
reps = 30; % number of reps

cmSens = [1, 165/255, 0]; % use orange as color for the sensitive optimum

Aglobal = 0;

% window to calculate SNR on. Use same window as for calculating responses
snrWin = [0.01, 0.3];

%% Outer loop through brown noise level

% Set seed depending on where brown coef range is starting, for reproducibility when parallelizing.
% This is configured in main
rng(brownCoefRange(1));

for brownCoef = brownCoefRange(1):0.1:brownCoefRange(2)

    fprintf('Brown noise coefficient = %0.1f\n', brownCoef);

    %% Inner loop through number of responsive channels
    for nresp = 0:5:40 % 0 to 80% responsiveness
        
        %% Create data
        
        fprintf('%d of %d channels responsive\n', nresp, nchs);
        
        outdir = fullfile('output', 'simLoopSNR', sprintf('noise%0.1f_nchs%d-%d', brownCoef, nchs, nresp));
        mkdir(outdir);
    
        chsResp = 1:nresp; % first x channels responsive, for ease. order doesn't matter. will sort and color them red
    
        %%
        for rr = 1:reps % multiple repetitions at each significance
            
            % i) artifact only
            V0 = zeros(length(tt), nchs);
            Aart = 50 + rand(nchs, 1)*5; % slightly different artifact amplitudes for each channel
            artifact = sin(2*pi*600*tt)';
            artifact(tt < 0 | tt > 0.002) = 0;
            V0 = V0 + artifact*Aart';
    
            % A) Add the evoked potentials
            A = 100;
            V1 = V0;
            sig = genRandSig(tt, length(chsResp), A);
            V1(:, chsResp) = V0(:, chsResp) + sig;
    
            % B) Option to add a global noise
            if Aglobal
                sigCommon = genRandSig(tt, 1, Aglobal);
            else
                sigCommon = 0;
            end
            V2 = V1 + sigCommon;
    
            % C) Add common noise to all channels at each trial
            V3 = repmat(V2', 1, 1, ntrs); % ch x time points x trial
            phLN = rand(ntrs, 3)*2*pi; % LN phases
            LN = zeros(length(tt), ntrs);
            for ii = 1:ntrs
                LN(:, ii) = 8*sin(2*pi*60*tt - phLN(ii, 1)) + 2*sin(2*pi*120*tt - phLN(ii, 2)) + 1*sin(2*pi*180*tt - phLN(ii, 3));
            end
    
            % brown noise shared across channels
            BN = cumsum(brownCoef*randn(2*length(tt), ntrs)); % variable brown noise coefficient now
            BN = ieeg_highpass(BN, srate, true);
            BN = BN((0.5*length(tt)+1) : 1.5*length(tt), :);
    
            noiseCommon = LN + BN;
            V3 = V3 + shiftdim(noiseCommon, -1);
    
            % D) add random brown noise
            noiseRand = cumsum(brownCoef*randn(nchs, 2*length(tt), ntrs), 2); % give double the number of time points so we can highpass it
            for ii = 1:nchs
                noiseRand(ii, :, :) = ieeg_highpass(squeeze(noiseRand(ii, :, :)), srate, true);
            end
            noiseRand = noiseRand(:, (0.5*length(tt)+1) : 1.5*length(tt), :);
            V4 = V3 + noiseRand;
    
            %% Calculate SNR
    
            if nresp == 0
                Psig = nan;
                PnoiseRand = squeeze(sum(noiseRand(chsResp, tt >= snrWin(1) & tt < snrWin(2), :).^2, 2) / diff(snrWin));
                snr = nan;
            else
        
                % power of signal (for each responsive channel created)
                Psig = sum(sig(tt >= snrWin(1) & tt < snrWin(2), :).^2)' / diff(snrWin); % express as per second
        
                % power of the aperiodic noise (common brown noise + random noise), calculated separately for each trial
                noiseSum = noiseRand + shiftdim(BN, -1);
                PnoiseRand = squeeze(sum(noiseSum(chsResp, tt >= snrWin(1) & tt < snrWin(2), :).^2, 2) / diff(snrWin));
    
                % calculate snr for each trial separately (same Psig for all trials at one channel)
                snr = repmat(Psig, 1, ntrs) ./ PnoiseRand;
        
                % for simplicity, we only vary and consider random noise in SNR, assuming that periodic noise can be mostly attenuated by filtering
    
            end
            fprintf('%0.2f ', geomean(snr, 'all')); % geometric mean since power is lognormal
    
            %% Apply CARLA and plot outputs
            
            [Vout, CAR, stats] = CARLA(tt, V4, srate, true); % get the sensitive output
    
            % number of channels used for the CAR
            nCAR = length(stats.chsUsed);
            [~, nCARglob] = max(mean(stats.zMinMean, 2)); % number of channels at global maximum
    
            % Plot average zmin across trials
            figure('Position', [200, 200, 400, 300], 'Visible', 'off'); hold on
            errorbar(mean(stats.zMinMean, 2), std(stats.zMinMean, 0, 2), 'k-o');
            plot(nCARglob, mean(stats.zMinMean(nCARglob, :), 2), 'b*'); % global max as blue
            if nCARglob ~= nCAR; plot(nCAR, mean(stats.zMinMean(nCAR, :), 2), '*', 'Color', cmSens); end
            yline(0, 'Color', 'k');
            saveas(gcf, fullfile(outdir, sprintf('zmin_rep%d', rr)), 'png');
            saveas(gcf, fullfile(outdir, sprintf('zmin_rep%d', rr)), 'svg');
    
            % Sort and plot channels by increasing covariance, draw line at cutoff
            V4MeanSorted = mean(V4(stats.order, :, :), 3);
            respBool = antifind(chsResp, nchs);
            respBool = respBool(stats.order); % logical array of where responsive channels are
            cm = zeros(nchs, 3);
            cm(respBool, 1) = 1; % make red
            figure('Position', [200, 200, 250, 600], 'Visible', 'off');
            yspace = 80;
            ys = ieeg_plotTrials(tt, V4MeanSorted', yspace, [], cm, 'LineWidth', 1);
            yline(ys(nCARglob)-yspace/2, 'Color', 'b', 'LineWidth', 1.5);
            if nCARglob ~= nCAR; yline(ys(nCAR)-yspace/2, 'Color', cmSens, 'LineWidth', 1.5); end
            xlim([-0.1, 0.5]); set(gca, 'xtick', [0, 0.5]);
            xlabel('Time (s)'); ylabel('Channels');
            saveas(gcf, fullfile(outdir, sprintf('chsSorted_rep%d', rr)), 'png');
            saveas(gcf, fullfile(outdir, sprintf('chsSorted_rep%d', rr)), 'svg');
            
            close all;
    
            % Accuracy values. positive means responsive/excluded from CAR
            % We keep these variables as named here, but note that FN and FP are now renamed RCM and NCM in the manuscript.
            TP = sum(find(respBool) > nCAR); % responsive channels successfully excluded from CAR (above the cutoff)
            TN = sum(find(~respBool) <= nCAR); % NR channels successfully below or at cutoff
            FN = sum(find(respBool) <= nCAR); % responsive channels incorrectly included in CAR. *This matters most
            FP = sum(find(~respBool) > nCAR); % NR channels incorrectly excluded from CAR
            
            % same for the global threshold
            TPglob = sum(find(respBool) > nCARglob);
            TNglob = sum(find(~respBool) <= nCARglob);
            FNglob = sum(find(respBool) <= nCARglob);
            FPglob = sum(find(~respBool) > nCARglob);
    
            fid = fopen(fullfile(outdir, sprintf('accuracy_rep%d.txt', rr)), 'w');
            fprintf(fid, 'TP\t%d\nTN\t%d\nFN\t%d\nFP\t%d\n', TP, TN, FN, FP);
            fprintf(fid, 'TPglob\t%d\nTNglob\t%d\nFNglob\t%d\nFPglob\t%d', TPglob, TNglob, FNglob, FPglob);
            fclose(fid);
    
            % save snr info
            save(fullfile(outdir, sprintf('snr_rep%d.mat', rr)), 'snr', 'Psig', 'PnoiseRand');
    
        end
        
        fprintf('\n');
        
    end
end

%return

%% Save a few examples of how the channels look with the same signal but variable SNR

outdir = fullfile('output', 'simLoopSNR');

rng(25); % a representative seed with the 3 channels showing comparable SNR to the population geomeans

nchsEx = 4;
chsResp = 1:3;

% i) artifact only
V0 = zeros(length(tt), nchsEx);
Aart = 50 + rand(nchsEx, 1)*5; % slightly different artifact amplitudes for each channel
artifact = sin(2*pi*600*tt)';
artifact(tt < 0 | tt > 0.002) = 0;
V0 = V0 + artifact*Aart';

% A) Add the evoked potentials
A = 100;
V1 = V0;
sig = genRandSig(tt, length(chsResp), A);
V1(:, chsResp) = V0(:, chsResp) + sig;

% no global noise, copy V1 over
V2 = V1;

% keep the line noise the same for all
phLN = rand(ntrs, 3)*2*pi; % LN phases
LN = zeros(length(tt), ntrs);
for ii = 1:ntrs
    LN(:, ii) = 8*sin(2*pi*60*tt - phLN(ii, 1)) + 2*sin(2*pi*120*tt - phLN(ii, 2)) + 1*sin(2*pi*180*tt - phLN(ii, 3));
end

% create the same series of random noise for all noise levels, just varying amplitude
BNBase = randn(2*length(tt), ntrs);
noiseRandBase = randn(nchsEx, 2*length(tt), ntrs);

for brownCoef = 0.4:0.1:1.0

    % C) Add common noise to all channels at each trial
    V3 = repmat(V2', 1, 1, ntrs); % ch x time points x trial
    
    % brown noise shared across channels
    BN = cumsum(brownCoef*BNBase); % variable brown noise coefficient now
    BN = ieeg_highpass(BN, srate, true);
    BN = BN((0.5*length(tt)+1) : 1.5*length(tt), :);
    
    noiseCommon = LN + BN;
    V3 = V3 + shiftdim(noiseCommon, -1);
    
    % D) add random brown noise
    noiseRand = cumsum(brownCoef*noiseRandBase, 2); % give double the number of time points so we can highpass it
    for ii = 1:nchsEx
        noiseRand(ii, :, :) = ieeg_highpass(squeeze(noiseRand(ii, :, :)), srate, true);
    end
    noiseRand = noiseRand(:, (0.5*length(tt)+1) : 1.5*length(tt), :);
    V4 = V3 + noiseRand;

    % Calculate the SNRs in the example
    Psig = sum(sig(tt >= snrWin(1) & tt < snrWin(2), :).^2)' / diff(snrWin); % express as per second
    noiseSum = noiseRand + shiftdim(BN, -1);
    PnoiseRand = squeeze(sum(noiseSum(chsResp, tt >= snrWin(1) & tt < snrWin(2), :).^2, 2) / diff(snrWin));
    snr = repmat(Psig, 1, ntrs) ./ PnoiseRand;
    snrByCh = geomean(snr, 2); % average across trials for each channel
    fprintf('%0.2f, ', snrByCh(:)); fprintf('\n');

    % Plot and save
    yspace = 200;
    figure('Position', [200, 200, 600, 600]);
    ys = ieeg_plotTrials(tt, mean(V4, 3)', yspace, [], 'k', 'LineWidth', 1);
    hold on
    for jj = 1:4 % plot the individual trials
        plot(tt, ys(jj) + squeeze(V4(jj, :, :)), 'Color', [0.5, 0.5, 0.5]);
    end
    hold off
    ylim([ys(end)-yspace, yspace]);
    kids = get(gca, 'Children');
    set(gca, 'Children', [kids(end-1:-2:end-7); kids(1:end-8); kids(end:-2:end-6)]); % first (top) the means, then the trials, lastly the horzline
    xlabel('Time from Stim. (s)');
    saveas(gcf, fullfile(outdir, sprintf('exampleSNRs_noise%0.1f_nchs%d-%d.png', brownCoef, nchsEx, length(chsResp))));
    saveas(gcf, fullfile(outdir, sprintf('exampleSNRs_noise%0.1f_nchs%d-%d.svg', brownCoef, nchsEx, length(chsResp))));
    close(gcf);
end
