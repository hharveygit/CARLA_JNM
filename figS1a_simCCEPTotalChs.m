%% This script is derived from simCCEPLoop, except an outer loop varies the number of total channels
% We simulate 30 reps for each combination of total channels and responsiveness and test CARLA's accuracy
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

ntrs = 12; % number of trials
reps = 30; % number of reps for each condition

cmSens = [1, 165/255, 0]; % use orange as color for the sensitive optimum

Aglobal = 0;

% Total number of channels to simulate for each condition
nChsToTest = [25, 20, 15, 10];


%% Outer loop through number of total channels

rng(nChsToTest(1)); % Set seed depending on what number of channels we want, for reproducibility

for nn = 1:length(nChsToTest)

    nchs = nChsToTest(nn);
    fprintf('Number of channels = %d\n', nchs);

    %% Inner loop through number of responsive channels
    for nresp = 0:0.2*nchs:0.8*nchs % 0 to 80% responsiveness, by 20% intervals
        
        %% Create data
        
        fprintf('%d of %d channels responsive\n', nresp, nchs);
        
        outdir = fullfile('output', 'simLoopNChs', sprintf('nchs%d-%d', nchs, nresp));
        mkdir(outdir);
    
        chsResp = 1:nresp; % first x channels responsive, for ease. order doesn't matter. will sort and color them red
    
        %%
        for rr = 1:reps % multiple repetitions at each significance
            fprintf('.');
            
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
            BN = cumsum(0.4*randn(2*length(tt), ntrs)); % variable brown noise coefficient now
            BN = ieeg_highpass(BN, srate, true);
            BN = BN((0.5*length(tt)+1) : 1.5*length(tt), :);
    
            noiseCommon = LN + BN;
            V3 = V3 + shiftdim(noiseCommon, -1);
    
            % D) add random brown noise
            noiseRand = cumsum(0.4*randn(nchs, 2*length(tt), ntrs), 2); % give double the number of time points so we can highpass it
            for ii = 1:nchs
                noiseRand(ii, :, :) = ieeg_highpass(squeeze(noiseRand(ii, :, :)), srate, true);
            end
            noiseRand = noiseRand(:, (0.5*length(tt)+1) : 1.5*length(tt), :);
            V4 = V3 + noiseRand;
    
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
    
        end
        
        fprintf('\n');
        
    end
end
