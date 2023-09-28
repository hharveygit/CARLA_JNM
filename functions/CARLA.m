function [Vout, CAR, stats] = CARLA(tt, V, srate, sens, nboot)
%
%   This function performs Common Average Re-referencing by Least Anticorrelation (CARLA).
%   Channels are ranked in order of increasing covariance across trials (or variance if single trial). 
%   Next, channels are iteratively added into a CAR, and the mean correlation between each unre-referenced channel and all re-referenced channels is 
%   calculated, for each bootstrapped mean signal. zmin denotes the mean correlation belonging, at each bootstrap, to the most anticorrelated channel.
%   The optimum CAR size corresponds to when there is least anticorrelation (when zmin takes on its least negative value)
%
%   Vout = CARLA(tt, V); [NOT RECOMMENDED]
%   Vout = CARLA(tt, V, srate);
%   [Vout, stats] = CARLA(tt, V, srate, nboot);
%       tt =        1 x t num. Time points matching V, in seconds
%       V =         n x t x k num, or n x t num. Signal data to be re-referenced, with n channels by t time points by k trials (2 dimensional if k=1).
%       srate =     num. Sampling frequency of data. If not given, this will be estimated from the time points (not recommended because of potential imprecision).
%       sens =      (optional) boolean. Determines whether the more sensitive cutoff will be used, corresponding to before the first significant decrease in zMinMean.
%                       Ignored if V is n x t (needs trials to calculate significance). Default = false (returns global maximum in zMinMean)
%       nboot =     (optional) integer. Number of bootstrapped mean signal samples to generate when calculating correlations. Default = 100. Ignored if data is n x t
%
%   RETURNS:
%       Vout =      n x t x k num, or n x t num. Re-referenced signal data.
%       CAR =       t x k num, or 1 x t num. The CAR at each trial
%       stats =     struct, containing fields:
%                       chsUsed =       m x 1 num, m < n. The sorted n channel numbers used in the optimum CAR
%                       vars =          n x 1 num. Trial covariance of each channel (variance if K = 1)
%                       order =         n x 1 num. The channels sorted in order of increasing covariance (variance if K = 1)
%                       zMin =          n x n x nboot num, or n x n num if k = 1. The z (Fisher z-transformed Pearson's r) for the channel with the most negative 
%                                           average z, at each bootstrap. Dimensions are (z to each channel) x (number of channels in CAR subset) x (bootstrap)
%                       zMinMean =      n x nboot num, or 1 x n num if k = 1. Mean zMin across the first dimension (averaged across all "target" rereferenced channels)
%                                           for each bootstrap if V was 3 dimensional
%
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
%   2023/04/27 by Harvey Huang                       
%    
    % Ensure ch x time points structure if k = 1
    if ismatrix(V)
        if size(V, 1) == length(tt), V = V'; end
    end
    
    if nargin < 5 || isempty(nboot), nboot = 100; end
    if nargin < 4 || isempty(sens), sens = false; end
    
    % Estimate sampling frequency from time points if necessary (not recommended)
    if nargin < 3 || isempty(srate)
        srate = (length(tt) - 1) / (max(tt) - min(tt));
        warning('sampling frequency was estimated to be %0.02f', srate);
    end
    
    nChs = size(V, 1);
    nTrs = size(V, 3);
    assert(size(V, 2) == length(tt), 'Error: second dimension does not match tt. Data should be channels x timepoints x trials');
    
    stats = struct();
    
    % Notch filter signals to remove line noise. This is not used in the final output
    Vclean = V;
    for ff = 60:60:180
        dNotch = designfilt('bandstopiir', 'FilterOrder', 4, ...
                                'DesignMethod', 'butter', ...
                                'HalfPowerFrequency1', ff-2, ... 
                                'HalfPowerFrequency2', ff+2, ...
                                'SampleRate', srate);
        for ii = 1:nTrs % filter per trial
            Vclean(:, :, ii) = filtfilt(dNotch, Vclean(:, :, ii)')';
        end
    end
    
    % extract data on response interval to perform analysis on
    twin = [0.01, 0.3]; % should this be an input parameter? For now keep hardcoded
    Vseg = Vclean(:, tt >= twin(1) & tt <= twin(2), :);
    
    % Rank channels. If single trial rank based on increasing variance. If more than 1 trial, rank by covariance (self-self excluded). Better estimator of consistent structure
    if nTrs == 1
        stats.vars = var(Vseg, 0, 2);
    else
        stats.vars = nan(nChs, 1);
        for ii = 1:nChs
            covCurr = cov(squeeze(Vseg(ii, :, :)));
            stats.vars(ii) = mean(covCurr(logical(triu(ones(size(covCurr)), 1))), 'all');
        end
    end
    [~, stats.order] = sort(stats.vars);

    % channels x CARsize x bootstrapped sample
    if nTrs == 1
        stats.zMin = nan(nChs, nChs); % no bootstrapping if only single trial
    else
        stats.zMin = nan(nChs, nChs, nboot);
    end
        
    % Pull subset U (order(1:ii)) from V, growing U each time by 1, to construct CAR
    for ii = 2:nChs % min CAR size is 2
        
        % Vseg after subtracting the putative CAR on subset U
        VsegReref = Vseg - mean(Vseg(stats.order(1:ii), :, :), 1);
        
        if nTrs == 1 % no bootstrapping, test correlation at the single trial level
            Useg = Vseg(stats.order(1:ii), :)'; % UNREREFERENCED channels in subset to correlate with the rereferenced signals
            UsegReref = VsegReref(stats.order(1:ii), :)'; % rereferenced signals in subset
            
            r = corr(Useg, UsegReref); % rows = which unrereferenced channel was used.
            r(1:(ii+1):end) = nan; % omit self-self correlations
            z = atanh(r); % fisher z-transform
            
            % choose "most responsive" channel with greatest anticorrelation, on average
            [~, kkMost] = min(mean(z, 2, 'omitnan'));
            stats.zMin(stats.order(1:ii), ii) = z(kkMost, :);
            
            continue;
            
        end
        
        % if nTrs > 1
        for jj = 1:nboot
            
            % bootstrap to calculate mean signal
            inds = datasample(1:nTrs, nTrs);
            Useg = mean(Vseg(stats.order(1:ii), :, inds), 3)';
            UsegReref = mean(VsegReref(stats.order(1:ii), :, inds), 3)';
            
            r = corr(Useg, UsegReref); % rows = which unrereferenced signal was used.
            r(1:(ii+1):end) = nan; % omit self-self correlations
            z = atanh(r); % fisher z-transform
            
            % choose "most responsive" channel with greatest anticorrelation, on average, across channels, for this bootstrapped sample
            [~, kkMost] = min(mean(z, 2, 'omitnan'));
            stats.zMin(stats.order(1:ii), ii, jj) = z(kkMost, :);
            
        end
        
    end
    
    % calculate mean across the unreferenced target channels
    stats.zMinMean = mean(stats.zMin, 1, 'omitnan');
    
    if sens && nTrs > 1 % find sensitive optimum
        
        nMin = ceil(0.1*nChs); % minimum number of channels to start at, for stability. Hardcoded at 10 %
        zMMxTrs = mean(stats.zMinMean(1, :, :), 3); % mean ZMinMean across bootstrapped samples for each n
        
        ii = nMin;
        while ii <= nChs
            
            if ii == nChs % at the end. no more comparisons needed
                nOptimum = ii;
                break
            end
            
            % loop until next local maximum found
            if zMMxTrs(ii + 1) > zMMxTrs(ii)
                ii = ii + 1;
                continue;
            end
            
            zMMxTrs(1:ii-1) = nan; % set all earlier values to nan, to avoid them being found as nextGreater or as trough
            
            % find index of the next greater zMinmean (ignores earlier nans), if it exists
            nextGreater = find(zMMxTrs > zMMxTrs(ii), 1, 'First');
            if isempty(nextGreater) % no more samples greater than current peak. Optimum is identical to global optimum. Done
                nOptimum = ii;
                break
            end
            zMMCurr = squeeze(stats.zMinMean(1, ii, :)); % all bootstrapped samples at current max
            [~, idx] = min(zMMxTrs(1:nextGreater-1)); % locate Zminmean trough before the next max (ignores earlier nans)
            zMMTrough = squeeze(stats.zMinMean(1, idx, :)); % all bootstrapped samples at trough before next greater value
            
            % all pairwise differences in the mean correlation estimate between the trough and current peak.
            % Note that zMinMean is the estimated correlation of the mean signal. We are testing the difference of zMinMean directly,
            % Not the mean of zMinMean, which would be redundant and inappropriate because the bootstrapped samples are not independent.
            [X, Y] = meshgrid(zMMTrough, zMMCurr);
            diffs = X - Y; diffs = diffs(:);
            conf = [-inf, prctile(diffs, 95)]; % left-tailed 95% confidence interval (H_a = difference is less than 0)
            
            % found significant decrease. Current max is optimal.
            if conf(2) < 0
                nOptimum = ii;
                break
            end
            
            % Difference was not significant. Start at next greater sample
            ii = nextGreater;
            
        end
        
        if nOptimum == nMin, warning('Optimum CAR was detected at minimal 10% floor'); end
        
    else
        % Find optimum n at max (least negative) zminMean, averaged across trials.
        [~, nOptimum] = max(mean(stats.zMinMean, 3));
    end
    
    % delete first single dimension if 3D, after finding nOptimum
    stats.zMinMean = squeeze(stats.zMinMean);
    
    % which channels were used to make CAR
    stats.chsUsed = sort(stats.order(1:nOptimum));
    
    CAR = mean(V(stats.chsUsed, :, :));
    
    % Using the original (unfiltered signals)
    Vout = V - CAR;
    
    CAR = squeeze(CAR); % squeeze out the first dimension to make t x k, if 3 dimensional
    
end
