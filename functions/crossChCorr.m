function rsqs = crossChCorr(tt, V, srate)
%
%   This function calculates how mean cross-channel R-squared change at each CAR size
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
    % Ensure ch x time points structure if k = 1
    if ismatrix(V)
        if size(V, 1) == length(tt), V = V'; end
    end
    
    % Estimate sampling frequency from time points if necessary (not recommended)
    if nargin < 3 || isempty(srate)
        srate = (length(tt) - 1) / (max(tt) - min(tt));
        warning('sampling frequency was estimated to be %0.02f', srate);
    end
    
    nChs = size(V, 1);
    nTrs = size(V, 3);
    assert(size(V, 2) == length(tt), 'Error: second dimension does not match tt. Data should be channels x timepoints x trials');
    
    stats = struct();
    
    % Create 60, 120, 180 Hz notch filters
    for ff = 60:60:180
        dNotch(ff/60) = designfilt('bandstopiir', 'FilterOrder', 4, ...
                                'DesignMethod', 'butter', ...
                                'HalfPowerFrequency1', ff-2, ... 
                                'HalfPowerFrequency2', ff+2, ...
                                'SampleRate', srate);
    end
    
    % notch filter signals to use for covariance calculating
    Vclean = V;
    for ff = 60:60:180
        for ii = 1:nTrs % filter per trial
            Vclean(:, :, ii) = filtfilt(dNotch(ff/60), Vclean(:, :, ii)')';
        end
    end
    
    % extract data on response interval to perform analysis on
    twin = [0.01, 0.3]; % should this be an input parameter? For now keep hardcoded
    Vseg = Vclean(:, tt >= twin(1) & tt <= twin(2), :);
    VorigSeg = V(:, tt >= twin(1) & tt <= twin(2), :);
    
    % Rank channels. If single trial rank based on increasing variance. If more than 1 trial, rank by covariance (self-self excluded). Better estimator of consistent structure
    stats.vars = nan(nChs, 1);
    for ii = 1:nChs
        covCurr = cov(squeeze(Vseg(ii, :, :)));
        stats.vars(ii) = mean(covCurr(logical(triu(ones(size(covCurr)), 1))), 'all');
    end
    [~, stats.order] = sort(stats.vars);
    
    %
    rsqs = nan(nChs, nChs, nChs+1); % first column corresponds to no rereferencing
        
    % Pull subset U (order(1:ii)) from V, growing U each time by 1, to construct CAR
    for ii = 0:nChs % min CAR size is 2
        
        % rereferenced signal is calculated from the non-notch data
        if ii == 0
            VsegReref = VorigSeg;
        else
            VsegReref = VorigSeg - mean(VorigSeg(stats.order(1:ii), :, :), 1);
        end

        % apply notch filter to the rereferenced signal before calculating correlations, to reduce effect of line noise
        for ff = 60:60:180
            for jj = 1:nTrs % filter per trial
                VsegReref(:, :, jj) = filtfilt(dNotch(ff/60), VsegReref(:, :, jj)')';
            end
        end

        VsegRerefMean = mean(VsegReref, 3); % mean across trials
        
        r = corr(VsegRerefMean');
        r(1:(nChs+1):end) = nan; % omit self-self correlations

        rsqs(:, :, ii+1) = r.^2; % convert to R^2 value
        
        
    end
    
end
