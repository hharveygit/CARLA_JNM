%% This script calculates and prints the number of channels used in the analysis for each subject ("good" in at least 1 run of data)
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
% path to subject data
dataPath = 'data';

% subjects and runs corresponding to each subject to load
subs = {'sub-1', 'sub-2', 'sub-3', 'sub-4'};
runs = {[1, 2, 3, 5], [12, 13, 14], [2, 3, 4], [1, 2]};

% iterate through each subject, then each run, and calculate number of good channels
for ss = 1:4

    sub = subs{ss};
   
    seegChStatus = {};
    
    % iterate across runs for each subject
    for rr = 1:length(runs{ss})
        
        run = runs{ss}(rr);
        
        channelsPath = fullfile(dataPath, sub, 'ses-ieeg01', 'ieeg', sprintf('%s_ses-ieeg01_task-ccep_run-%02d_channels.tsv', sub, run));
        channels = ieeg_readtableRmHyphens(channelsPath);
        seegChStatus = [seegChStatus, channels.status(strcmp(channels.type, 'SEEG'))];
    
    end
    
    % Good channels are defined as any that are not bad across all runs (and therefore not used in any analysis)
    goodChs = any(strcmp(seegChStatus, 'good'), 2);
    fprintf('Number of good sEEG channels across all runs in %s = %d\n', sub, sum(goodChs));

end
