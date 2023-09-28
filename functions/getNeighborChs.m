%% Returns the neighbors within n contacts for input channels
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
function neighbors = getNeighborChs(chs, n, selfincl)
    if nargin < 3, selfincl = true; end % whether or not the input channel(s) should be included as neighbors
    if nargin < 2 || isempty(n), n = 1; end
    
    isDigit = @(x) x > 47 & x < 58; % returns true for char array elements that are digits (0 - 9)
    
    if ischar(chs), chs = {chs}; end % if only 1 channel is requested
    
    neighbors = {};
    for ii = 1:length(chs)
        ch = strip(upper(chs{ii}));
       
        % split channel name into lead and contact
        lead = ch(~isDigit(ch));
        contact = str2double(ch(isDigit(ch)));
        
        % all neighoring contact integers within n
        neighborcontacts = contact + (-n:n);
        neighborcontacts(neighborcontacts < 1 | neighborcontacts > 18) = [];
                
        neighbors = [neighbors, arrayfun(@(x) sprintf('%s%d', lead, x), neighborcontacts, 'UniformOutput', false)];
    end
    
    % remove redundant entries
    neighbors = unique(neighbors, 'stable');
    
    % remove input channels if necessary
    if ~selfincl
        neighbors = setdiff(neighbors, chs);
    end
    
end