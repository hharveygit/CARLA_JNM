%% Creates a 2-column cell array where the first column indicates unique occurreences in the input, and second col tracks indices of each occurrence
% Modeled after Python groupby function
%
%   2024/02/14
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
function groups = groupby(strs)
    
    gs = cellfun(@(x) find(strcmp(x, strs)), unique(strs), 'UniformOutput', false);
    groups = [unique(strs), gs];
    
end