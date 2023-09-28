function dataOut = rmBadTrialsAnnots(dataIn, chs, annots)
%% Removes bad trials according to text in a cell array of chars, annotations
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
    assert(size(dataIn, 1) == length(chs), 'Error: First dimension of input data does not match channel names');
    assert(size(dataIn, 3) == length(annots), 'Error: Third dimension of input data does not match annotations length');
    
    dataOut = dataIn;
    
    annots = upper(strip(annots)); % clean up, make uppercase
    
    % Remove trials that are bad across all channels
    trs2Rm = strcmp(annots, 'ALL');
    fprintf('Removing %d/%d trials\n', sum(trs2Rm), length(annots));
    dataOut(:, :, trs2Rm) = [];
    annots(trs2Rm) = [];
    
    % find channels bad in individual trials and set as nan, to be removed later
    for ii = 1:length(annots)
        try
            annot = annots{ii};
            
            % Channels mentioned in description. matches strs starting with L/R, followed by 1 or more A-Z letters, then 1 or more numerical digits
            out = regexp(annot, '[LR][A-Z]+\d+', 'match');
            if isempty(out), continue; end

            % identify which channels are bad at this trial
            idxbad = ismember(chs, out);
            fprintf('Setting %d bad channels at trial number %d to NaN values: %s\n', sum(idxbad), ii, string(join(out, ', ')));
            dataOut(idxbad, :, ii) = nan;
        catch
        end
    end
    
end