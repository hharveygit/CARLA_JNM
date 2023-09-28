%% Older, alternative version of generating random signals. Used for schematic purposes, not for quantification
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
function sig = genRandSigOrig(tt, nChs, amp)
    % tt = time points
    % amp = average amplitude of s1
    % nChs = number of channels to generate
    % Adds s1 and s2. s1 represents faster, direct response. s2 represents slower, indirect response
    
    % signal parameters
    A = 0.8*amp + 0.4*amp*rand(nChs, 1); % amplitudes
    
    tau1 = 0.05 + rand(nChs, 1)*0.06; % time constant of slower exp decay
    tau2 = ones(nChs, 1)*0.01; % time constant of second exp decay (shorter, fixed)
    tau3 = 0.025; % time constant to dampen the nested, faster sinusoid (s1), corresponds to half the half-period of
    
    f1 = 8 + rand(nChs, 1)*4; % freq (hz) of first sinusoid (mean period = 100 ms)
    f2 = 2 + rand(nChs, 1)*2; % freq of second sinusoid, (mean period = 250 ms)
    
    ph1 = rand(nChs, 1)*2*pi; % phase of first sinusoid
    ph2 = rand(nChs, 1)*2*pi; % phase of second sinusoid
    
    sig = zeros(length(tt), nChs);
    for ii = 1:nChs
            sig(:, ii) = A(ii)*( (exp(-tt/tau1(ii)) - exp(-tt/tau2(ii))) .* (exp(-tt/tau3).*sin(2*pi*f1(ii)*tt - ph1(ii)) + sin(2*pi*f2(ii)*tt - ph2(ii))));
    end
    sig(tt < 0, :) = 0; % make causal
    
end