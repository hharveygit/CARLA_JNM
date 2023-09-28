%% File to simulate an evoked potential as part of the code to test CARLA
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
function sig = genRandSig(tt, nChs, amp)
    % tt = time points
    % amp = average amplitude of s1
    % nChs = number of channels to generate
    % Adds s1 and s2. s1 represents faster, direct response. s2 represents slower, indirect response
    
    % signal parameters
    A = 0.8*amp + 0.4*amp*rand(nChs, 1); % amplitudes
    
    tau1 = 0.01 + rand(nChs, 1)*0.02; % time constant to dampen the faster sinusoid (s1)
    tau2 = 0.06 + rand(nChs, 1)*0.08;
    tauOn1 = 0.005; % time constant of subtracted exponential decay for both
    tauOn2 = 0.025; % time constant of subtracted exponential decay for both
    
    f1 = 8 + rand(nChs, 1)*4; % freq (hz) of first sinusoid (mean period = 100 ms)
    f2 = 1 + rand(nChs, 1)*2;  % freq of second sinusoid, (mean period = 500 ms)
    
    ph1 = rand(nChs, 1)*2*pi; % phase of first sinusoid
    ph2 = rand(nChs, 1)*2*pi; % phase of second sinusoid
    
    sig = zeros(length(tt), nChs);
    for ii = 1:nChs
            sig(:, ii) = A(ii)*( (exp(-tt/tau1(ii)) - exp(-tt/tauOn1)).*sin(2*pi*f1(ii)*tt - ph1(ii)) + (exp(-tt/tau2(ii)) - exp(-tt/tauOn2)).*sin(2*pi*f2(ii)*tt - ph2(ii)) );
    end
    sig(tt < 0, :) = 0; % make causal
    
end