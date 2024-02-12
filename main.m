%% Use each section of this main script to output all results and figures for the manuscript.
%
%   2023/09/27
%
%    If this code is used in a publication, please cite the manuscript:
%    "CARLA: Adjusted common average referencing for cortico-cortical evoked potential data"
%    by H Huang, G Ojeda Valencia, NM Gregg, GM Osman, MN Montoya,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    CARLA manuscript package: main script
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Unless specified, each section of code below can be run directly from this script to generate all outputs.
%   The scripts listed may also be opened so that individual sections within can be separately examined and run.
%   Note that in output filenames, "zminmean" refers to zeta in the manuscript
%
%   Subject CCEP and anatomical data will be available in BIDS format upon manuscript publication.
%   In the meantime, all results on simulated data can be reproduced with the appropriate section below.
%
%% Configure paths and dependencies

% Make sure you are currently working in the root manuscript directory (CARLA)
addpath('functions');

% Data should be downloaded from OpenNeuro and copied to a folder here (./data)
dataPath = 'data';

% Add paths to other dependencies. Change the dummy paths listed here to where they are located on your system.
addpath(genpath('path/to/mnl_ieegBasics'));
addpath(genpath('path/to/matmef'));
addpath('path/to/mnl_seegview');
addpath(genpath('path/to/vistasoft'));
addpath('path/to/spm12'); addpath(fullfile('path/to/spm12', 'config'));
addpath(genpath('path/to/freesurfer/7.1.1/matlab'));
addpath(genpath('path/to/gifti'));

% To save figures as vectorized objects
set(0, 'DefaultFigureRenderer', 'painters');

%% i) Generate T1 MRI slices with stimulation site electrodes overlaid (Figures 4A, 6C)

% This was done using sliceGUI.m in the mnl_seegview repository:
% https://github.com/MultimodalNeuroimagingLab/mnl_seegview

% Define paths to defaced T1 MRI and matching electrodes table
niiPath = fullfile(dataPath, 'derivatives', 'freesurfer', 'sub-1', 'sub-1_ses-mri01_T1w_acpc_deFaced.nii');
electrodes = fullfile(dataPath, 'sub-1', 'ses-ieeg01', 'ieeg', 'sub-1_ses-ieeg01_electrodes.tsv');

% Run sliceGUI, and manually interact with it
sliceGUI(niiPath, electrodes);

%% Table 1: Calculate the number of measurement electrodes (channels) used in each subject

% outputs are directly printed to command window
tab1_dispNumElectrodes

%% Figure 1: Save simulated CCEP channels and their theoretical components

% outputs saved to .\output\simComponents
fig1_simCCEPComponents

%% Figure 2: Methods illustration of how CARLA works, on single-trial simulated CCEP data

% outputs saved to .\output\simulationTrial
fig2_CARLAonSingleTrial

%% Figure 3: Save example of enveloped sinusoids used to simulate evoked potentials

% outputs saved to .\output\EPConstruct
fig3_EPConstruct

%% Figure 4: Simulate 50 channels with responsiveness varying from 0 to 50, with 30 repetitions each, and apply CARLA

% First, run this to generate outputs for each number of responsive channels
% outputs saved to .\output\simLoop
Aglobal = 0; % no global signal desired (0 amplitude)
simCCEPLoop;

% Then, keeping Aglobal in the workspace, run this to make summary figures (4C, D), as well as identify median candidates for 4A, B
% .\output\simLoop will be renamed .\output\simLoopNoGlob, where outputs are saved.
summarizeSimCCEPLoop

%% Figure S1: Simulate datasets like Fig 4, but with fewer total channels (25, 20, 15, 10) and calculate CARLA accuracy

% First, run this to generate outputs for each total channel size and fraction of responsive channels
% outputs are saved to .\output\simLoopNChs
figS1a_simCCEPTotalChs

% Then run this to make summary figure (S1A). S1B was constructed with examples most similar to the median of this summary output
figS1b_summarizeSimCCEPTotalChs

%% Figure S2: Simulate datasets with lower signal-to-noise ratio (SNR) and calculate CARLA accuracy

% Brown noise coefficients to save outputs for. Split up into chunks for more efficient parallelization
% SNR is varied by changing the Brown Noise coefficient. The population mean SNR is estimated at the end with figS2b_summarizeSimCCEPSNRs, below.
brownCoefRanges = [0.4, 0.5;
                   0.6, 0.7;
                   0.8, 0.9;
                   1.0, 1.0];

% First, run this to generate outputs for each SNR level and number of responsive channels. This also saves the plots for Figure S2A
% outputs are saved to .\output\simLoopSNR
% this loop can be parallelized if desired
for ii = 1:4
    brownCoefRange = brownCoefRanges(ii, :);
    figS2a_simCCEPSNRs
end

% Then run this to make summary figure, S2B.
% This also calculates the geomean (sigma) SNR for each Brown Noise coefficient condition (variable snrGeomean in workspace)
figS2b_summarizeSimCCEPSNRs

%% Figure 5: Apply CARLA to real data (sub 1, stim site 1)

% Configure the subject, data run, and stim site
sub = 'sub-1';
run = 2;
site = 'RMO8-RMO9';

% which channels to plot in panel F
chs2Plot = [100, 203];

% outputs saved to .\output\realCCEPs\. These correspond to figures 5B-F. To generate 5A (stim site overlaid on T1 MRI), see section i), above.
applyCARLARealCCEPs

%% Figure S3: Average channel responses for sub 1, stim site 1 after applying CAR vs. CARLA

% This supplemental figure is a direct addendum to Figure 5, showing the mean signal for the first 50 channels after rereferencing by CARLA or CAR
% Run after the previous section, keeping all workspace variables
figS3_compareCARvCARLA

%% Figure S4: Apply CARLA to real data at another 4 stim sites in subjects 1-4

% config parameters for each subject shown in S1
config = struct('sub', {'sub-1', 'sub-2', 'sub-3', 'sub-4'}, ...
                'run', {1, 12, 3, 2}, ...
                'site', {'ROC11-ROC12', 'RZ4-RZ5', 'LA2-LA3', 'RO3-RO4'});
chs2Plot = 1; % dummy value, no individual channels are plotted in Figure S1

for ii = 1:length(config)
    sub = config(ii).sub;
    run = config(ii).run;
    site = config(ii).site;
    
    % outputs saved to .\output\realCCEPs\.
    % Omission of the middle section of sorted channels (subjects 1-3) in Figure S1 was done manually in illustrator.
    applyCARLARealCCEPs
end

%% Figure S5: Save simulated CCEP channels, this time with a global signal, and their theoretical components

% outputs saved to .\output\simComponentsGlobal
figS5_simCCEPComponentsGlobalSig

%% Figure 6A, B Simulate 50 channels with responsiveness from 0 to 50, like Figure 4, but now with a global signal added

% First run this to generate outputs for each number of responsive channels
% outputs saved to .\output\simLoop
Aglobal = 25; % global signal with mean amplitude 25
simCCEPLoop;

% Then, keeping Aglobal in the workspace, run this to make summary figure (6B), and to identify a trial candidate for 6A
% .\output\simLoop will be renamed .\output\simLoopWithGlob, where outputs are saved.
summarizeSimCCEPLoop

% To generate figure 6C, see section i) on using sliceGUI, above.

%% Figure 6D-F: Apply CARLA to real data containing global signal (sub 1, stim site 2)

% Configure the subject, data run, and stim site
sub = 'sub-1';
run = 3;
site = 'RK5-RK6';

% which channels to plot in panel F
chs2Plot = [23, 51, 205];

% outputs saved to .\output\realCCEPs\
applyCARLARealCCEPs

%% Figure 7: Loop through stimulation sites in each subject and compute optimal CARLA cutoff
% This also calculates the cross-channel R-squared values for all possible adjusted common average sizes, to be summarized in Figure 8

% First run this script. Outputs are saved to .\output\realCCEPsLoop\.
applyCARLARealCCEPsLoop

% Then run this script to generate outputs for Figure 7.
% Note that this step uses Freesurfer outputs (in data/derivatives/freesurfer/.) to localize each stim site to a tissue type
fig7_summarizeCARLARealCCEPs

%% Figure 8: Calculate cross-channel R-squared values for no CAR vs different types of fixed CARs vs CARLA

% outputs saved to .\output\realCCEPsLoop\.
fig8_summarizeInterChCorr

