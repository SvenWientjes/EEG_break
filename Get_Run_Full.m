%% Simulate EEG
% This script simulates one 'participant' performing some kind of cognitive
% experiment. It switches strategy (superimposed ERPs) halfway through the
% experiment.
% -----------------------------------------------------------------------------------------
% This script relies on the SEREEGA toolbox for Simulating Event Related
% EEG Activity. This toolbox can be downloaded from:
% https://github.com/lrkrol/SEREEGA
% and needs to be added to the Matlab Path.
%
% SEREEGA relies on functions from the EEGLAB toolbox. It can be downloaded
% from:
% https://sccn.ucsd.edu/eeglab/download.php
% The folder called 'popfunc' needs to be added to the Matlab Path for
% SEREEGA to function in this script. It is best to NOT add this folder
% permanently to the path, so addpath is used. Alter the path to popfunc on
% your computer for this script to work.
%
% This script relies on Nick Yeung's EEG noise simulator. It can be
% downloaded from:
% https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
% and needs to be added to the Matlab Path.
% -----------------------------------------------------------------------------------------
% This script has several alterable variables:
% Hz      = Simulated EEG sampling rate
% N       = Number of trials per strategy
% Lms     = Length of trials in milliseconds
% nLoc    = Number of ERP sources per strategy
% nChange = Number of ERP sources that change on strategy shift
% SNR     = Signal to Noise Ratio of ERP vs EEG noise
%
% dat_i   = Determines how many participants we will simulate. They
% will be saved in separate folders with the names 'Run_[dat_i]'.
%
% under '%Save relevant datas', you can change the path where the simulated
% participants will be saved. To also save the EEGLAB data format, also
% change the path in the pop_saveset function.
% -----------------------------------------------------------------------------------------
% Written by Sven Wientjes at the UvA
% Master Brain & Cognitive Sciences research intern at the Van Maanen lab for Mathematical Psychology and
% Cognitive Modeling (MPCM)

%% Define preliminary variables
Hz      = 500;
N       = 10;
Lms     = 3500;
nLoc    = 7;
nChange = 6;
SNR     = 1;

addpath('C:\Program Files\MATLAB\R2018b\toolbox\eeglab14_1_2b\functions\popfunc')

%% Start the loop for generating data
for dat_i = 1:30 
    [EEG_raw, info] = simulate_eeg(Hz,  N, Lms, nLoc, nChange, SNR);

    % Save relevant datas
    mkdir(['D:\StageData\SNR-1_Tr-10_L-7_S-abo\Run_', int2str(dat_i)])
    save(['D:\StageData\SNR-1_Tr-10_L-7_S-abo\Run_', int2str(dat_i), '\Data'])
    
    % Prepare for eeglab data: AMICA requires
    EEG                   = struct()
    EEG.icawinv           = 0;
    EEG.icaweights        = zeros(64);
    EEG.icasphere         = zeros(64);
    EEG.trials            = N*2;
    EEG.xmax              = 3500;
    EEG.xmin              = 0;
    EEG.icaact            = [];
    EEG.chanlocs          = struct();
    EEG.filename          = 'placeholder'
    EEG.option_computeica = 0;
    
    % Check file correctness and save as .set / .fdt combo
    EEG = pop_editset(EEG, 'setname', 'TheData', 'data', EEG_raw, 'dataformat', 'array', 'srate', 500, 'nbchan', 64, 'pnts', Lms/1000*Hz*N*2)
    pop_saveset(EEG, 'filename', ['Run_', int2str(dat_i)], 'filepath', [pwd, '\SNR-1_Tr-10_L-7_S-abo\Run_', int2str(dat_i)])
end

eeg_checkset(EEG)
