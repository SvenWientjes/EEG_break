function [EEGmat, infoStruct] = simulate_eeg(Hz, N, Lms, nLoc, nChange, SNR)
% This function generates a non-staionary EEG signal (different ERPs on
% different trials) - the EEG noise is simulated using the process of Rafal
% Bogacz and Nick Yeung from 2002. The EEG contains 64 electrodes.
% -----------------------------------------------------------------------------------------
% Inputs:
%  Hz      = Sampling rate of the simulated EEG
%  N       = How many trails per strategy (for now only two strategies exist, so
%            total number of trials is always N*2)
%  Lms     = How many milliseconds per trial (every trial is equal length)
%  nLoc    = How many ERP sources are active per trial
%  nChange = How many of these ERP sources are different between strategies
% -----------------------------------------------------------------------------------------
% Outputs:
%  EEGmat     = An (64 by Hz*Lms/1000 by N*2) Matrix containing the signals
%               with EEG noise over them.
%  
%  infoStruct = A struct containing info about the process which simulated
%               EEGmat and contains all chosen input variables and the following additional
%               variables:
%      ConcMat : A 2D matrix of the full experiment: ERP signals with
%                sensor noise
%      CompMat : A matrix of the pure ERP signals, not containing any EEG or sensor noise
%      forMat  : A matrix of linear weights describing how each pure ERP
%               signal propagates to the electrodes
%      erp     : A struct describing the simulation process of
%               strategy 1 pure ERPs
%      erp_2   : A struct describing the simulation process of strategy
%               2 pure ERPs
%      erp_stay: Index describing with ERP processes of strategy 1 are
%                'carried over' to strategy 2.
% -----------------------------------------------------------------------------------------
% Function depends on the SEREEGA toolbox for simulating ERPs. This
% includes a dependency on EEGlab. Make sure these are properly installed.
% -----------------------------------------------------------------------------------------
% Written by Sven Wientjes at the UvA
% Master Brain & Cognitive Sciences research intern at the Van Maanen lab for Mathematical Psychology and
% Cognitive Modeling (MPCM)
% -----------------------------------------------------------------------------------------
% Some of the paths used in this function must be rewritten to fit the
% user's computer, these paths are here:

% Path to 'new york head' lead field
nyPath = 'C:\Program Files\MATLAB\R2018b\toolbox\SEREEGA-master\sa_nyhead.mat';

% Path to Yeung & Bogacz noise simulator
addpath('C:\Users\wient\OneDrive\Documents\Uni-SvenStudie\Brain & Cognition year 1\Stage\NoiseSim')

% Path to EEGlab adminfunc folder
addpath('C:\Program Files\MATLAB\R2018b\toolbox\eeglab14_1_2b\functions\adminfunc')

%% Start the simulation process
    
% Initialize epochs
epochs        = struct();
epochs.n      =      N;   %Trials
epochs.srate  =  Hz;      %Hz
epochs.length = Lms;      %Trial in ms

% Load the Lead Field (lf)
matfile(nyPath);
lf = lf_generate_fromnyhead('montage', 'S64'); %64 electrode setup

% Draw source locations ; at least 50mm apart
sourcelocs  = lf_get_source_spaced(lf, nLoc, 50);

% Simulate the pure ERPs, with:
% 1, 2 or 3 peaks // centered between 200 - Lms minus 200 ms // 25 to 200 ms wide
% ampliture between -1 and 1 (with specific spacing)
erp = erp_get_class_random([1:3], [200:Lms-200], [25:200], ...
		[-1:.1:-.5, .5:.1:1], 'numClasses', nLoc, 'latencyRelDvs', [0.01:0.02:0.3],...
        'widthRelDvs', [0.01:0.02:0.3], 'amplitudeRelDvs', [0.01:0.02:0.3]);
    
% Make a 'brain component'
c = utl_create_component(sourcelocs, erp, lf);

% Make scalp data (TestMat) and the pure ERPs (CompMat_1) of strategy 1
[TestMat CompMat_1] = generate_scalpdata(c, lf, epochs, 'sensorNoise', 0.0001);

% _______________________________STRATEGY 2_______________________________%
if nChange > 0
    % Independently draw new source locations: As many 'new' ones as new ERPs
    % will be made for strategy 2 (nChange)
    sourcelocs_2 = lf_get_source_spaced(lf, nChange, 50);
    
    % Simulate new pure ERP waveforms: also nChange // same details as strategy
    % 1 erp
    erp_2 = erp_get_class_random([1:3], [200:Lms-200], [25:200], ...
        [-1:.1:-.5, .5:.1:1], 'numClasses', nChange, 'latencyRelDvs', [0.01:0.02:0.3],...
        'widthRelDvs', [0.01:0.02:0.3], 'amplitudeRelDvs', [0.01:0.02:0.3]);
    
    erp_stay = 0;
    
    % Only concatenate stuff if not everything changes
    if nChange ~= nLoc
        %Sample which old ERPs stay active -> index of strategy 1 erp struct
        erp_stay = datasample(1:nLoc, nLoc-nChange, 'Replace', false);
    
        %Concatenate new and old ERPs
        erp_2 = [erp(erp_stay) erp_2];
    
        %Concatenate new and old source locations
        sourcelocs_2 = [sourcelocs(erp_stay) sourcelocs_2];
    end
    
    %Get the new brain components
    c2_1 = utl_create_component(sourcelocs_2, erp_2, lf);
    
    %Get the sensor data (TestMat_2) and the pure ERP waveforms (CompMat_2) of the new brain components
    [TestMat_2 CompMat_2] = generate_scalpdata(c2_1, lf, epochs, 'sensorNoise', 0.0001);
    
    %Concatenate with the old sensor data -> now we have a structural break!
    TestMat_Full = cat(3, TestMat, TestMat_2);
    
    %Concatenate trials to get 2D only time/electricity
    ConcMat = reshape(TestMat_Full, 64, Hz*Lms/1000*N*2);                                             %Sensor data
    CompMat = cat(3, CompMat_1, CompMat_2);                                                                  %Pure ERP waves

% If we wish not to have a strategy shift:
elseif nChange == 0
    % Copy sourcelocs to sourcelocs_2 and c to c_2 for forMat // copy erp
    % to erp_2
    sourcelocs_2 = sourcelocs;
    c2_1 = c;
    erp_2 = erp;
    erp_stay = 0;
    % Extend the scalp data of the same ERPs
    [TestMat_2 CompMat_2] = generate_scalpdata(c, lf, epochs, 'sensorNoise', 0.0001);
    
    %Concatenate all the scalp data
    TestMat_Full = cat(3, TestMat, TestMat_2);
    
    %Reshape to 2D to get only electricity over time
    ConcMat = reshape(TestMat_Full, 64, Hz*Lms/1000*N*2);                                            %Sensor data
    CompMat = cat(3, CompMat_1, CompMat_2);                                                                 %Pure ERP waves
end

%% Add separate Yeung noise to every channel
ConcMat = ConcMat/max(max(abs(ConcMat))) * 20 * SNR;      %Scale to maximum 1 -> maximum 20 (or 10 for SNR = 0.5, or 40 for SNR = 2)
EEGmat = ConcMat;

for i = 1:64
    noiseMat = noise(Lms/1000*Hz, N*2, Hz);
    noiseMat = noiseMat/max(abs(noiseMat))*20;           %Scale to maximum 1 -> maximum 20
    EEGmat(i,:) = ConcMat(i,:) + noiseMat;
end

%% Extract the forward matrix (weight propagations) of each ERP component
forMat = zeros(64, nLoc,2);
allSourcelocs = [sourcelocs;sourcelocs_2];
all_c = [c c2_1];
for j = 1:2
    for i = 1:nLoc
        forMat(:,i,j) = reshape(lf.leadfield(:, allSourcelocs(j,i), :), 64, 3) * all_c(i+nLoc*(j-1)).orientation'; %Add nComp if second run
    end
end

%% Create the information struct

% Initialize struct
infoStruct = struct();
infoStruct.Hz = Hz;
infoStruct.N = N;
infoStruct.Lms = Lms;
infoStruct.nLoc = nLoc;
infoStruct.nChange = nChange;
infoStruct.forMat = forMat;
infoStruct.ConcMat = ConcMat;
infoStruct.CompMat = CompMat;
infoStruct.erp = erp;
infoStruct.erp_2 = erp_2;
infoStruct.erp_stay = erp_stay;
infoStruct.SNR = SNR;

end
