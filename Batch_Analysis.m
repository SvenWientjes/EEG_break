%% Batch Analysis
% This script analyzes the data of AMICA model probabilities by running
% them through a univariate structural break detection model. The model is
% developed by P. Perron and J. Bai. It is implemented in the m_break
% package by Yohei Yamamoto. The m_break package can be downloaded at:
% http://people.bu.edu/perron/code/m-break-matlab.zip
% and needs to be added to the Matlab Path for this script to function.
% -----------------------------------------------------------------------------------------
% Inputs:
%
% Generative parameters (e.g. Hz, Lms, etc.) must be identical to those
% used in generating the data using Get_Run_Full.
%
% The path at the start of the 'run = Runs' for-loop must lead to the
% correct folder containing the Run_i folders with AMICA output.
%
% Runs   = Which runs do you wish to combine in the eval output table
% Smooth = Sets the type of probability smoothing
%      'Within'  = Not crossing over trial borders
%      'Overall' = Allow crossing over trial borders
% windMs = How many samples around any point are used for smoothing (Hamming Window)
% 
% -----------------------------------------------------------------------------------------
% Output: 
%
% eval = Table containing evaluative statistics of the break detection
% procedure. Rows are 'participants'.
%         Seq_N = Number of breaks estimated by sequential procedure
%         Seq_L = Location of breaks estimated by sequential procedure
%         Bic_N = Number of breaks estimated by Bayesian Information Criterion
%         Bic_L = Location of breaks estimated by Bayesian Information Criterion
% -----------------------------------------------------------------------------------------
% Written by Sven Wientjes at the UvA
% Master Brain & Cognitive Sciences research intern at the Van Maanen lab for Mathematical Psychology and
% Cognitive Modeling (MPCM)

%% Load toolboxes
addpath('E:\eeglab14_1_2b\plugins\amica1.5')
%% Start Batch Analysis
% Settings for analysis
Runs   = [1:30];
Smooth = 'Overall'; % Can be 'Within' or 'Overall'
windMs = 5000;

Hz      = 500;
N       = 50;
Lms     = 3500;
nLoc    = 7;
nChange = 6;
SNR     = 2;

% Cell array for Run names in table
for i = 1:length(Runs)
    RunNames{i} = ['Run_', int2str(Runs(i))];
end
% Preallocate table for keeping evaluation statistics
eval = table(zeros(length(Runs),1),zeros(length(Runs),1),zeros(length(Runs),1),zeros(length(Runs),1), ...
    'VariableNames', {'Seq_N', 'Seq_L', 'Bic_N', 'Bic_L'},...
    'RowNames', RunNames);
binSeqL = [];
binBicL = [];

for run = Runs
    Path = ['E:\SNR-2_Tr-50_L-7_S-abo\Run_', int2str(run), '\Run_', int2str(run)];
    
    % Load Data
    Model = loadmodout15(Path);
    
    % Convert Log10 odds to probability
    modOdds = 10.^Model.v;
    
    % Shape into (Models x Samples x Trials)
    trialModOdds = reshape(modOdds, Model.num_models, Lms/1000*Hz, N*2);
    
    % Run Probability Smoother
    smoothOdds = probSmooth(trialModOdds, Smooth, windMs, Hz, Lms, N);
    
    % Run Break Detection
    breakTest = smoothOdds(1,:);
    testPoint = length(breakTest)/2;
    
    [seqN, seqL, bicN, bicL] = Structural_Break(breakTest, 1, 1);
    
    % Allocate break statistics to eval table
    eval.Seq_N(['Run_', int2str(run)]) = seqN;
    eval.Seq_L(['Run_', int2str(run)]) = seqL(1);
    eval.Bic_N(['Run_', int2str(run)]) = bicN;
    eval.Bic_L(['Run_', int2str(run)]) = bicL(1);
    
    % For more than 1 detected break
    if length(seqL)>1
        binSeqL = [binSeqL, seqL(2:end)];
    end
    if length(bicL)>1
        binBicL = [binBicL, bicL(2:end)];
    end

end

writetable(eval, 'E:\SNR-2_Tr-50_L-7_S-abo\EvalStat.txt', 'Delimiter', ',')
