function smoothTs = probSmooth(trialModOdds, Smooth, windMs, Hz, Lms, N)
% Takes the AMICA model probability timeseries (2 models) and smoothes it
% using a Hamming window over windMs samples. Half backwards in time, half
% forwards in time. 
% -----------------------------------------------------------------------------------------
% Inputs:
%
% trialModOdds = timeseries of AMICA model probabilities (2d)
% Smooth       = indicator for border crossing over trials or not
%       'Within'  = only use samples within each trial itself to smooth
%       'Overall' = allow smoothing to happen across neighbouring trials
% windMs       = how many milliseconds of smoothing
% Hz           = original sampling rate of EEG
% Lms          = length per trial (EEG recording)
% N            = number of trials per strategy (total = N*2)
% -----------------------------------------------------------------------------------------
% Outputs:
%
% smoothTs     = AMICA model probabilities smoothed over time (2d)
% -----------------------------------------------------------------------------------------
% Written by Sven Wientjes at the UvA
% Master Brain & Cognitive Sciences research intern at the Van Maanen lab for Mathematical Psychology and
% Cognitive Modeling (MPCM)

%% Define the starting variables
windHz = round(windMs/1000*Hz);   %Translate to Hz -> check for downsampling
smoothProb = zeros(2, Lms/1000*Hz, N*2);

%Round to nearest odd integer
if mod(windHz,2)<1
    windHz = windHz+1;
end
wind = hamming(windHz)';
halfS = median(1:windHz)-1;

%% Smoothing within trials
if strcmp(Smooth, 'Within')
    for trial = 1:N*2
        for t = 1:Lms/1000*Hz
            for model = 1:2
                %Out of bounds at start
                if t-halfS < 1
                    halfS_t = halfS -1;
                    while t - halfS_t < 1
                        halfS_t = halfS_t-1;
                    end
                    winLen_t = length(t-halfS_t:t+halfS)-1;
                    smoothProb(model, t, trial) = sum(wind(end-winLen_t:end) .* trialModOdds(model,t-halfS_t:t+halfS,trial));
                    %Out of bounds at the end
                elseif t+halfS > length(trialModOdds)
                    halfS_t = halfS -1;
                    while t + halfS_t > length(trialModOdds)
                        halfS_t = halfS_t-1;
                    end
                    winLen_t = length(t-halfS:t+halfS_t);
                    smoothProb(model, t, trial) = sum(wind(1:winLen_t) .* trialModOdds(model,t-halfS:t+halfS_t,trial));
                else
                    smoothProb(model, t, trial) = sum(wind .* trialModOdds(model,t-halfS:t+halfS,trial));
                end
            end
            smoothProb(1,t,trial) = smoothProb(1,t,trial)/sum(smoothProb(:,t,trial));
            smoothProb(2,t,trial) = 1 - smoothProb(1,t,trial);
        end
    end
    smooth_trial_M_p = sum(smoothProb, 2);
    smooth_trial_M_p = reshape(smooth_trial_M_p, 2, N*2);
    smooth_trial_M_p = smooth_trial_M_p ./ sum(smooth_trial_M_p);
    
    smoothTs = smooth_trial_M_p;
    
%% Smoothing over trials
elseif strcmp(Smooth, 'Overall')
    prob_fullexp = reshape(smoothProb, 2, length(smoothProb)*N*2);
    modod_fullexp = reshape(trialModOdds, 2, length(smoothProb)*N*2);
    
     for t = 1:length(prob_fullexp)
        for model = 1:2
            %Out of bounds at start
            if t-halfS < 1
                halfS_t = halfS -1;
                while t - halfS_t < 1
                    halfS_t = halfS_t-1;
                end
                winLen_t = length(t-halfS_t:t+halfS)-1;
                smoothProb_fullexp(model, t) = sum(wind(end-winLen_t:end) .* modod_fullexp(model,t-halfS_t:t+halfS));
           %Out of bounds at the end
            elseif t+halfS > length(modod_fullexp)
                halfS_t = halfS -1;
                while t + halfS_t > length(modod_fullexp)
                    halfS_t = halfS_t-1;
                end
                winLen_t = length(t-halfS:t+halfS_t);
                smoothProb_fullexp(model, t) = sum(wind(1:winLen_t) .* modod_fullexp(model,t-halfS:t+halfS_t));
            else
                smoothProb_fullexp(model, t) = sum(wind .* modod_fullexp(model,t-halfS:t+halfS));
            end
        end
    smoothProb_fullexp(1,t) = smoothProb_fullexp(1,t)/sum(smoothProb_fullexp(:,t));
    smoothProb_fullexp(2,t) = 1 - smoothProb_fullexp(1,t);
     end
    smoothProb_fullexp = reshape(smoothProb_fullexp, 2, Lms/1000*Hz, N*2);
    smooth_trial_M_p_fullexp = sum(smoothProb_fullexp, 2);
    smooth_trial_M_p_fullexp = reshape(smooth_trial_M_p_fullexp, 2, N*2);
    smooth_trial_M_p_fullexp = smooth_trial_M_p_fullexp ./ sum(smooth_trial_M_p_fullexp);
    
    smoothTs = smooth_trial_M_p_fullexp;
end

end