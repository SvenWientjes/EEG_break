function [nbr, datese, mBIC, lBIC] = structural_break(modelProbs, CI_lv, doBIC)
% This function is a simple procedure to automize getting output from the
% univariate structural break model by J. Bai and P. Perron, as implemented
% by Yohei Yamamoto. It requires the m_break toolbox.
% link
% -----------------------------------------------------------------------------------------
% Inputs:

% modelProbs = row vector of data (in our case AMICA model probabilities or
% probability weights for every trial)
% CI_lv            = confidence level at which to run Sequential Procedure
%       1 = 10%
%       2 = 5%
%       3 = 2.5%
%       4 = 1%
% -----------------------------------------------------------------------------------------
% Outputs:

% nbr       = Number of breaks estimated by Sequential Procedure
% datese = Location of break estimated by Sequential Procedure
% mBIC    = Number of breaks estimated by BIC
% lBIC      = Location of breaks estimated by BIC (Global Optimization Procedure)
% -----------------------------------------------------------------------------------------
% Written by Sven Wientjes at the UvA
% Master Brain & Cognitive Sciences research intern at the Van Maanen lab for Mathematical Psychology and
% Cognitive Modeling (MPCM)

%% Start Structural Break procedure
% Add the m_break folder to the path
addpath('E:\m_break')

yyy = modelProbs';
bigt=length(yyy);
y=yyy;

z=ones(bigt,1);

x=[];

q=1;                      % number of regressors z
p=0;                      % number of regressors x 
m=5;                      % maximum number of structural changes allowed
eps1=.15;                 % value of the trimming (in percentage) for the construction 
                          % and critical values of the supF type tests (used in the 
                          % supF test, the Dmax, the supF(l+1|l) and the sequential 
                          % procedure). If these tests are used, h below should be set 
                          % at int(eps1*bigt). But if the tests are not required, 
                          % estimation can be done with an arbitrary h. There are five 
                          % options: eps1= .05, .10, .15, .20, or .25. For each option, 
                          % the maximal value of m above is: 10 for eps=.05, 8 for 
                          % eps1=.10, 5 for eps1=.15, 3 for eps1=.20, and 2 for eps1=.25.

h=round(eps1*bigt);       % minimal length of a segment (h >= q). Note: if robust=1, h 
                          % should be set at a larger value. 

% The followings are options if p > 0 -------------------------------                          
fixb=0;                   % set to 1 if use fixed initial values for beta
betaini=0;                % if fixb=1, load the initial value of beta
maxi=20;                  % maximum number of iterations for the nonlinear procedure to 
                          % obtain global minimizers
printd=1;                 % set to 1 if want the output from the iterations to be printed
eps=0.0001;               % criterion for the convergence
% --------------------------------------------------------------------
robust=1;                % set to 1 if want to allow for heterogeneity and autocorrelation 
                         % in the residuals, 0 otherwise. The method used is Andrews(1991) 
                         % automatic bandwidth with AR(1) approximation and the quadratic 
                         % kernel. Note: Do not set to 1 if lagged dependent variables are 
                         % included as regressors.
prewhit=1;               % set to 1 if want to apply AR(1) prewhitening prior to estimating 
                         % the long run covariance matrix.
hetdat=1;                % option for the construction of the F tests. Set to 1 if want to
                         % allow different moment matrices of the regressors across segments. 
                         % If hetdat=0, the same moment matrices are assumed for each segment 
                         % and estimated from the ful sample. It is recommended to set 
                         % hetdat=1 if p>0.
hetvar=1;                % option for the construction of the F tests.Set to 1 if want to allow 
                         % for the variance of the residuals to be different across segments. 
                         % If hetvar=0, the variance of the residuals is assumed constant 
                         % across segments and constructed from the ful sample. This option 
                         % is not available when robust=1.  
hetomega=1;              % used in the construction of the confidence intervals for the break 
                         % dates. If hetomega=0, the long run covariance matrix of zu is 
                         % assumed identical across segments (the variance of the errors u 
                         % if robust=0)
hetq=1;                  % used in the construction of the confidence intervals for the break 
                         % dates. If hetq=0, the moment matrix of the data is assumed identical 
                         % across segments.
doglobal=1;              % set to 1 if want to cal the procedure to obtain global minimizers
dotest=1;                % set to 1 if want to construct the supF, UDmax and WDmax tests 
                         % doglobal must be set to 1 to run this procedure.
dospflp1=1;              % set to 1 if want to construct the supF(l+1|l) tests where under
                         % the null the l breaks are obtained using global minimizers. 
                         % doglobal must be set to 1 to run this procedure.
doorder=1;               % set to 1 if want to cal the procedure that selects the number of
                         % breaks using information criteria. doglobal must be set to 1 to 
                         % run this procedure.
dosequa=1;               % set to 1 if want to estimate the breaks sequentialy and estimate 
                         % the number of breaks using supF(l+1|l) test   
dorepart=1;              % set to 1 if want to modify the break dates obtained from the 
                         % sequential method using the repartition method of Bai (1995),
                         % Estimating breaks one at a time. This is needed for the confidence 
                         % intervals obtained with estim below to be valid.
estimbic=1;              % set to 1 if want to estimate the model with the number of breaks 
                         % selected by BIC.
estimlwz=0;              % set to 1 if want to estimate the model with the number of breaks  
                         % selected by LWZ
estimseq=1;              % set to 1 if want to estimate the model with the number of breaks
                         % selected using the sequential procedure
estimrep=0;              % set to 1 if want to estimate the model with the breaks selected
                         % using the repartition method
estimfix=0;              % set to 1 if want to estimate the model with a prespecified number
                         % of breaks equal to fixn set below
fixn=0;

% Get sequential procedure nBreaks & Locations
[nbr, datese]=sequa(m,CI_lv,q,h,bigt,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1);
if isempty(datese)
    datese = 0;
end

if doBIC == 1
% Estimate models (for application of information criteria: BIC)
 if p==0                     %Partial check
    [glob,datevec,bigvec]=dating(y,z,h,m,q,bigt);
 else
     disp('This is a partial structural change model and the following')
     disp('specifications were used:')
     disp(['The number of regressors with fixed coefficients, p: ' num2str(p)])
     disp(['Whether (1) or not (0) the initial values of beta are fixed: ' num2str(fixb)])
     disp(['If so, at the following values: ' num2str(betaini)])
     disp(['The convergence criterion is: ' num2str(eps)])
     disp(['Whether (1) or not (0) the iterations are printed: ' num2str(printd)])
     disp('----------------------')
    [glob,datevec,bigvec]=nldat(y,z,x,h,m,p,q,bigt,fixb,eps,maxi,betaini,printd);
 end

if p==0
     zz=z;         %Partial check
 else
     zz=[z x];
end
ssrzero=ssrnul(y,zz); %Idk check

% Get BIC nBreaks
[mBIC, whatever] = order(ssrzero,glob,bigt,m,q); %Where do ssrzero and globl come from?

% Get BIC Locations
if mBIC ~= 0
    lBIC = datevec(1:mBIC,mBIC);
else
    lBIC= 0;
end

else
    mBIC = 'Not Applicable';
    lBIC = 'Not Applicable';
end

end