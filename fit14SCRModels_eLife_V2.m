%% These functions fit the 9 models (8 versions of Rescorla-Wagner(RW) and hybrid Pearce-Hall(PH) to SCR data for Abend et al., 2021


% Below is an example of how one would fit 8 RW models to a dataset:


% For the 8 RW models, a matrix was used to define the presence of
% additional model features.  A 1 indicates that the model feature is
% included and a 0 indicates that it is not included

% Column 1 = habituation params
% Column 2 = bayesian learning rate
% Column 3 = learning inertia

% model = [0 0 0; 
%           0 0 1;
%           0 1 0;
%           1 0 0;
%           0 1 1;
%           1 0 1;
%           1 1 0;
%           1 1 1];

% A single model was chosen for each run of the fit. The function fitSCRGen
% takes the model information as an input (e.g.: modelinfo=model(1,:) would
% be the simplest RW model, model 1)

% Then use fminsearch to minimize the negative log-likelihood of the fit to
% the model (set options for limits on fminsearch):

% [fitparams, ll] = fminsearch(@(params) fitSCRGen(params, scr, us, csType, vi, modelinfo, n, mu), params, options);

%inputs:
%params: starting values for parameters
%scr: skin conductance response to conditioned stimulus
%us: skin conductance response to unconditioned stimulus
%vi: initial value estimate
%modelinfo: 3 element vector defining which model should be fit

%outputs: 
%fitparams: values of the parameters from the particular fit
%ll: log likelihood of fit

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Model Functions %%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wrapper function for fitSCRfGen, for compatibility with fminsearch :
function [ll] = fitSCRGen(params, scr, us, csType, vi, modelV, n, mu)

[ll, ~] = fitSCRfGen(params, scr, us, csType, vi, modelV, n, mu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For fitting 8 RW models to acquisition data:
function [ll, v] = fitSCRfGen(params, scr, us, csType, vi, modelV, n, mu)

v   = zeros(2, length(scr) + 1);
v(1:2, 1) = vi;
rpe = zeros(length(scr), 1);

for trial = 1 : length(scr)
    rpe(trial) = us(trial) - v(csType(trial), trial);
    
    if modelV(2) == 1
        lrp1 = params(1)/(trial^1/2);
        lrp2 = params(2)/(trial^1/2);
    else
        lrp1 = params(1);
        lrp2 = params(2);
    end
    
    if modelV(3) == 1
        ftrial = min(trial-1,2);
        if trial > 1
            srpe = sum(rpe(trial:-1:trial-ftrial));
        else
            srpe = rpe(trial);
        end
    else
        srpe = rpe(trial);
    end
    
    if csType(trial) == 1
        v(1, trial+1) = v(1, trial) + lrp1*srpe + mu(trial);
        v(2, trial+1) = v(2, trial) + lrp2*srpe + mu(trial);
    else
        v(1, trial+1) = v(1, trial) + mu(trial);
        v(2, trial+1) = v(2, trial) + mu(trial);
    end
end

z1 = find(csType == 1);
z2 = find(csType == 2);

if modelV(1) == 1
    exp1 = params(3);
    exp2 = params(4);
    
    if isnan(n)
        vec = [0 0];
        n = 2;
    else
        vec = zeros(n,1)';
    end
    
    v(1, z1) = v(1, z1).*exp(-exp1*([vec 1:(length(z1) - n)]));
    v(2, z2) = v(2, z2).*exp(-exp2*([vec 1:(length(z2) - n)]));
end

ll = 0;

for trial = 1 : length(scr)
    ll = ll + abs(v(csType(trial), trial) - scr(trial)).^2;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Wrapper function for fitSCRfGenExt, for compatibility with fminsearch :
function [ll] = fitSCRGenExt(params, scr, us, csType, vi, modelV, mu)

    
[ll, ~] = fitSCRfGenExt(params, scr, us, csType, vi, modelV, mu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For fitting 8 RW models to extinction data:
function [ll, v] = fitSCRfGenExt(params, scr, us, csType, vi, modelV, mu)

v   = zeros(2, length(scr) + 1);
v(1:2, 1) = vi;
rpe = zeros(length(scr), 1);

for trial = 1 : length(scr)
    rpe(trial) = us(trial) - v(csType(trial), trial);
    
    if modelV(2) == 1
        lrp1 = params(1)/(trial^1/2);
        lrp2 = params(2)/(trial^1/2);
    else
        lrp1 = params(1);
        lrp2 = params(2);
    end
    
    if modelV(3) == 1
        ftrial = min(trial-1,1);
        if trial > 1
            srpe = sum(rpe(trial:-1:trial-ftrial));
        else
            srpe = rpe(trial);
        end
    else
        srpe = rpe(trial);
    end
    
    csn_seq = find(csType == 2);
    first_csn = csn_seq(1);
    if trial < first_csn || first_csn == 1
        if csType(trial) == 1
            v(1, trial+1) = v(1, trial) + lrp1*srpe + mu(trial);
            v(2, trial+1) = v(2, trial) + mu(trial);
        else
            v(1, trial+1) = v(1, trial) + mu(trial);
            v(2, trial+1) = v(2, trial) + lrp2*srpe + mu(trial);
        end
    else
        if csType(trial) == 1
            v(1, trial+1) = v(1, trial) + lrp1*srpe + mu(trial);
            v(2, trial+1) = v(2, trial) + lrp2*srpe + mu(trial);
        else
            v(1, trial+1) = v(1, trial) + mu(trial);
            v(2, trial+1) = v(2, trial) + mu(trial);
        end
    end
end

z1 = find(csType == 1);
z2 = find(csType == 2);

if modelV(1) == 1
    exp1 = params(3);
    exp2 = params(4);
    
    v(1, z1) = v(1, z1).*exp(-exp1*([0 1:(length(z1)-1)]));
    v(2, z2) = v(2, z2).*exp(-exp2*([0 1:(length(z2)-1)]));
end

ll = 0;

for trial = 1 : length(scr)
    ll = ll + abs(v(csType(trial), trial) - scr(trial)).^2;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For fitting Pearce-Hall model to data:
function [ll, v] = fitSCRfPHG(params, scr, us, csType, vi, mu)

v(1:2,1) = vi;

alphpos = 1;
alphneg = 1;

for trial = 1 : length(scr)
    
    rpe = us(trial) - v(csType(trial), trial);
    
    if csType(trial) == 1
        v(1, trial+1) = v(1, trial) + params(1)*alphpos*(rpe) + mu(trial);
        v(2, trial+1) = v(2, trial) + params(2)*alphneg*(rpe) + mu(trial);
        
        alphpos = (1-params(3))*alphpos + abs(rpe)*params(3);
        alphneg = (1-params(4))*alphneg + abs(rpe)*params(4);
    else
        v(1, trial+1) = v(1, trial) + mu(trial);
        v(2, trial+1) = v(2, trial) + mu(trial);
    end
    
end

ll = 0;

for trial = 1 : length(scr)
    ll = ll + abs(v(csType(trial), trial) - scr(trial)).^2;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For fitting Li et al., models
function [ll, v, alpha_all, vals_all, a_pos, a_neg, warningTracker, b] = fitSCRfPHG2(params, scr, us, csType, vi, mu, modelinfo)
%saves associability parameters trial by trial 9.27.21. DB
%modified to match Li et al., 2011 models 02.02.2022
v(1:2,1) = vi;

alphpos = 1;
alphneg = 1;
a_pos=zeros(length(scr),1);
a_neg=zeros(length(scr),1);
alpha_all=zeros(length(scr),1);

for trial = 1 : length(scr)
        
    %RPE(trial) = us(trial) - v(csType(trial), trial); %added for troubleshooting
    
    %rpe = us(trial) - v(csType(trial), trial); %original
    
    %making us binary 02.15.2022 for Li et al. models:
    if us(trial)>0
        r=1;
    else
        r=0;
    end
    
    rpe = r - v(csType(trial), trial);
    
    if csType(trial) == 1
        v(1, trial+1) = v(1, trial) + params(1)*alphpos*(rpe) + mu(trial); %params(1)=kappa for cs+
        v(2, trial+1) = v(2, trial) + params(2)*alphneg*(rpe) + mu(trial); %params(2)=kappa for cs-
        
        alphpos = (1-params(3))*alphpos + abs(rpe)*params(3); %params(3)=gamma cs+
        alphneg = (1-params(4))*alphneg + abs(rpe)*params(4); %params(4)=gamma cs-
    else
        v(1, trial+1) = v(1, trial) + mu(trial);
        v(2, trial+1) = v(2, trial) + mu(trial);
    end
    
    a_pos(trial,1)=alphpos; %trial by trial update of each alpha for cs+ and cs-
    a_neg(trial,1)=alphneg;
end

ll = 0;

%Interleaving Vs and Alphas for CS+ and CS- trials to do just 1 regression :
for trial = 1 : length(scr)
    if csType(trial)==2
        alpha_all(trial,1)=a_neg(trial);
    elseif csType(trial)==1
        alpha_all(trial,1)=a_pos(trial);
    end
    vals_all(trial,1)=v(csType(trial), trial);
    
end

%with regression for both (use one of these at a time):
if modelinfo==1
    [b, bint, r]=regress(scr, [ones(length(scr),1), vals_all]); %V
    %return betas separately
elseif modelinfo==2
    [b, bint, r]=regress(scr, [ones(length(scr),1), alpha_all]); %a
elseif modelinfo==3
    [b, bint, r]=regress(scr, [ones(length(scr),1), alpha_all, vals_all]); %V+a
end
ll=sum(r.^2); %minimize the square of the residuals if using regression in function


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For fitting Tzovara et al., models
%wrapper function for compatibility with fminsearch:
function [ll] = fitSCRUnc2(params, scr, us, ucData, csType, vi, habituationFlag)

[ll, ~] = fitSCRfUnc2(params, scr, us, ucData, csType, vi, habituationFlag);
end

function [ll, yhat] = fitSCRfUnc2(params, scr, scr_all, ucData, csType, vi, habituationFlag)

%Tzovara models with and without habituation
%updated 4.14.22 to only have one set of betas

%scr is only CS+ trials (length 10)
%scr_all is CS+ and CS- trials (length 20)
%csType is length 20 for all trials 1=CS+ 2=CS-

if habituationFlag==1
    exp1=params(3); %CS+ habituation    
    exp2=params(4); %CS- habituation
end

%CS+ params
v_CSP   = zeros(1, length(scr) + 1);
E_theta_CSP = zeros(1,length(scr) + 1);
v_CSP(1) = 0;
beta_one = params(1);
beta_knot = params(2);
non_us = ucData == 0; %ucData is binary so for non_us 1=CS+ reinforced with UCS 0=not reinforced 80/20 schedule
yhat_CSP = zeros(length(scr)+1,1);
alpha_CSP = zeros(length(scr),1); alpha_CSP(1) = 1;
beta_CSP = zeros(length(scr),1); beta_CSP(1) = 1;

%CS- params
v_CSN   = zeros(1, length(scr) + 1);
E_theta_CSN = zeros(1,length(scr) + 1);
v_CSN(1) = 0;
non_stim=ones(1,length(ucData)); %cs- never reinforced
yhat_CSN = zeros(length(scr)+1,1);
alpha_CSN = zeros(length(scr),1); alpha_CSN(1) = 1; %not sure if these should be initialized at 1 or 0 (Bayesian learner in Tzovara 2018)
beta_CSN = zeros(length(scr),1); beta_CSN(1) = 1;

%CS+ (csType=1), compute alpha and beta:
for trial = 1 : length(scr)  
    alpha_CSP(trial+1) = alpha_CSP(trial) + ucData(trial); %goes up every time ucs was delivered
    beta_CSP(trial+1) = beta_CSP(trial) + non_us(trial); %beta doesn't change when ucs delivered, goes up by 1 otherwise

    E_theta_CSP(trial + 1) = (alpha_CSP(trial))/(alpha_CSP(trial) + beta_CSP(trial));
    v_CSP(trial+1) = -log(alpha_CSP(trial) + beta_CSP(trial));
end

%for CS- only increment beta, it was never reinforced with ucs
%alpha stays at 1 but beta increments up
%CS- (csType=2)
for trial = 1 : length(scr)  
    alpha_CSN(trial+1) = 1; %cs- never reinforced
    beta_CSN(trial+1) = beta_CSN(trial) + non_stim(trial); %just increments every trial
    E_theta_CSN(trial + 1) = (alpha_CSN(trial))/(alpha_CSN(trial) + beta_CSN(trial));
    v_CSN(trial+1) = -log(alpha_CSN(trial) + beta_CSN(trial));
end


if habituationFlag==1
    %add habituation to CS+
    
    vec = [0 0]; n = 2; %don't habituate on first 2 trials?
    %vec = zeros(n,1)'; %other option for this
    
    v_CSP = v_CSP.*exp(-exp1*([vec 1:(length(v_CSP) - n)])); %habituation CS+
    %v=v*exp(alpha(k-ko)^t) ko=2. ko not fit to individual subjects
    
    v_CSN = v_CSN.*exp(-exp2*([vec 1:(length(v_CSN) - n)])); %habituation CS-
end

ll = 0;

%compute yhat for CS+ and CS- combined (modified to regress all at once)

% %original, separate regression:
% for trial = 1:length(scr)
%     z_CSP(trial) = v_CSP(trial) + E_theta_CSP(trial);
%     z_CSN(trial) = v_CSN(trial) + E_theta_CSN(trial);
%     yhat_CSP(trial) = beta_one_CSP * z_CSP(trial) + beta_knot_CSP;
%     yhat_CSN(trial) = beta_one_CSN * z_CSN(trial) + beta_knot_CSN;
% end
% 
% yhat_CSP=yhat_CSP(1:10);
% yhat_CSN=yhat_CSN(1:10);
% 
% %combine CSP and CSN yhats into one yhat, properly interleaved:
% csp_idx=find(csType==1);
% csn_idx=find(csType==2);
% yhat(csp_idx,1)=yhat_CSP;
% yhat(csn_idx,1)=yhat_CSN;

%%%%%%%%% new version with CS+ and CS- regression together:
for trial = 1:length(scr)
    z_CSP(trial) = v_CSP(trial) + E_theta_CSP(trial);
    z_CSN(trial) = v_CSN(trial) + E_theta_CSN(trial);
end
z_CSP=z_CSP(1:10);
z_CSN=z_CSN(1:10);

csp_idx=find(csType==1);
csn_idx=find(csType==2);

z_CSP_CSN(csp_idx)=z_CSP;
z_CSP_CSN(csn_idx)=z_CSN;

yhat= beta_one*z_CSP_CSN + beta_knot; %single regression


for trial=1:length(scr_all)    
    %calculate ll
    ll = ll + ((yhat(trial) - scr_all(trial)).^2);
end

yhat(1) = vi;
end


