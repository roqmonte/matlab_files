function res_pca = PCA_reg(y,x,z,setup,print,labels,dates)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 31/Jan/2020
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Function performs PCA regression. The code standarizes the 
% data  in the first step. Note that by default, loadings are normalized 
% such that L'*L/n = In. 
%
% Use fac*load' to get normalized data
% Use fac*load'*std + mu to get original data
% Use sdata*load/n to get contribuions to factors.
%
% Inputs:
%   y          : Target variable.
%   x          : Independent variables (cte excluced) for pca.
%   z          : Additional regressors.
%   setup:
%   -.k        : Number of factors.
%   -.cv       : k-fold cross-validation, if two cut sample on two (Default: 5).
%   -.mcreps   : Monte-Carlo repetitions for cross-validation (Default: 1).
%   print      : (1) Do charts and print results on screen.
%   labels     : Column vector with labels for dep. and indep variables.
%   dates      : Vector with info for dates: [year,month,freq].
%                Freq:(1) monthly data; (2) quaterly data.
%
% Outputs:
%   res_psl:
%   -.facts     : Factors/scores.
%   -.loadings  : Loadings for x's variables into factor.
%   -.betasf    : Betas to factors (cte term incl first element).
%   -.betas     : Betas to stand. x's variable (cte term incl first element).
%   -.sdata     : Standarized data.
%   -.yyhat     : Standarized target variable and model fit.
%   -.mu        : Sample means.
%   -.sigma     : Sample standard deviations.
%   -.varxp     : % var. explained by each factor.
%   -.cs_vxp    : Cumulative sum var. explained by each factor.
%   -.kopt      : Optimal# of factors
%   -.MSE       : Matris with MSE from cross-validation and MC repetitios.
%   -.ols       : Results from OLS regression.
%
% Index.
% 1. Setup and PCA reg estimation
% 2. Functions.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 1. Setup and PCA reg estimation
% Checking unput from code
if isfield(setup,'cv') == 0
    setup.cv = 5;
end
if isfield(setup,'mcreps') == 0
    setup.mcreps = 1;
end
% Checking labels.
if exist('labels','var') == 0 || isempty(labels)
    labels = [];
end
% Checking print.
if exist('print','var') == 0
    print = 0;
end
% Checking dates.
if exist('dates','var') == 0 || isempty(dates)
    dates = [1900,1,1];
end
    
% Standarize the data
data = [y x];
T    = size(data,1);
[sdata,mu,sigma] = standardise(data);

% Standarize z variables
nz = size(z,2);
if nz == 0
    z   = [];
    sz  = [];
    zmu = [];
    zsg = [];
else
    [sz,zmu,zsg] = standardise(z);
end

% Model selection
% Simple PCA reg estimation
if setup.cv == 1
    % Do PCA and keep k factors.
    res_aux = PCAf(sdata(:,2:end),0);
    facts   = res_aux.pca(:,1:setup.k);
    loads   = res_aux.loads(:,1:setup.k);
    % Do regression.
    res_temp = OLSest(sdata(:,1),[ones(T,1) facts sz],'OLS',[],print,labels,dates);
    betasf  = res_temp.b;
    kopt = setup.k;
    
% Cross-validation and Monte-carlo replications.
else
    % Cross-validation, no Monte-Carlo repetitions
    if setup.mcreps == 1
        [facts,loads,betasf,kopt,MSE] = cv_kfold(sdata,sz,setup.cv,setup.k);
        
    % Cross-validation and Monte-Carlo repetitions
    else
        h_wait = waitbar(0,'Doing cross-validation with MC replications...');
        % Do MC repititions.
        MSE = [];
        for i0 = 1:setup.mcreps
            [~,~,~,~,aux] = cv_kfold(sdata,sz,setup.cv,setup.k);
            MSE = [MSE aux];
            waitbar(i0/setup.mcreps,h_wait);
            clear aux;
        end
        close(h_wait);
        % Find optimal k
        [~,kopt] = min(mean(MSE,2));
        kopt = kopt - 1;
        
        % Do PCA and keep k factors.
        if kopt > 0
            res_aux = PCAf(sdata(:,2:end),0);
            facts   = res_aux.pca(:,1:kopt);
            loads   = res_aux.loads(:,1:kopt);
            % Do regression.
            res_temp = OLSest(sdata(:,1),[ones(T,1) facts sz],'OLS');
            betasf  = res_temp.b;
        else
            facts = 0;
            loads = 0;
            betasf= 0;
        end
    end
    % Saving additional results.
    res_pca.kopt = kopt;
    res_pca.MSE  = MSE;    
end


% Do table/charts
if kopt > 0
    res_temp = OLSest(sdata(:,1),[ones(T,1) facts sz],'OLS',[],print,labels,dates);
    % Getting contribution from each x's into target variable.
    n = size(loads,1);
    if nz == 0
        betas  = [betasf(1); loads*betasf(2:end)/n];
    else
        betas  = [betasf(1); loads*betasf(2:end-nz)/n; betasf(end-nz+1:end)];
    end

    % Results.
    res_pca.facts    = facts;
    res_pca.loadings = loads;
    res_pca.betasf   = betasf;
    res_pca.betas    = betas;
    res_pca.sdata    = sdata;
    res_pca.yyhat    = [res_temp.dta(:,1) res_temp.yhat];
    res_pca.mu       = [mu zmu];
    res_pca.sigma    = [sigma zsg];
    % Fraction of the variance explained by each factor
    temp1 = (var(facts) .* (res_pca.betasf(2:end-nz,1).^2)')';
    temp2 = 1 - sum(temp1);
    res_pca.varxp   = [temp1; temp2];
    res_pca.cs_vxp  = cumsum(res_pca.varxp);
    res_pca.ols = res_temp;
else
    res_pca = [];
end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 2. Functions.

% Funtion perform k-fold cross-validation.
function [facts,loads,betasf,kopt,MSE] = cv_kfold(data,sz,cv,k)
% Inputs: 
%   data    : Data.
%   sz      : Additional regressors (z-score).
%   cv      : k-fold cross-validation.
%   k       : Number of factors.
% Outputs:
%   facts   : PSL fators.
%   loads   : Loadings (x into factors).
%   betas   : Target variable betas to factors.
%   kopt    : Optimal number of factos.
%   MSE     : Matrix with MSE from cross-validation.

% Sample size and xs varianles.
T  = size(data,1);
nz = size(sz,2);
% Selection variable
if cv > 2
    id1 = randi([1,cv],T,1);
else
    temp = randi([round(T*0.2) round(T*0.8)],1,1);
    id1 = [ones(temp,1); 2*ones(T-temp,1)];
    clear temp
end
id2 = 1:cv;

% Do PCA and keep k factors.
res_aux = PCAf(data(:,2:end),0);
facts   = res_aux.pca(:,1:k);
% New data, target variable and factors.
data0   = [data(:,1) facts sz];

% Do cross-validation
for i0 = 1:cv
    % Variable selection
    temp_test = (id1==id2(i0));
    temp_est  = 1-temp_test;
    
    % Training and testing sample
    tr_data = data0(temp_est==1,:);
    ts_data = data0(temp_est==0,:);
    % Regressors.
    if nz == 0
        tr_x = tr_data(:,2:end);
        ts_x = ts_data(:,2:end);
    else
        tr_x = tr_data(:,2:end-nz);
        ts_x = ts_data(:,2:end-nz);
        tr_z = tr_data(:,end-nz+1:end);
        ts_z = ts_data(:,end-nz+1:end);
    end
    
    % Sample size trainin sample.
    Tn = size(tr_data,1);
    
    % Loop over the number of factors.
    for i1 = 0:k
        if i1 == 0
            % Model evaluation.
            uhat = (ts_data(:,1) - mean(tr_data(:,1)));
        else
            % Model estimation & evaluation
            if nz ==  0 
                res_aux = OLSest(tr_data(:,1),[ones(Tn,1) tr_x(:,1:i1)],'OLS');
                uhat = (ts_data(:,1) - [ones(size(ts_data,1),1) ts_x(:,1:i1)]*res_aux.b);
            else
                res_aux = OLSest(tr_data(:,1),[ones(Tn,1) tr_x(:,1:i1) tr_z],'OLS');
                uhat = (ts_data(:,1) - [ones(size(ts_data,1),1) ts_x(:,1:i1) ts_z]*res_aux.b);
            end
        end

        % Matrix with MSE.
        MSE(i1+1,i0)  = uhat'*uhat/Tn;
    end
    % Cleaning memory
    clear temp_test temp_est tr_data ts_data tr_cte ts_cte i1 res_aux uhat;
end
% Choosing optimal number of factors in the model.
[~,kopt] = min(mean(MSE,2));
kopt = kopt - 1;
if kopt > 0
    % Do PCA and keep k factors.
    res_aux = PCAf(data(:,2:end),0);
    facts   = res_aux.pca(:,1:kopt);
    loads   = res_aux.loads(:,1:kopt);
    % Do regression.
    res_temp = OLSest(data(:,1),[ones(T,1) facts sz],'OLS');
    betasf  = res_temp.b;
else
    facts = 0;
    loads = 0;
    betasf= 0;
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%