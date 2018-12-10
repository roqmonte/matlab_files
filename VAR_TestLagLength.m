function results = VAR_TestLagLength(data,info,exo)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 29/Jan/2017
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Tests for appropriate lag length in a VAR (with/without
% block exogeneity). Likelihood ratio tests inclusion of lag order j-1 vs. j. 
% H0: j should not be included as additional lag. 
% Input:
%   data:
%   -.all           : Data all variables in the model.
%   -.exo           : Data for y1t block (T x n_1).
%   -.endo          : Data for y2t block (T x n_2).
%   info:
%   -.max_lag       : Maximum lag order to be tested.
%   -.alpha         : Significance level.
%   -.block         : (1) VAR block exogeneity; (0) otherwise (Default: 0).
%   exo             : Exogenous variables (optional).
%
% Output:
%   results:
%   -.nobs          : Number of observations T-max_lag.
%   -.alpha         : Significance level.
%   -.c_val         : Critical value of likelihood ratio test.
%   -.lik_ratio     : Likelihood ratio (1 x z+1 vector).
%   -.P_val         : P-values for likelihood ratio tests (1 x z+1 vector).
%   -.aic           : Akaike information criterion (1 x z+1 vector).
%   -.hqc           : Hannan-Quinn information criterion (1 x z+1 vector).
%   -.sic           : Schwarz information criterion (1 x z+1 vector).
%   -.lags          : Lags tested (1 x z+1 vector).
%   -.table         : Table with results.
%
% Index:
% 1. Initial Setup.
% 2. Estimations and Results.
% 3. Printing results.
% 4. Functions.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 1. Initial Setup.
% Model with/without block exogeneity.
if isfield(info,'block') == 0
    do_block = 0;
else
    do_block = info.block;    
end;
if do_block == 0
    fname = 'EstimateVAR';    
elseif do_block == 1
    fname = 'EstimateBlockVAR';    
end
% Model with/without exo variables.
if exist('exo','var') == 1
    mod_id = 0; add_r  = 1;
else 
    mod_id = 1; add_r = 0; exo = [];
end
% Get input
alpha   = info.alpha;
max_lag = info.max_lag;
[T,n]   = size(data.all);

% Allocate memory
Sig         = zeros(n,n,max_lag+add_r);
lik_ratio   = NaN(1,max_lag+add_r);
P_val       = NaN(1,max_lag+add_r);
LogL        = zeros(1,max_lag+add_r);
aic         = zeros(1,max_lag+add_r);
hqc         = zeros(1,max_lag+add_r);
sic         = zeros(1,max_lag+add_r);
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 2. Estimations and Results.
j0 = 1;
for i = mod_id:max_lag
    % Adjust sample to use same data in each run
    begin = 1+max_lag-i;
    if do_block == 1
        estim_data.exo = data.exo(begin:end,:);
        estim_data.endo = data.endo(begin:end,:);
        estim_data.all = data.all(begin:end,:);
    else
        estim_data = data.all(begin:end,:);
    end;
    % Estimate VAR models.
    if i == 0 && size(exo,2) == 0
        Sig(:,:,j0) = NaN;
        LogL(j0)    = NaN;
        aic(j0)     = NaN;
        sic(j0)     = NaN;
        hqc(j0)     = NaN;    
    else
        eval(['output = ' fname '(estim_data,' num2str(i) ',exo(begin:end,:));']);
        % Store relevant reuslts.
        Sig(:,:,j0) = output.Sig;
        LogL(j0)    = output.LogL;
        aic(j0)     = output.aic;
        sic(j0)     = output.sic;
        hqc(j0)     = output.hqc;
    end
    m(j0) = size(output.F,2) + size(exo,2);
    j0 = j0 + 1;
end;

% Degrees of freedom (number of restrictions)
if do_block == 1
    df = (size(output.A_1(:,:,1),1)^2 + size(output.B_1(:,:,1),1)*size(output.B_1(:,:,1),2) + size(output.B_2(:,:,1),1)^2);
else
    df = (size(output.Y,2)^2);
end

% Critical value LR test
c_val = chi2inv(1-alpha,df);
% LR test statistic and P-value
for i = 2:max_lag+1
    lik_ratio(i) = (T-m(i))*(log(det(Sig(:,:,i-1)))-log(det(Sig(:,:,i))));
    % P-value LR test
    P_val(i) = 1-chi2cdf(lik_ratio(i),df);
    if mod_id == 1 && i == max_lag
        break
    end;
end;

% Results
results.nobs        = T-max_lag;
results.alpha       = alpha;
results.c_val       = c_val;
results.LogL        = LogL;
results.lik_ratio   = lik_ratio;
results.P_val       = P_val;
results.aic         = aic;
results.sic         = sic;
results.hqc         = hqc;
results.lags        = mod_id:max_lag;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 3. Printing results.
% Prepare matrix to be printed
nobs = results.nobs ;
fid = 1;
fprintf(fid,'\n');
fprintf(fid,'*****************************************************************************************\n');
fprintf(fid,'VAR lag order selection criteria (see Canova, 2007) \n');
fprintf(fid,'Included observations: %1.0f \n',nobs);
temp_0 = {'Lag' 'LogL' 'LR' 'P-val' 'AIC' 'SIC' 'HQ'};
temp_1 = [num2cell((0:max_lag)') num2cell(LogL') num2cell(lik_ratio') num2cell(P_val') num2cell(aic') num2cell(sic') num2cell(hqc')];
table = [temp_0; temp_1];
display(table);
fprintf(fid,'LogL  : log likelihood \n');
fprintf(fid,'LR    : Sequential likelihood ratio test statistic \n');
fprintf(fid,'P-val : P-value likelihood ratio test \n');
fprintf(fid,'AIC   : Akaike information criterion \n');
fprintf(fid,'SC    : Schwarz information criterion \n');
fprintf(fid,'HQ    : Hannan-Quinn information criterion \n');
fprintf(fid,'*****************************************************************************************\n');
results.table = table;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%