function results = VAR_TestBlockExo(data,info,exo)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 03/Feb/2017
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Tests block exogeneity of y1t wrt y2t. The null hypothesis
% is no block exogeneity.
% Input:
%   data:
%   -.exo           : Data for y1_t block (T x n_1).
%   -.endo          : Data for y2_t block (T x n_2).
%   info:
%   -.p             : Lag order.
%   -.alpha         : Significance level.
%   exo             : Exogenous variables (optional).
%
% Output:
%   results:
%   -.T             : Number of observations (T-p).
%   -.alpha         : Significance level.
%   -.c_val         : Critical value of likelihood ratio test.
%   -.lik_ratio     : Likelihood ratio test statistic.
%   -.P_val         : P-value for likelihood ratio test.
%   -.LogL          : log likelihood.
%   -.aic           : Akaike information criterion.
%   -.sic           : Schwarz information criterion.
%   -.hqc           : Hannan-Quinn information criterion.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Getting info from the code.
p     = info.p;
alpha = info.alpha;
T     = size(data.endo,1)-p;
% Model estimation.
if exist('exo','var') == 0
    exo = [];
end
results_u = EstimateVAR(data.all,p,exo);
results_r = EstimateBlockVAR(data,p,exo);

% Calculate degrees of freedom (= number of restrictions)
df = size(data.exo,2)*size(data.endo,2)*info.p;

% Number of parameters in longest unrestricted equation (see Enders, 1995)
% Note: set npar = 0 for no df adjustment (see Hamilton, 1994)
npar = results_u.k;
% Calculate test statistic
Sigr      = results_r.Sig;
Sigu      = results_u.Sig;
lik_ratio = (T-npar)*(log(det(Sigr))-log(det(Sigu)));
% Obtain critical value
c_val = chi2inv(1-alpha,df);
% Calculate P-value
P_val = 1-chi2cdf(lik_ratio,df);

% Printing results.
fid = 1;
fprintf(fid,'\n');
fprintf(fid,'**********************************************************\n');
fprintf(fid,'VAR block exogeneity test (see Hamilton, 1994) \n');
fprintf(fid,'Included observations: %1.0f \n',T);
fprintf(fid,'----------------------------------------------------------\n');
fprintf(fid,'Null hypothesis: block exogeneity \n');
fprintf(fid,'----------------------------------------------------------\n');
fprintf(fid,'Significance level (in percent) = %0.0f \n',alpha*100);
fprintf(fid,'Critical value Chi-square distribution = %0.4f \n',c_val);
fprintf(fid,'Test statistic (likelihood ratio) = %0.4f \n',lik_ratio);
fprintf(fid,'P-value = %0.4f \n',P_val);
fprintf(fid,'**********************************************************\n');

% Results
results.T           = T;
results.alpha       = alpha;
results.c_val       = c_val;
results.lik_ratio   = lik_ratio;
results.P_val       = P_val;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%