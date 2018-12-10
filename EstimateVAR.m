function results = EstimateVAR(data,p,exo)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 29/Jan/2017
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Reduced-form VAR estimation by OLS.
% Inputs:
%   data        : Data.
%   p           : Lag order.
%   exo         : Exogenous variables (optional).
%
% Outputs:
%   results:
%   -.C         : Coefficients exogenous variables.
%   -.A         : Coefficients lagged variables.
%   -.F         : Companion form.
%   -.fitted    : Fitted values (T x n matrix).
%   -.resid     : Regression residuals (T x n matrix).
%   -.Sig       : Error covariance matrix.
%   -.LogL      : Log-likelihood
%   -.aic       : Akaike information criterion.
%   -.sic       : Schwarz information criterion.
%   -.hqc       : Hannan-Quinn information criterion.
%   -.Y         : Dependent variables of the system (lhs variables).
%   -.X         : Independent variables of the system (rhs variables).
%   -.k         : Number of  parameters per equation.
%   -.kn        : Total Number parameters of the model.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Getting info for estimation.
if exist('exo','var') == 0
    exo = [];
end

% Set up dependent variables
Y = data(1+p:end,:);
% Generating lags of dependent variables.
X = [];
i = 1;
while i <= p
    X = [X data(1+p-i:end-i,:)];
    i = i + 1;
end;

% Adding exogenous variables
X = [X exo(1+p:end,:)];

% Var estimation.
n = size(data,2);
T = length(data)-p;
coef = (X'*X\X'*Y)';    % OLS equation per equation.
E    = (Y-X*coef');     % Computing residuals.
Sig  = (E'*E)/(T-n);    % Var/cov matrix.

% Total number of parameters
k = sum(sum(coef~=0));

% Get coefficients on exogenous variables
k_x = 0; C = [];
if size(exo,1) > 0
    k_x = size(exo,2);
    C = coef(:,end-k_x+1:end);
end;
% Get coefficients for lagged variables.
A = reshape(coef(:,1:end-k_x),n,n,p);

% Companion form VAR(p) model.
n = size(A,1);
if p == 0
    F = [];
elseif p == 1
    F = A;
elseif p > 1
    F = [reshape(A,n,n*p); [eye(n*(p-1)) zeros(n*(p-1),n)]];
end

% Results.
results.C       = C;
results.A       = A;
results.F       = F;
results.fitted  = X*coef';
results.resid   = E;
results.Sig     = Sig;
% Log likelihood
results.LogL    = -T*n/2*(1+log(2*pi)) - T/2*log(det(Sig));
% Information criteria
results.aic     = log(det(Sig)) + 2*k/T;
results.sic     = log(det(Sig)) + k*log(T)/T;
results.hqc     = log(det(Sig)) + 2*k*log(log(T))/T;
% Data
results.Y       = Y;
results.X       = X;
% Parameters
results.k       = size(X,2);
results.kn      = k;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%