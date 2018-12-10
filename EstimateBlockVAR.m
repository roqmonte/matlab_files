function results = EstimateBlockVAR(data,p,exo)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 12/Oct/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Estimate VAR model with block exogeneity (y1_t wrt y2_t).
% The original regression model is given by the two-dimensional system
%   (1) y1_t = C_1 x_t + A_11 y1_t-1 + ... + A_1p y1_t-p + e_1t,
%   (2) y2_t = C_2 x_t + B_11 y1_t-1 + ... + B_1p y1_t-p
%                      + B_21 y2_t-1 + ... + B_2p y2_t-p + e_2t,
% Where cov([e_1t; e_2t]) = [Sig_11 Sig_12; Sig_21 Sig_22].
%
% A reparameterized regression model is given by
%   (1a) y1_t = C_1 z_t + A_11 y1_t-1 + ... + A_1p y_1t-p + e_1t,
%   (2a) y2_t = D z_t + D_0 y1_t
%                     + D_11 y1_t-1 + ... + D_1p y1_t-p
%                     + D_21 y2_t-1 + ... + D_2p y2_t-p + v_2t,
% where cov(v_2t) = H.
% The parameters of the original regression model (full information
% likelihood estimates) are recovered from the estimates of the
% reparameterized model as follows:
%  (i)   Sig_21 = D_0 * Sig_11
%  (ii)  C_2    = D   + Sig_21 * inv(Sig_11) * C_1
%  (iii) B_1    = D_1 + Sig_21 * inv(Sig_11) * A_1
%  (iv)  B_2    = D_2
%  (v)   Sig_22 = H   + Sig_21 * inv(Sig_11) * Sig_12
% Inputs:
%   data:
%   -.exo       : Data for y1_t block.
%   -.endo      : Data for y2_t block
%   p           : Lag order p.
%   exo         : Matrix, exogenous variables (optional).
%
% Outputs:
%   results:
%   -.A_1       : Lag matrices for exo block.
%   -.B_1       : Lag mats. endo block (exo vars; n2 x n1 x p).
%   -.B_2       : Lag mats. endo block (end vars; n2 x n2 x p).
%   -.C_1       : Coefficients exogenous variables for y1t.
%   -.C_2       : Coefficients exogenous variables for y2t.
%   -.F         : Companion form.
%   -.fitted    : Fitted values.
%   -.resid     : Regression residuals.
%   -.Sig       : Matrix, error covariance matrix.
%   -.LogL      : Log-likelihood.
%   -.aic       : Akaike information criterion.
%   -.hqc       : Hannan-Quinn information criterion.
%   -.sic       : Schwarz information criterion.
%   -.Y         : Dependent variables of the system (lhs variables).
%   -.X         : Independent variables of the system (rhs variables).
%   -.kn        : Number of parameters of the VAR model, block1, block2.
%
% Index:
% 1. Initial Setup of the model.
% 2. Model estimation.
% 3. Results.
% 4. Additional Results.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 1. Initial Setup of the model.
% Getting info for estimation.
if exist('exo','var') == 0
    exo = [];
end

% Provide inputs
exo_data  = data.exo;
endo_data = data.endo;
data      = [exo_data endo_data];
n_1       = size(exo_data,2);
n_2       = size(endo_data,2);
n         = n_1 + n_2;
T         = length(exo_data)-p;

% Setup endogenous variables for each block.
Y_1 = exo_data(1+p:end,:);
Y_2 = endo_data(1+p:end,:);

% Lags for exo block and endo block (respectively).
X_1 = [];
X_2 = exo_data(1+p:end,:);
i = 1;
while i <= p
    X_1 = [X_1 exo_data(1+p-i:end-i,:)];
    X_2 = [X_2 exo_data(1+p-i:end-i,:)];
    i = i + 1;
end;
i = 1;
while i <= p
    X_2 = [X_2 endo_data(1+p-i:end-i,:)];
    i = i + 1;
end;

% Adding exogenous variables
X_1 = [X_1 exo(1+p:end,:)];
X_2 = [X_2 exo(1+p:end,:)];
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 2. Model estimation.
% OLS estimates equation by equation
coef_1 = (X_1'*X_1\X_1'*Y_1)';
coef_2 = (X_2'*X_2\X_2'*Y_2)';

% Computing residuals for each block.
E_1 = (Y_1-X_1*coef_1');
V_2 = (Y_2-X_2*coef_2');

% Compute residual covariances for each block.
Sig_11 = E_1'*E_1/T;
H = V_2'*V_2/T;

% Get coefficients on exogenous variables
k_x = 0;
C_1 = []; D = [];
if size(exo,1) > 0
    k_x    = size(exo,2);
    C_1    = coef_1(:,end-k_x+1:end);
    D      = coef_2(:,end-k_x+1:end);
    coef_1 = coef_1(:,1:end-k_x);
    coef_2 = coef_2(:,1:end-k_x);
end;


% Get coefficients on contemporaneous endogenous variables
D_0    = coef_2(:,1:n_1);
coef_2 = coef_2(:,n_1+1:end);

% Get coefficients on lagged endogenous variables
A_1 = reshape(coef_1,n_1,n_1,p);
D_1 = reshape(coef_2(:,1:n_1*p),n_2,n_1,p);
D_2 = reshape(coef_2(:,n_1*p+1:end),n_2,n_2,p);

% Recover parameters of original model
Sig_21 = D_0*Sig_11;
Sig_12 = Sig_21';
Sig_22 = H + Sig_21/Sig_11*Sig_12;
Sig    = [Sig_11 Sig_12; Sig_21 Sig_22];
for j = 1:p
    B_1(:,:,j) = D_1(:,:,j) + Sig_21/Sig_11*A_1(:,:,j);
    B_2(:,:,j) = D_2(:,:,j);
end;
C_2 = [];
if size(exo,1) > 0
    C_2 = D + Sig_21/Sig_11*C_1;
end;

% Case for lag lenght test only.
if p == 0
    B_1 = [];
    B_2 = [];
end;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 3. Results.
% Dependent variables
Y = [Y_1 Y_2];
% Set up independent variables
X = []; i = 1;
while i <= p
    X = [X data(1+p-i:end-i,:)];
    i = i + 1;
end;
X = [X exo(1+p:end,:)];

% Computing fitted values
coef = [];
for i = 1:p
    tmp = [A_1(:,:,i) zeros(n_1,n_2); B_1(:,:,i) B_2(:,:,i)];
    coef = [coef tmp];
end;
tmp  = [C_1; C_2];
coef = [coef tmp];

% Companion form VAR(p) model.
n_1 = size(A_1,1);
n_2 = size(B_1,1);
n = n_1 + n_2;
if p == 0
    F = [];
elseif p == 1
    F1 = [A_1(:,:,1) zeros(n_1,n_2)];
    F2 = [B_1(:,:,1) B_2(:,:,1)];
    F  = [F1; F2];
elseif p > 1
    F1 = [];
    F2 = [];
    for i = 1:p
        F1 = [F1 A_1(:,:,i) zeros(n_1,n_2)];
        F2 = [F2 B_1(:,:,i) B_2(:,:,i)];
    end;
    F3 = [F1; F2];
    F  = [F3; [eye(n*(p-1)) zeros(n*(p-1),n)]];
end;
clear i F1 F2 F3;

% Total number of parameters
k = sum(sum(coef~=0));
% Results
results.A_1     = A_1;
results.B_1     = B_1;
results.B_2     = B_2;
results.C_1     = C_1;
results.C_2     = C_2;
results.F       = F;
results.fitted  = X*coef';
results.resid   = Y - results.fitted;
results.LogL    = -T*n/2*(1+log(2*pi))-T/2*log(det(Sig));
results.aic     = log(det(Sig)) + 2*k/T;
results.sic     = log(det(Sig)) + k*log(T)/T;
results.hqc     = log(det(Sig)) + 2*k*log(log(T))/T;
results.Sig     = Sig;
% Data
results.Y       = Y;
results.X       = X;
% Parameters
results.kn      = k;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%