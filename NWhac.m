function result = NWhac(y,qn)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Computes Newey & West HAC matrix.
% Inputs:
%   y   : Vector with data.
%   qn  : Truncation lag.
%
% Outputs:
%   results : Newey and West HAC covariance estimator.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Deviation to mean.
T  = size(y,1); 
dy = y - repmat(mean(y),T,1);
% Case on which no truncation lag is specified.
if nargin == 1
   qn = floor(T^(1/4));  
end

% Variance at t.
G0 = dy'*dy/T;

% Adding auto-covariance.
for j = 1:qn-1
   % Covariance.
   gamma = (dy(j+1:T,:)'*dy(1:T-j,:))./(T-1);
   % Summation.
   G0 = G0 + (gamma+gamma').*(1-abs(j/qn));
end
% HAC matrix according to Newey and West.
result = G0;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%