function [new_data,nout] = data_process(data,trans,check)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Transforms data based on transformation code. 
% Inputs:
%   data  : Data (TxN)
%   trans : Transformation code, seven options.
%           (1) Level (i.e. no transformation): x(t)
%           (2) First difference: x(t)-x(t-1)
%           (3) Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
%           (4) Natural log: ln(x)
%           (5) First difference of natural log: ln(x)-ln(x-1)
%           (6) Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
%           (7) First difference of percent change: (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
%   check : (0) no change in data; (1) fix outliers
%
% Outputs:
%   new_data : Transformed data.
%   nout     : Number of outliers.
%
% Index:
% 1. Code.
% 2. Functions.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
% 1. Code.
% Setup
T  = size(data,1);
N  = size(data,2);
yt = NaN(T,N);

% Perform transformation.
for i = 1:N
    dum = transxf(data(:,i),trans(i));
    yt(:,i) = dum;
end

% Checking for Outliers.
[data_fix,nout] = outliers(yt);
if max(nout) > 0
    disp('1. Check for outliers');
elseif max(nout) == 0
    disp('1. No outliers detected in any series');
    nout = [];
end

% Results.
if check == 1
    disp('2. Outliers were replaced');
    new_data = data_fix;
elseif check == 0
    disp('2. Outliers were not replaced');
    new_data = yt;
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%   
% 2. Functions.
% Function transforms a vector according to the transformation code.
function y = transxf(x,tcode)
% input:
%   x      : Series (Tx1).
%   tcode  : Transformation code (1-7)
% output:   
%   y      : Transformed series (as a column vector)
%
% Initial setup.
n     = size(x,1);
small = 1e-6;
y     = NaN*ones(n,1);

% Transformations.
% Level (i.e. no transformation): x(t)
if tcode == 1
    y = x;
% First difference: x(t)-x(t-1)
elseif tcode ==  2
    y(2:n) = x(2:n,1) - x(1:n-1,1);
% Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))  
  elseif tcode == 3
    y(3:n) = x(3:n) - 2*x(2:n-1) + x(1:n-2);
% Natural log: ln(x)
  elseif tcode == 4
    if min(x) < small 
        y=NaN; 
    else
        y=log(x);
    end;
% First difference of natural log: ln(x)-ln(x-1)  
  elseif tcode == 5
    if min(x) > small
        x = log(x);
        y(2:n) = x(2:n) - x(1:n-1);
    end;
% Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))  
  elseif tcode == 6
    if min(x) > small
        x = log(x);
        y(3:n) = x(3:n) - 2*x(2:n-1) + x(1:n-2);
    end;
% First difference of percent change: (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)  
  elseif tcode == 7
  y1(2:n) = (x(2:n) - x(1:n-1)) ./ x(1:n-1);
  y(3:n)  = y1(3:n) - y1(2:n-1);
end;

% Function replaces outliers with NaN values. Outlier definition: data 
% point x of a series X(:,i) if abs(x-median)>10*interquartile_range.
function [Y,n] = outliers(X)
% Input:
%   x   : data (TxN).
% Output:   
%   Y   : dataset with outliers replaced with NaN 
%   n   : number of outliers found in each series
%
% Computing median of each series
median_X     = nanmedian(X,1);
median_X_mat = repmat(median_X,size(X,1),1);
% Computing quartiles and interquartile range (IQR).
Q   = prctile(X,[25, 50, 75],1);
IQR = Q(3,:)-Q(1,:);
% Repeat IQR of each series over all data points in the series
IQR_mat = repmat(IQR,size(X,1),1);

% Determine outliers 
Z = abs(X - median_X_mat);
outlier = Z > (10*IQR_mat);

% Replace outliers with NaN and num of outliers.
Y = X;
Y(outlier) = NaN;
n = sum(outlier,1);
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%   