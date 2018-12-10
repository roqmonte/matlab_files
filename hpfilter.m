function [s,desvabs] = hpfilter(y,w,print)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Hodrick-Prescott filter. Recomended values for w based on the
% frequency of y: 129,600(M), 1600(Q), % 6.25(A).
% Inputs:
%  y        : Original series.
%  w        : Smoothing parameter.
%  print    : (1) Do graph.
%
% Outputs:
%   s       : Trend HP.
%   desvabs : cyclical component HP.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
% Checking inputs.
if nargin < 2
    error('Requires at least two arguments.');
end
% Data info.
[m,n] = size (y);
if m < n
    y = y';
    m = n;
end

% Do filter
d = repmat([w -4*w ((6*w+1)/2)], m, 1);
d(1,2) = -2*w;      d(m-1,2) = -2*w;
d(1,3) = (1+w)/2;   d(m,3) = (1+w)/2;
d(2,3) = (5*w+1)/2; d(m-1,3) = (5*w+1)/2;
B = spdiags(d, -2:0, m, m);
B = B+B';
s = B\y;

% Cyclical component.
desvabs = mean(abs(y-s)./s);

% Do graph.
if nargin == 3
    if print == 1
        t = size(y,2);
        for i = 1:t
            figure(i)
            plot(s(:,i),'r');   
            grid on;   
            hold on;   
            plot(y(:,i));   
            title(['HP trend vs Series #',num2str(i)]);
        end
    end
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%