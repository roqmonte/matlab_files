function [mat1,mat2] = LagN(y,n)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Generates matrix with lags values of y. Note: output contains 
% vector y in the first column.
% Inputs:
%   y    : Vector of data
%   n    : Number of Lags
%
% Outputs:
%   mat1 : Matrix with Lags of y, excluding first n observations.
%   mat2 : Matrix with Lags of y, including first n observations (NaN).
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
% Case 1: y is a vector.
colsy = size(y,2);

if size(y,2) == 1
    z = y;
    for i1 =1:n
        zt = NaN(i1,colsy);
        zz = [zt; trimr(y,0,i1)];
        z  = [z zz];
    end
    % Computing mat1
    mat1 = z(n+1:end,:);
    % Computing mat2
    mat2 = z;
    
% Case 2: y is a matrix.
else
    zaux = [];
    for i1 =1:n
        zt = NaN(i1,colsy);
        zz = [zt; trimr(y,0,i1)];
        zaux = [zaux zz];
    end
    zaux = [y zaux];
    % Computing mat1
    mat1 = zaux(n+1:end,:);
    % Computing mat2
    mat2 = zaux;
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%