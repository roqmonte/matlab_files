function irf = IRFar_build(rho,hmax,label)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 28/Oct/2017
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Computes IRF based on a set of coefficitnes.
% Inputs:
%   rho       : Vector with coefficient from AR part of the model.
%   hmax      : Horizon for IRF.
%   label     : Column vector with label dep variable.
%
% Outputs:
%   irf       : IRF for AR(p) model to hmax.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Label for chart.
if isempty(label)
    label = 'Dep. variables';
end  

% Building IRF using Companion form of the model
p = size(rho,1);

A = rho';
if p == 1
    F = A;
elseif p > 1
    F = [A; [eye(p-1) zeros(p-1,1)]];
end

% Storing IRF up to h-horizon.
J   = [eye(1) zeros(1,p-1)];
irf = zeros(hmax,1);
for i = 0:35
    irf(i+1,1) = J*(F)^(i)*J';
end

% IRF plot.
plot((0:12),irf(1:13),'-b','LineWidth',1);
% Labels.
title('Impulse response function','FontSize',10);
legend(label,'Location','northwest')
xlim([0 12]);
box off
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%