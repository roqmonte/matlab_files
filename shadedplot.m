function [ha,hb,hc] = shadedplot(x,y1,y2,varargin)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Draws two lines on a plot and shades area between them. 
% Vector x contains x-axis values, and y1:y2 contain y-axis values. 
% The arguments areacolor and linecolor allow user to set color of shaded 
% area and boundary lines. 
% Inputs:
% 	x        : X-axis values.
%   y1       : Upper bound.
%   y2       : Lower bound.
%   varargin : Color setup, face color and line color.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Plot the shaded area
y = [y1; (y2-y1)]'; 
ha = area(x, y);
% This makes the bottom area invisible
set(ha(1), 'FaceColor', 'none') 
set(ha, 'LineStyle', 'none')

% Plot the line edges
switch length(varargin)
    case 2
        hold on
        hb = plot(x, y1, 'LineWidth', 1);
        hc = plot(x, y2, 'LineWidth', 1);
        hold off
end

% Set the line and area colors if they are specified
switch length(varargin)
    case 0
    case 1
        set(ha(2), 'FaceColor', varargin{1})
    case 2
        set(ha(2), 'FaceColor', varargin{1})
        set(hb, 'Color', varargin{2})
        set(hc, 'Color', varargin{2})
    otherwise
end

% Put grid on top of the colored area
set(gca, 'Layer', 'top')
grid off
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%