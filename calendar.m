function [xTick,xTickLabel,dates_all] = calendar(fyds,fmds,nds,f)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Generats dates to plot figures using matlab functions.
% Inputs:
%   fyds : Initial year.
%   fmds : Initial month/quarter.
%   nds  : Number of observations.
%   f    : Frequency; (m) for monthly and (q) for quarterly data.
%
% Outputs:
%   xTick      : Vector with tickers for dates
%   xTickLabel : Vector with labels for dates.
%   dates_all  : Code for all dates
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
% Type of data
% Monthly data.
if f == 'm'
    % Checking initial month
    if fmds > 13
        error('Check initial month');
    end
    tid = 12;
    % Dates for monthly data.
    i = 1;
    for i0 = 0.01:0.01:0.12
       aux(i,1) = i0;
       i = i + 1;
    end
    % Initial year.
    year = fyds;
    % Dates.
    final_dates = zeros(nds,1);
    i = 1;
    for i0 = 1:nds + (fmds-1)
        final_dates(i0,1) = year + aux(i,1);
        i = i + 1;
        if i == 13
           year = year + 1;
           i = 1;
        end
    end
end

% Quaterly data.
if f == 'q'
    if fmds > 4
        error('Check initial quarter');
    end
    tid = 4; 
    % Dates for monthly data.
    i = 1;
    for i0 = 0.01:0.01:0.04
        aux(i,1) = i0;
        i = i + 1;
    end
    % Initial year.
    year = fyds;
    % Dates.
    final_dates = zeros(nds,1);
    i = 1;
    for i0 = 1:nds + (fmds-1)
        final_dates(i0,1) = year + aux(i,1);
        i = i + 1;
        if i == 5
           year = year + 1;
           i = 1;
        end
    end
end

% Final results
calds = [round(final_dates(fmds:end,1)) round((final_dates(fmds:end,1)-round(final_dates(fmds:end,1)))*100)];
time1 = calds(:,1) + (calds(:,2) - ones(size(calds,1),1))/tid;
time2 = final_dates(fmds:end,1);
aa    = 0.25;

% Building dates
ndates = time2;
Tf     = size(time1,1);
label_r = (1:round(Tf*aa):Tf)';

% Adding last label 
if label_r(end) < Tf
    label_r = [label_r; Tf];
end

% Setting labels for dates
new_label = [];
new_dates = ndates(label_r,1);
for i0 = 1:size(new_dates,1)
    % Year an dmonth id.
    year = num2str(round(new_dates(i0)));
    month_id = (new_dates(i0) - round(new_dates(i0)));    
    if f == 'm'
        % Check whether months are identified correctly; otherwise correct
        if round(month_id*100) == 1
            month = 'Jan';
        elseif round(month_id*100) == 2
            month = 'Feb';
        elseif round(month_id*100) == 3
            month = 'Mar';
        elseif round(month_id*100) == 4
            month = 'Apr';
        elseif round(month_id*100) == 5
            month = 'May';
        elseif round(month_id*100) == 6
            month = 'Jun';
        elseif round(month_id*100) == 7
            month = 'Jul';
        elseif round(month_id*100) == 8
            month = 'Aug';
        elseif round(month_id*100) == 9
             month = 'Sep';
        elseif round(month_id*100) == 10
            month = 'Oct';
        elseif round(month_id*100) == 11
            month = 'Nov';
        elseif round(month_id*100) == 12
            month = 'Dec';
        end
    elseif f == 'q'
        % Check whether quaters are identified correctly; otherwise correct
        if round(month_id*100) == 1
            month = 'Q1';
        elseif round(month_id*100) == 2
            month = 'Q2';
        elseif round(month_id*100) == 3
            month = 'Q3';
        elseif round(month_id*100) == 4
            month = 'Q4';
        end
    end
    new_label = [new_label; [year month]];
end
% new labels Lables
xTick     = label_r;
xTickLabel= new_label;
dates_all = final_dates;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    