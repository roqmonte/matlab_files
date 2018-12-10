function VAR_Plot_fit(results,info,fid,vars)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 29/Jan/2017
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Plots data, fitted values and residuals from VAR model
% estimation.
% Input:
%   results:
%   -.Y           : Dependent variables of the system (lhs variables).
%   -.fitted      : Fitted values.
%   -.resid       : Regression residuals.
%   -.Sig         : Error covariance matrix.
%   info:
%   -.dates_xTick : Xtick for dates.
%   -.dates_label : Labels for dates.
%   -.names       : Labels with variable names.
%   -.widths      : vector with line widths.
%   -.fsizes      : vector with font sizes.
%   fid           : Figure number for plots.
%   vars          : Data selection for charts
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Geting inputs
if exist('fid','var') == 0
    fid = 1;
end;
grey            = [0.3,0.3,0.3];
std_error       = sqrt(diag(results.Sig));
data            = results.Y;
fitted          = results.fitted;
resid           = results.resid;
line_width      = info.widths(1);
line_width_alt  = info.widths(2);
fsize           = info.fsizes(1);
fsize_alt       = info.fsizes(2);
names           = info.names;
xTick           = info.dates_xTick;
xTickLabel      = info.dates_label;

% Data selection for charts
if exist('vars','var') == 1
    std_error = std_error(vars,1);
    data      = data(:,vars);
    fitted    = fitted(:,vars);
    resid     = resid(:,vars);
    names     = info.names(:,vars);
end
[nobs,nvars]    = size(data);

% Number of variables and graph setup.
aux = size(data,2);
if aux <= 3
    k1 = 1; k2 = aux;
elseif aux <= 4
    k1 = 2; k2 = 2;
elseif aux > 4 && aux <= 6
    k1 = 3; k2 = 2;
elseif aux > 6 && aux <= 9
    k1 = 3; k2 = 3;
elseif aux > 9 && aux <= 16
    k1 = 4; k2 = 4;
elseif aux > 16 && aux <= 24
    k1 = 4; k2 = 6;
elseif aux > 24 && aux <= 30
    k1 = 5; k2 = 6;
elseif aux > 30
    error('Max number of variables reached.');
end

% Plot data against fitted values
figure(fid);
k = 1;
for j = 1:nvars
    subplot(k2,k1,k);
    plot([NaN(xTick(end)-nobs,1); data(:,j)],'-b','LineWidth',line_width);
    hold on;
    plot([NaN(xTick(end)-nobs,1); fitted(:,j)],'--r','LineWidth',line_width);
    title(names(j),'FontSize',fsize);
    set(gca,'FontSize',fsize_alt);
    xlim([xTick(1) xTick(end)]);
    set(gca,'XTick',xTick);
    set(gca,'XTickLabel',xTickLabel);
    k = k + 1;
end;
legend1 = legend('Data','Fitted values');
set(legend1,'FontSize',fsize_alt,'Orientation','horizontal',...
'Position',[0.327835 0.00079499 0.3672619 0.05176190476],'Box', 'off');
clear k j fid1;

% Plot Regression residuals
figure(fid+1);
k = 1;
for j = 1:nvars
    subplot(k2,k1,k);
    h(1) = plot([zeros(xTick(end)-nobs,1); zeros(nobs,1)],'-','Color',grey,'LineWidth',line_width_alt);
    hold on;
    h(2) = plot([NaN(xTick(end)-nobs,1); resid(:,j)],'-b','LineWidth',line_width);
    h(3) = plot(repmat(std_error(j),nobs+(xTick(end)-nobs),1),':','Color',grey,'LineWidth',line_width_alt);
    h(4) = plot(repmat(-std_error(j),nobs+(xTick(end)-nobs),1),':','Color',grey,'LineWidth',line_width_alt);
    title(names(j),'FontSize',fsize);
    set(gca,'FontSize',fsize_alt);
    xlim([xTick(1) xTick(end)]);
    set(gca,'XTick',xTick);
    set(gca,'XTickLabel',xTickLabel);
    k = k + 1;
end;
legend1 = legend(h([2 3]),'Residuals','+/-1 Standard Deviation');
set(legend1,'FontSize',fsize_alt,'Orientation','horizontal',...
'Position',[0.327835 0.00079499 0.3672619 0.05176190476],'Box', 'off');
clear aux1 aux2 k j fid2 T_aux aux1 aux2 fid j k T_aux;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%