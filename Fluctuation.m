function results = Fluctuation(DeltaL,dates,ws,hac,alpha,print)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Computes fluctuation test of Giacomini & Rossi (2010).
% Inputs:
%   DeltaL : Square diff of forecast errors.
%   dates  : Vector with info for dates: [year,month,freq].
%            Freq:(1) monthly data; (2) quaterly data.
%   ws     : Rolling window size for estimation (Default: 25%T).
%   hac    : Variance estimation: (0) OLS ; (1) HAC matrix (Default).
%   alpha  : Significance level 10%/5% (Default).
%   print  : (1) Do graph.
%
% Outputs:
%   results :  
%   -.F     : Fluctuation test.
%   -.cvlow : Lower bound for the test.
%   -.cvup  : High bound for the test.
%   -.tests : Date for the one time reversal test.
%   -.dates : Dates consistent for graphing. 
%
% Index:
% 1. Test.
% 2. Functions.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
% 1. Test.
% Checking if variables are defined, setting default variables if they are
% not.
if exist('dates','var') == 0 || isempty(dates)
    dates = [1900,1,1]; 
end
if exist('ws','var') == 0 || isempty(ws)
    ws = round(0.25*size(DeltaL,1));
end
if exist('hac','var') == 0 || isempty(hac)
    hac = 1;
end
if exist('alpha','var') == 0 || isempty(alpha)
    alpha = 0.05;
end

% Filter of the data.
y1 = filter(ones(1,ws)/ws,1,DeltaL); 
y  = y1(ws:end,:); 
% Sample size and kernel truncation parameter.
T  = length(y); 
qn = floor(T^(1/4)); 

% HAC estimation for the variance.
if hac == 1
    sigma = sqrt(NWhac(DeltaL,qn));
else
    sigma = sqrt(cov(DeltaL));
end 

% Statistic for the Fluctuation test.
F = sqrt(ws)*(y./sigma);

% Critical values for the fluctuation test.
%       ws/T  0.05   0.10
cvtable=[0.1, 3.393, 3.170;
         0.2, 3.179, 2.948;
         0.3, 3.012, 2.766;
         0.4, 2.890, 2.626;
         0.5, 2.779, 2.500;
         0.6, 2.634, 2.356;
         0.7, 2.560, 2.252;
         0.8, 2.433, 2.130;
         0.9, 2.248, 1.950];
% Finding the critical value.
if alpha==0.05
    j = 2; 
elseif alpha==0.10
    j = 3;
else
    error('Wrong input for alpha.');
end
mu = (ceil((ws/T)*10))/10; 
i  = mu*10; 
cv = cvtable(i,j);

% Graph fluctuation test.
if nargin == 6
    if print == 1
        dta1 = flucgraph(F,cv,dates,size(DeltaL,1));
    end    
end

% Results
results.F     = F;
results.cvup  = repmat(cv,size(F,1),1);
results.cvlow = -1*repmat(cv,size(F,1),1);
results.dates = dta1;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
% 2. Functions.
% Function graphs the fluctiation test.
function dt1 = flucgraph(series,cv,dates,T)
% Inputs:
%   series: Fluctuation test.
%   cv    : Critical value test.
%   dates : Dates for the graph.
%   T     : Original sample size

% Making dates
T2 = size(series,1);
if dates(end) == 1
    freq = 'm';
else
    freq = 'q';
end
[dt1,dt2,f] = calendar(dates(1),dates(2),T,freq);
% Fixing time lenght.
dt1 = dt1(end-T2+1:end);
dt2 = dt2(end-T2+1:end);
% Confidence intervals
ub = repmat(cv,size(dt1,1),1);
lb = -1*ub;

% Figure
figure(1);
% figure(grid)
plot(dt1,zeros(size(series,1),1), '--k', 'LineWidth',2)
hold on
plot(dt1,series,'LineWidth',2)
hold on
plot(dt1,lb,'-r','LineWidth',2)
hold on
plot(dt1,ub,'-r', 'LineWidth',2)
% Format.
xlim([dt1(1) dt1(end)]);
ylim([(min([series;lb]) - 0.2) (max([series;ub]) + 0.2)])
set(gca,'XTick',dt1(1:f:T2));
set(gca,'XTickLabel',dt2(1:f:T2));
title('Fluctuation test','FontSize',14);

% Labels
xlabel('Time', 'FontSize',12);
ylabel('rMSFE','FontSize',12)
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%