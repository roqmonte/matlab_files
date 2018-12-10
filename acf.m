 function [results,mat] = acf(x,print)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Computes autocorrelation and partial auto-correlation fns.

% Inputs:
%  x     : Time-series vector.
%  print : (1) do table and graph.
%
% Outputs:
% results: Table Autocorrelation and Partial Autocorrelation functions 
%          (AC function, CI, Qstat and Pvalues).
% mat    : Matrix with results Autocorrelation and Partial Autocorrelation functions 
%          (AC function, CI, Qstat and Pvalues).
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Error checking on inputs
if (size(x,2) > 1)
    error('acf cannot handle a matrix -- only vectors');
end

% Sample information
nobs = size(x,1);

% Truncation order for HAC matrix.
k = round(nobs/4)+10;
% De mean of the variable
xc = x - mean(x);
v0 = (xc'*xc)/nobs; 

% Computing the autocorrelations.
for i=1:k
    yt  = xc(1:nobs-i,1);
    ytk = xc(1+i:nobs,1);
    ck  = (yt'*ytk)/nobs;
    acf(i,1) = ck/v0;   
end
% Significance of xk (95% confidence intervals)
lb = -2/sqrt(nobs)*ones(k,1);
ub = 2/sqrt(nobs)*ones(k,1);  
% Diagnostics (Ljung & Box, 1978) Q-Statistic
qstat=zeros(k+1,1);
for j=2:k+1
    qstat(j,1) = (acf(j-1,1)^2)/(nobs-(j-1))+ qstat(j-1,1);
    prob(j,1)=  1-chis_prb(abs((nobs*(nobs+2))*qstat(j,1)),j);
end
% Formating results
qstat = trimr((nobs*(nobs+2))*qstat,1,0);
prob  = trimr(prob,1,0);

% Partial correlation function
[~,data] = LagN(x,k+1);
pac   = zeros(k,1);
ytemp = data(:,1);
for order = 2:k+1
    [Q,R] = qr([ones(nobs-order+1,1) data(order:end,2:order)],0); 
    b     = R\(Q'*ytemp(order:end));
    pac(order-1,1) = b(end);
end

% Results.
res1 = [acf pac lb ub qstat prob];
% Labels
labels = {'Lag' 'AC fun' 'PAAC fun' 'Lb' 'Ub' 'Qstat' 'Pvalue'};
% Autocor function.
aux1 = (1:k)';
aux2 = [aux1 res1];
res2 = num2cell(aux2(1:12,:));
results = [labels; res2];
mat  = res1;
% Do graph
if nargin == 2
    if print == 1
        % Print table        
        display(results);
        
        % Ploting reults
        subplot(2,1,1)
        figure(1)
        bar(res1(:,1));
        hold on
        plot(res1(:,3), '--r');
        hold on
        plot(res1(:,4), '--r');
        xlim([1 k]);
        ylim([-1 1]);
        % Labels
        title('Autocorrelation funtion','FontSize',11);

        subplot(2,1,2)
        bar(res1(:,2));
        hold on
        plot(res1(:,3), '--r');
        hold on
        plot(res1(:,4), '--r');
        xlim([1 k]);
        ylim([-1 1]);
        % Labels
        title('Partial-Autocorrelation funtion','FontSize',11);
    end
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%