function results = MSarp_flex(y,x,x1,x2,x3,Sg2,Pr,print,labels,dates)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 29/Nov/2018
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Estimation of a MS-AR(p) model with different specification
% across regimes and with transition probability given by:
% P, where P is given by:
%
%     |  P11  P21 ... Pj1  |
%     |  P12  P22 ... Pj2  |
% P = |   :    :  ...  :   |
%     |   :    :  ...  :   |
%     |  P1j  P2j ... Pjj  |
%
% Where Pkj is the probability of moving from regime k to j.
%
% Imputs:
%   y       : Endogenous variable.
%   x       : Independet variables.
%   x1      : Independet variables for regime 1.
%   x2      : Independet variables for regime 2.
%   x3      : Independet variables for regime 3.
%   Sg2     : Variance specification, (1) same across regimes; (2) regime dependent.
%   Pr      : Set restrictions for P matrix (unit entries).
%   print   : (1) Do charts and print results on screen.
%   labels  : Column vector with labels for dep. and indep variables.
%   dates   : Vector with info for dates: [year,month,freq].
%             Freq:(1) monthly data; (2) quaterly data.
%
% Outputs:
%   results   :
%   -.logf     : Loglikelihood.
%   -.Fil      : Filter probability.
%   -.Smt      : Smoothed probability .
%   -.H        : Inverse Hessian matrix (var/cov matrix).
%   -.niter    : Number of iterations.
%   -.dta      : Data estimation.
%   -.yhat     : Fit of the model.
%   -.yhat2    : Fit of each regime.
%   -.uhat     : in-sample residuals model.
%   -.uhat2    : in-sample residuals from each regime.
%   -.thetaf   : Theta hat.
%   -.Param    : Parameters of the model.
%   -.Pval     : Asymtotic Pvalues.
%   -.Sg2      : variance.
%   -.SSR      : Sum square resids.
%   -.R2       : R-square.
%   -.R2adj    : Adjusted R-square.
%   -.AIC      : AIC.
%   -.HQC      : HQC.
%   -.BIC      : BIC.
%   -.T        : Sample size.
%   -.k        : Number of parameters.
%   -.table    : Table with print results.
%
% Index:
% 1. Initial Setup.
% 2. Estimation of MSAR models.
% 3. Additional results.
% 4. Function.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%  
% 1. Initial Setup.
% Checking labels.
if exist('labels','var') == 0 || isempty(labels)
    labels = [];
end
% Checking print.
if exist('print','var') == 0
    print = 0;
end
% Checking cong model
if exist('x1','var') == 0 || exist('x2','var') == 0
    error('No data for one of the Regime');
end

% Checking inputs number of regimes
St = (1-isempty(x1))+(1-isempty(x2))+(1-isempty(x3));
% Checking restriction P matrix
if exist('Pr','var') == 0 || isempty(Pr)
    if St == 2
       Pr = ones(2,2); 
    elseif St == 3
       Pr = ones(3,3);
    end
end
% No matrix Pr.
if size(Pr,1) ~= St || size(Pr,2) ~= St
    error('Check dim of Pr matrix.');
end

% Checking labels.
if exist('labels','var') == 0 || isempty(labels)
    % Label for dep. variable
    lab_y = {'dep var'};
    % Labels for x's variables.
    lab_x = {};
    for i0 = 1:size(x,2)
        aux  = strcat('xvar ',num2str(i0));
        lab_x = [lab_x; aux];
    end
    % Labels for x1's variables.
    lab_x1 = {};
    for i0 = 1:size(x1,2)
        aux    = strcat('Reg1 xvar',num2str(i0));
        lab_x1 = [lab_x1; aux];
    end
    % Labels for x2's variables.
    lab_x2 = {};
    for i0 = 1:size(x2,2)
        aux    = strcat('Reg2 xvar',num2str(i0));
        lab_x2 = [lab_x2; aux];
    end
    % Labels for x3's variables.
    if St == 3
        lab_x3 = {};
        for i0 = 1:size(x3,2)
            aux    = strcat('Reg3 xvar',num2str(i0));
            lab_x3 = [lab_x3; aux];
        end 
    else
        lab_x3 = {};
    end
    labels = [lab_y; lab_x; lab_x1; lab_x2; lab_x3];
end
clear aux i0 lab_y lab_x lab_x1 lab_x2 lab_x3;

% Initial values.
if isempty(x) == 1
    B0 = [];
else
    B0 = (x'*x)\(x'*y);
end
% Initial betas for Reg. 1 & 2.
B1 = (x1'*x1)\(x1'*y);
B2 = (x2'*x2)\(x2'*y);   
sigma(1,1) = sqrt((y - x1*B1)'*(y - x1*B1)/(size(y,1)-size(B1,1)));
sigma(2,1) = sqrt((y - x2*B2)'*(y - x2*B2)/(size(y,1)-size(B2,1)));
% Initial betas for Reg. 3.
B3 = [];
if St == 3
    B3 = (x3'*x3)\(x3'*y);
    sigma(3,1) = sqrt((y - x3*B3)'*(y - x3*B3)/(size(y,1)-size(B3,1)));
end
% Iniial betas.
betas = [B0; B1; B2; B3];

% Condition for Sg2
if Sg2 == 1
    sigma = mean(sigma);
end
l = sigma;
% Selecting the number of states and defining initial values for the prob.
if St == 2
    p_l = ones(2,1);
elseif St == 3
    p_l = [ones(sum(Pr(:,1))-1,1); ones(sum(Pr(:,2))-1,1); ones(sum(Pr(:,3))-1,1)];   
else
    error('Wrong number of regimes.');
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 2. Estimation of MSAR model.
% Infor for the opt. routine.
theta = [betas;l;p_l];
T = size(y,1);

% Hessian Matrix.
H = fdhess('MSarp_evalf_flex',theta,y,x,x1,x2,x3,Pr);
H0inv = H\eye(size(H));
% Checking the consistency of the Hessian Matrix.
if max(max(isnan(H0inv))) == 1 || max(max(isinf(H0inv))) == 1 || isreal(H0inv) == 0 || isreal(H) == 0
    H = eye(size(theta,1));
else
    [V,Diag] = eig(H0inv);
    Diag     = abs(Diag);
    H        = V*Diag*V';
end;
% Numerical Optimización.
tol = 1e-08; n_ite = 2000; 
[~,thetaf,~,H,~,~,niter] = csminwel('MSarp_evalf_flex',theta,H,[],tol,n_ite,y,x,x1,x2,x3,Pr);

% Filter and Smoothed probabilities.
[lf,e0,sm,param,yhat,yhat2,uhat] = MSarp_evalf_flex(thetaf,y,x,x1,x2,x3,Pr);
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 3. Results.
% Pvalues of parameters.
k   = size(thetaf,1);
nb  = size(betas,1);
vb  = diag(H(1:nb,1:nb));
bms = thetaf(1:nb);
tts = bms ./ vb.^(1/2);

% Results.
results.logf     = -1*lf;
results.Fil      = e0;
results.Smt      = sm;
results.H        = H;
results.niter    = niter;

% Estimation results.
results.dta    = [y x x1 x2 x3];
results.yhat   = yhat;
results.yhat2  = yhat2;
results.uhat   = uhat;
results.uhat2  = repmat(y,1,St) - results.yhat2;
results.thetaf = thetaf;
results.Param  = param;

% Pvalues.
Pvalb = 2*(1-tcdf(abs(tts),T - size(H,1)));
for i0 = 1:length(Pvalb)
   if Pvalb(i0,1) < 0.001
        Pvalb(i0,1) = 0;
   end
end
results.Param.Pvalb = Pvalb;

% More results
results.SSR    = uhat'*uhat;
results.Sg2    = results.Param.Sg2;
results.R2     = 1 - results.SSR / ((y-mean(y))'*(y-mean(y)));
results.R2adj  = 1 - (1 - results.R2)*(T - 1)/(T - k);
results.AIC = -(2/T)*results.logf + 2*(k/T);
results.HQC = -(2/T)*results.logf + 2*(k/T)*log(log(T));
results.BIC = -(2/T)*results.logf + (k/T)*log(T);
results.T   = T;
results.k   = k;
% Table.
results.table = print_res(results,labels,St,print);

% Do charts
if print == 1       
    if exist('dates','var') == 0 || isempty(dates)
        dates = [1900,1,1];
    end
    print_charts(results,dates,labels,St);
end    
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 4. Function.
% Funtion do table.
function Estimation_results = print_res(results,labels_all,St,print)

% Building parameters.
Bx = results.Param.Bx; nx  = size(Bx,1);
B1 = results.Param.B1; nx1 = size(B1,1);
B2 = results.Param.B2; nx2 = size(B2,1);
Pvx  = results.Param.Pvalb(1:nx);
Pvx1 = results.Param.Pvalb(nx+1:nx+nx1);
Pvx2 = results.Param.Pvalb(nx+nx1+1:nx+nx1+nx2);

% Labels
ylab   = labels_all(1);
labtmp = labels_all(2:end);
lab_x  = labtmp(1:nx);
lab_x1 = labtmp(nx+1:nx+nx1);
lab_x2 = labtmp(nx+nx1+1:nx+nx1+nx2);

% Regime 1
reg1   = [[lab_x; lab_x1] num2cell([[Bx; B1] [Pvx; Pvx1]])];
labels = {'Reg. 1' 'Param' 'Pvalue'};
temp_1 = ['Sg2' num2cell(results.Sg2(1)) ' '];
Reg_1  = [labels; reg1; temp_1];
clear labels temp_1 reg1;
% Regime 2
reg2   = [[lab_x; lab_x2] num2cell([[Bx; B2] [Pvx; Pvx2]])];
labels = {'Reg. 2' 'Param' 'Pvalue'};
temp_1 = ['Sg2' num2cell(results.Sg2(2)) ' '];
Reg_2  = [labels; reg2; temp_1];
clear labels temp_1 reg2;
% Regime 3
if St == 2
    Reg_3  = [];
elseif St == 3
    B3     = results.Param.B3; nx3 = size(B3,1);
    Pvx3   = results.Param.Pvalb(nx+nx1+nx2+1:nx+nx1+nx2+nx3);
    lab_x3 = labtmp(nx+nx1+nx2+1:nx+nx1+nx2+nx3);
    reg3   = [[lab_x; lab_x3] num2cell([[Bx; B3] [Pvx; Pvx3]])];
    labels = {'Reg. 3' 'Param' 'Pvalue'};
    temp_1 = ['Sg2' num2cell(results.Sg2(3)) ' '];
    Reg_3  = [labels; reg3; temp_1];
    clear labels temp_1 reg3;
end

% Formating table
nmax  = max([size(Reg_1,1) size(Reg_2,1) size(Reg_3,1)]);
Part_1= repmat({'' '' '' '' '' '' '' '' ''},nmax,1);
Part_1(1:size(Reg_1,1),1:3) = Reg_1;
Part_1(1:size(Reg_2,1),4:6) = Reg_2;
if St == 3
    Part_1(1:size(Reg_3,1),7:9) = Reg_3;
end

% Second part
temp_1 = {'AIC'; 'HQC';  'BIC'};
temp_2 = num2cell([results.AIC;results.HQC;results.BIC]);
temp_3 = [{'SSR'; 'R2'; 'R2-adj'} num2cell([results.SSR;results.R2;results.R2adj])];
temp_4 = [{'T'} num2cell(results.T); {'k'} num2cell(results.k); {'Iterations'} num2cell(results.niter)];
temp_5 = [{'#Regime'} num2cell(St) {'logf'} num2cell(results.logf) {' '} {' '}];
stats  = [[{'' '' '' '' '' ''}; temp_1 temp_2 temp_3 temp_4; temp_5] repmat({'' '' ''},5,1)];
clear temp_1 temp_2 temp_3 temp_4 temp_5;
% Transition matrix
aux = results.Param.P;
for i0 = 1:St
    for i1 = 1:St
        if aux(i0,i1) < 0.001
            aux(i0,i1) = 0;
        end
    end 
end
results.Param.P = aux;
if St == 2
    Pmat  = [[{'P'}  num2cell((1:St)) ;num2cell([(1:St)' results.Param.P])] repmat({'' '' ''},3,2)];
elseif St == 3
    Pmat  = [[{'P'}  num2cell((1:St)) ;num2cell([(1:St)' results.Param.P])] repmat({''},4,5)];
end

% Print results.
Estimation_results = [Part_1; stats; repmat({'' '' ''},1,3); Pmat];
if St == 2
    Estimation_results = Estimation_results(:,1:6);
end
% Do table.
if print == 1
    fid = 1;
    fprintf(fid,'************************************************************************************\n');    
    display(Estimation_results);
    fprintf(fid,'************************************************************************************\n');
end

% Funtion print charts.
function charts = print_charts(results,dates,labels,St)
% Label for chart.
lab_dep_var = char(labels(1,:));
% Dates
if dates(end) == 1
    freq = 'm';
else
    freq = 'q';
end
[xTick,xTickLabel] = calendar(dates(1),dates(2),results.T,freq);

% Filter and Smoothed probabilities.
if St == 2
    figure(1)
    subplot(2,3,3)
    plot(results.Fil(:,1),'Color',[0 0 1]);
    hold on;
    plot(results.Fil(:,2),'Color',[1 0.1, 0]);
    ylim([0 1]);
    % xTickLabel Labels for chart.
    xlim([xTick(1) xTick(end)]);            
    set(gca,'xTick',xTick);
    set(gca,'xTickLabel', xTickLabel);
    title(('Filter Probability'),'FontWeight','bold','FontSize',11,'FontName','Arial');
    legend('Regime 1','Regime 2','Location','northwest');
    subplot(2,3,6)
    plot(results.Smt(:,1),'Color',[0 0 1]);
    hold on;
    plot(results.Smt(:,2),'Color',[1 0.1, 0]);
    ylim([0 1]);
    % xTickLabel Labels for chart.
    xlim([xTick(1) xTick(end)]);            
    set(gca,'xTick',xTick);
    set(gca,'xTickLabel', xTickLabel);
    title(('Smoothed Probability'),'FontWeight','bold','FontSize',11,'FontName','Arial');
    legend('Regime 1','Regime 2','Location','northwest');
    
elseif St == 3
    figure(1)
    subplot(2,3,3)
    plot(results.Fil(:,1),'Color',[0 0 1]);
    hold on;
    plot(results.Fil(:,2),'Color',[1 0.1, 0]);
    ylim([0 1]);
    plot(results.Fil(:,3),'Color',[1 1, 0.25]);
    ylim([0 1]);
    % xTickLabel Labels for chart.
    xlim([xTick(1) xTick(end)]);            
    set(gca,'xTick',xTick);
    set(gca,'xTickLabel', xTickLabel);
    title(('Filter Probability'),'FontWeight','bold','FontSize',11,'FontName','Arial');
    legend('Regime 1','Regime 2','Regime 3','Location','northwest');
    subplot(2,3,6)
    plot(results.Smt(:,1),'Color',[0 0 1]);
    hold on;
    plot(results.Smt(:,2),'Color',[1 0.1, 0]);
    hold on;
    plot(results.Smt(:,3),'Color',[1 1, 0.25]);
    ylim([0 1]);
    % xTickLabel Labels for chart.
    xlim([xTick(1) xTick(end)]);            
    set(gca,'xTick',xTick);
    set(gca,'xTickLabel', xTickLabel);  
    title(('Smoothed Probability'),'FontWeight','bold','FontSize',11,'FontName','Arial');    
    legend('Regime 1','Regime 2','Regime 3','Location','northwest');
end

% Do graph.
subplot(2,3,1)
plot(results.dta(:,1), '-k','LineWidth',0.75);
hold on
plot(results.yhat2(:,1),'Color',[0 0 1],'LineWidth',0.9);
hold on
plot(results.yhat2(:,2),'Color',[1 0.1, 0],'LineWidth',0.9);
if St == 3
    hold on
    plot(results.yhat2(:,3),'Color',[1 1, 0.25],'LineWidth',0.9);
end
% xTickLabel Labels for chart.
xlim([xTick(1) xTick(end)]);            
set(gca,'xTick',xTick);
set(gca,'xTickLabel', xTickLabel);
title('In sample fit of the model versus data','FontSize',11);
if St == 2
    legend(lab_dep_var,'Regime 1','Regime 2','Location','northwest');
elseif St == 3
    legend(lab_dep_var,'Regime 1','Regime 2','Regime 3','Location','northwest');
end
subplot(2,3,4)
plot(results.uhat, '-b');
hold on
plot(repmat(2*std(results.uhat),results.T,1), ':k');
hold on
plot(repmat(-2*std(results.uhat),results.T,1), ':k');
hold on
plot(zeros(results.T,1), '--k');
% xTickLabel Labels for chart.
xlim([xTick(1) xTick(end)]);            
set(gca,'xTick',xTick);
set(gca,'xTickLabel', xTickLabel);
title('Residuals +/- 2 S.D.','FontSize',11);

% Chart IRF.
subplot(2,3,2)

% Distribution of the residuals.
subplot(2,3,5)
h(1) = histogram(results.uhat2(:,1),31,'FaceColor',[0 0 1],'EdgeAlpha',1);
hold on
h(2) = histogram(results.uhat2(:,2),31,'FaceColor',[1 0 0],'EdgeAlpha',1);
legend('Regime 1','Regime 2','Location','northwest');
if St == 3
    hold on
    h(3) = histogram(results.uhat2(:,3),31,'FaceColor',[1 1, 0.25],'EdgeAlpha',1);
    legend('Regime 1','Regime 2','Regime 3','Location','northwest');
end
title('Histogram for residuals','FontSize',11);
charts = [];
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%