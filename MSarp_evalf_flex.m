function [logf,e0,sm,param,yhat,yhat2,uhat] = MSarp_evalf_flex(theta,y,x,x1,x2,x3,Pr)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 29/Nov/2018
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Evaluation of the log-likelihood of an MRAR(p) model.
%
% Imputs: 
%   theta   : Value for the vector of parameter.
%   y       : Endogenous variable.
%   x       : Independet variables.
%   x1      : Independet variables for regime 1.
%   x2      : Independet variables for regime 2.
%   x3      : Independet variables for regime 3.
%   Pr      : Matrix with unit entries fro restriction P matrix.
%
% Outputs:
%   logf : Log-like function.
%   e0   : Filter probabilities.
%   sm   : Smoothed probabilities.
%   param: parameters for each regime
%   yhat : fit of the model.
%   yhat2: fit from each regime.
%   uhat : in-sample residuals.
% 
% Index:
% 1. Initial Setup
% 2. Evaluation.
% 3. Filter and Smoothed Probabilities.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 1. Initial Setup
St = (1-isempty(x1))+(1-isempty(x2))+(1-isempty(x3));
T  = size(y,1);

% Computing Betas and Sg2 for each version of the model.
nx = size(x,2); nx1 = size(x1,2); nx2 = size(x2,2); nx3 = size(x3,2); 
Bx = theta(1:nx);
B1 = theta(nx+1:nx+nx1);
B2 = theta(nx+nx1+1:nx+nx1+nx2);
B3 = [];
% Getting B3 and theta2.
if St == 3
    B3 = theta(nx+nx1+nx2+1:nx+nx1+nx2+nx3);
end

% Transition Matrix definition.
p = zeros(St,St);
% No restrictions on P matrix
if Pr == ones(St);
    if St == 2
        p(1,1) = theta(end-1).^(2) / (1 + theta(end-1).^(2));
        p(2,1) = 1 - p(1,1);
        p(1,2) = theta(end).^(2) / (1 + theta(end).^(2));
        p(2,2) = 1 - p(1,2);
        sg2    = (theta(nx+nx1+nx2+1:end-2).^(2))';
    elseif St == 3
        p(1,1) = theta(end-5).^(2) / (1 + theta(end-5).^(2) + theta(end-2).^(2));
        p(2,1) = theta(end-2).^(2) / (1 + theta(end-5).^(2) + theta(end-2).^(2));
        p(3,1) = 1 - (p(1,1) + p(2,1));
        p(1,2) = theta(end-4).^(2) / (1 + theta(end-4).^(2) + theta(end-1).^(2));
        p(2,2) = theta(end-1).^(2) / (1 + theta(end-4).^(2) + theta(end-1).^(2));
        p(3,2) = 1 - (p(1,2) + p(2,2));
        p(1,3) = theta(end-3).^(2) / (1 + theta(end-3).^(2) + theta(end).^(2));
        p(2,3) = theta(end).^(2) / (1 + theta(end-3).^(2) + theta(end).^(2));
        p(3,3) = 1 - (p(1,3) + p(2,3));
        sg2    = (theta(nx+nx1+nx2+nx3+1:end-6).^(2))';
    end

% Restrictions on P matrix
else        
    temp  = sum(Pr);    
    if St == 2
        temp0 = 2;
        P_aux = theta(end-1:end);
    elseif St == 3
        temp0 = sum(sum(Pr))-St;
        P_aux = theta(end-temp0+1:end);
     end
    % Building theta.
    theta_m = [];
    j = 1;
    for i0 = 1:St
        if temp(i0) == 1
            P_aux2(1) = 1;
        elseif temp(i0) == 2
            P_aux2(1) = P_aux(j).^(2) / (1 + P_aux(j).^(2));
            P_aux2(2) = 1 - P_aux2(1);
            j = j + 1;
        elseif temp(i0) == 3
            P_aux2(1) = P_aux(j).^(2)   / (1 + P_aux(j).^(2) + P_aux(j+1).^(2));
            P_aux2(2) = P_aux(j+1).^(2) / (1 + P_aux(j).^(2) + P_aux(j+1).^(2));
            P_aux2(3) = 1 - (P_aux2(1) + P_aux2(2));
            j = j + 2;
        end
        theta_m = [theta_m; P_aux2'];
        clear P_aux2;
    end
    % Getting P matrix.
    j = 1;
    for i0 = 1:St;
        for i1 = 1:St
            if Pr(i1,i0) == 1
                p(i1,i0) = theta_m(j);
                j = j + 1;
            end
        end
    end
    clear j i0 i1;
    
    % Getting sg2
    sg2 = (theta(nx+nx1+nx2+nx3+1:end-temp0).^(2))';
end
    
% Getting Sg2
if size(sg2,2) == 1
    sg2 = repmat(sg2,1,St);
end

% Ergodic probabilities, initial values for the filter prob.
A = [eye(size(p,1)) - p; ones(1,size(p,1))] ;
e = (A'*A)\(A')*eye(St+1);
e = e(:,end);
clear A;

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 2. Evaluation.
% Restricting the log-like.
if max(abs(p(:))) > 1 || min(p(:)) < 0 || round(sum(sum(p))) ~= St
    logf = inf;

else    
    % Log-Likelihood evaluation.
    part1  = repmat((2*pi*sg2).^(-1/2),size(y,1),1);
    part21 = exp( - (y-[x x1]*[Bx; B1]).^2 .* repmat(0.5./sg2(1),size(y,1),1));
    part22 = exp( - (y-[x x2]*[Bx; B2]).^2 .* repmat(0.5./sg2(2),size(y,1),1));
    if St == 3
        part23 = exp( - (y-[x x3]*[Bx; B3]).^2 .* repmat(0.5./sg2(3),size(y,1),1));
    else
        part23 = [];
    end
    part2  = [part21 part22 part23]; 
    ll    = part1 .* part2;
    clear part1 part21 part22 part23 part2;
    % Conditional density and uptating of e(j,t).
    e0 = zeros(T+1,St);
    e0(1,:) = e;
    logf = 0;
    for i = 1:T
        % Updating (forecasting) the state (t+1|t).
        if i == 1
            e1 = e0(i,:);
            fi = e1*ll(i,:)';
        else
            e1 = e0(i,:)*p';
            fi = e1*ll(i,:)';
        end
        % Conditional density of yt
        logf = logf + log(fi);
        % Filter Probability (t|t)
        e0(i+1,:) = (e1.*ll(i,:))/fi;
    end
    % log-likelihood.
    logf = -1*logf;
end;
clear i e1 fi;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 3. Filter and Smoothed Probabilities.
if nargout > 1
    % Filter probability
    e0 = e0(2:end,:);
    % Forecast of the state (t+1|t).
    e1 = zeros(T,St);
    for i = 1:T
        e1(i,:) = (p*e0(i,:)')';
    end;
    % Smoothed probabilities.    
    sm(T,:)=e0(T,:);
    for i = T-1:-1:1
        sm(i,:)= e0(i,:)'.*((p'*(sm(i+1,:)./ e1(i,:))'));
    end; 
    % Parameters.
    param.Bx  = Bx;
    param.B1  = B1;
    param.B2  = B2;
    param.B3  = B3;
    param.Sg2 = sg2;
    param.P   = p;
    % Computing residuals.
    if St == 2
        yhat = sum(([([x x1]*[Bx; B1]) ([x x2]*[Bx; B2])]).*e1,2);
        yhat2= ([([x x1]*[Bx; B1]) ([x x2]*[Bx; B2])]).*e1;
        uhat = sum(([(y-[x x1]*[Bx; B1]) (y-[x x2]*[Bx; B2])]).*e1,2);
    elseif St == 3
        yhat = sum(([([x x1]*[Bx; B1]) ([x x2]*[Bx; B2]) ([x x3]*[Bx; B3])]).*e1,2);
        yhat2= ([([x x1]*[Bx; B1]) ([x x2]*[Bx; B2]) ([x x3]*[Bx; B3])]).*e1;
        uhat = sum(([(y-[x x1]*[Bx; B1]) (y-[x x2]*[Bx; B2]) (y-[x x3]*[Bx; B3])]).*e1,2);
    end
end;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%