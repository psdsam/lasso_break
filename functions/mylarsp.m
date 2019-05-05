function [history, breaks_hat, Lambda, p] = mylarsp(YY, XX, MSE_hat, alpha, numb)
%This function combines the Least Angle Regression (LARS) proposed
%by Bradley Efron, Trevor Hastie, Iain Johnstone and Robert Tibshirani in 
%"Least Angle Regression". Annals of Statistics. 32 (2): pp. 407499;
%Matlab code is modified from Sung Soo Kim's work at 
%https://www.mathworks.com/matlabcentral/fileexchange/23186-lars-algorithm
%Covariance Test by Richard Lockhart, Jonathan Taylor, Ryan J. Tibshirani, 
%Robert Tibshirani in "A significance test for the lasso", Annals of
%Statistics. 42(2): pp. 413468;
%FDR control by Max Grazier G'Sell, Stefan Wager, Alexandra Chouldechova 
%and Robert Tibshirani in "Sequential %selection procedures and false 
%discovery rate control" Journal of the Royal Statistical Society. 78(2): 
%pp. %423-444.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalize
sx              = sum(XX.^2).^(1/2);
x               = XX ./ repmat(sx,size(XX,1),1);
y               = YY;
n               = size(x,1);        % # of samples
m               = size(x,2);        % # of predictors
all_candidate   = 1:m;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization

active      = [];               % active set
inactive    = all_candidate;    % inactive set
mu_a        = zeros(n,1);       % current estimate (eq. 2.8)
beta        = zeros(1,size(x,2));

history.active_set          = [];
history.add                 = [];
history.drop                = [];
history.MSE                 = sum(y.^2)/length(y);
history.lambda              = inf;
history.T                   = 0;
drop                        = [];               % used for 'lasso'
k                           = 1;                % iteration index
s                           = 0;
lambda                      = inf;              %tunning initilized at infty
TS                          = 0;                %test statistic
p                           = zeros(m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(numb)
    nn = inf;
else
    nn = numb;
end

while TS<alpha && k<=nn+2
    
    if k >= m-1
        break;
    end
    c                   = x'*(y-mu_a);                      % eq 2.8
    [C_max,C_max_ind]   = max(abs(c(inactive)));            % eq 2.9
    C_max_ind           = inactive(C_max_ind);
    
     %%lambda
    if ~isempty(active)
        Pa    = x(:, C_max_ind)'*(eye(n) -x(:,active)*((x(:,active)'*x(:,active))\x(:,active)'))*y;
        Pa_1  =  x(:, C_max_ind)'*x(:,active)*(((x(:,active)'*x(:,active)))\s);
        Pb    = ((x(:,active)'*x(:,active))\x(:,active)')*y;
        Pb_1  = (x(:,active)'*x(:,active))\s;
    else
        Pa    = x(:, C_max_ind)'*y;
        Pa_1  = 0;
        Pb    = 0;
        Pb_1  = 1;
    end

    lambda_2  = max(Pb./Pb_1.*(Pb./Pb_1<lambda-1e-8));
    lambda_1  = max(Pa/(1-Pa_1)*(Pa/(1-Pa_1)<lambda-1e-8), Pa/(-1-Pa_1)*(Pa/(-1-Pa_1)<lambda-1e-8));

    lambda    = max(lambda_1,lambda_2);
    T         = (lambda==lambda_1);

    active          = sort(union(active,C_max_ind));  
    inactive        = setdiff(all_candidate, active);   % eq 2.9
    
    if ~isempty(drop)
        C_max_ind= [];
    end
    active          = setdiff(active,drop);             % eq 3.6
    inactive        = sort(union(inactive,drop));       % eq 3.6
    s               = sign(c(active));                  % eq 2.10
    
    xa              = x(:,active).*repmat(s',n,1);  % eq 2.4
    ga              = xa'*xa;                       % eq 2.5
    
    invga           = ga\eye(size(ga,1));           % eq 2.5
    aa              = sum(sum(invga))^(-1/2);       % eq 2.5
    wa              = aa*sum(invga,2);              % eq 2.6
    ua              = xa*wa;                        % eq 2.6
    
    a               = x'*ua;                        % eq 2.11
    tmp_1           = (C_max - c(inactive))./(aa - a(inactive));
    tmp_2           = (C_max + c(inactive))./(aa + a(inactive));
    tmp_3           = [tmp_1, tmp_2];
    tmp             = tmp_3((tmp_3>0));
    gamma           = min(tmp);                     % eq 2.13
    if isempty(gamma) % if this is the last step (i.e. length(active)==maxKernels)
        gamma       = C_max/aa;                     % eq 2.19, eq 2.21 and 5 lines below eq 2.22
    end
    
    d               = zeros(1,m);
    d(active)       = s.*wa;
    
    tmp             = zeros(1,m);
    tmp(active)     = -1*beta(active)./d(active);   % eq 3.4
    tmp2            = tmp((tmp>0));
    
    %lasso drop
    drop            = [];
    gamma_tilde     = inf;                          % eq 3.5
    if ~isempty(tmp2) && gamma >= min(tmp2)
        gamma_tilde = min(tmp2);                    % eq 3.5
        drop        = find(tmp==gamma_tilde);       % eq 3.6
    end
    mu_a_plus       = mu_a + min(gamma, gamma_tilde)*ua;    % eq 3.6
    beta_new        = beta + min(gamma, gamma_tilde)*d;     % eq 3.3
    active          = setdiff(active,drop);                 % eq 3.6
    inactive        = setdiff(all_candidate,active);
    beta_new(drop)  = 0;
        
        
     
    mu_a_OLS        = mu_a + C_max/aa*ua;          % eq 2.19, 2.21
    MSE             = sum((y - mu_a_OLS).^2)/length(y);  
    s               = sign(c(active));  %update s after dropping

    
    % update and save
    mu_a = mu_a_plus;
    beta = beta_new;
    
    k = k+1;
    history(k).active_set   = active;
    history(k).drop         = drop;
    history(k).add          = C_max_ind;
    history(k).MSE          = MSE;
    history(k).lambda       = lambda;
    history(k).T            = T;
    
    
    %test
    if ~isempty(history(k-1).add)
        XOLS = XX(:,history(k-1).active_set');
        ss_test = sign((XOLS'*XOLS)\XOLS'*YY);
        s_test = ss_test((history(k-1).active_set'==history(k-1).add));
        
        
        if ~isempty(history(k-2).active_set)
            [B_2,~] = lasso(XX(:,[history(k-2).active_set, history(k-2).drop]), YY,'Lambda',history(k).lambda/sqrt(n), 'RelTol',1e-8);
            A_activeset = [history(k-2).active_set, history(k-2).drop];
            J_add = history(k-1).add;
            sA = sign(B_2);
            
            XX1 = x(:,[A_activeset,J_add]);
            XX2 = x(:,A_activeset);
            C1 = XX1*((XX1'*XX1)\[sA;s_test]);
            C2 = XX2*((XX2'*XX2)\sA);
            C_test = (C1-C2)'*(C1-C2);
            
        else
            C_test = 1;
        end
        %covariance test
        T = C_test*history(k-1).lambda*(history(k-1).lambda-history(k).lambda)/MSE_hat;
        p(k-2) = (1-fcdf(T,2,n-m));
        TS = -1/(k-2)*sum(log(1-p(1:k-2)));
    end    
end

khat = k-2;
if khat>0
    Lambda         =history(khat+1).lambda;
    breaks_hat     = history(khat).active_set; 
    p              = p(khat);
else
    Lambda         = [];
    breaks_hat     = [];
    p              = [];
end
    
    
    
    
    
    
    
    


