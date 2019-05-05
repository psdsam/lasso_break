function [KK, Lambda, history, b, mse] = lassocovtest(X,Y,ss0,alpha, flag, numb, trim, X_pre, Y_pre)
%This function is to implement the Least Angle Regression (LARS) proposed
%by Bradley Efron, Trevor Hastie, Iain Johnstone and Robert Tibshirani in
%"Least Angle Regression". Annals of Statistics. 32 (2): pp. 407–499 and
%Covariance Test by Richard Lockhart, Jonathan Taylor, Ryan J. Tibshirani,
%Robert Tibshirani in "A significance test for the lasso", Annals of
%Statistics. 42(2): pp. 413–468.

%The input (X, Y) are observations.
%ss0 is the inpute for MSE. If MSE is not known just leave it empty.
%alpha is the controled FDR rate
%flag choose type of nonparametric fit for x.
% 0 for no transform
% 1 for linear
% 2 for qudratic
% 3 3rd order poly
% 4 4th order poly
%numb is the maxium number of discountities allowed. Specify a large number
%to save the iteration in LARs
%trim a paramter to trim the support

% output:
%KK: the location of breaks;
%Lambda: the laombda that resulting the selected active set;\
%history: the LAR history

if trim==0
    dummy_idx = sort(X);
else
    dummy_idx = sort(X(X>quantile(X, trim) & X< quantile(X, 1-trim)));
end
p         = length(dummy_idx)-2;
n         = length(X);
X_hat   = zeros(n,p);
for i = 1:p
    X_hat(:,i) = (X>dummy_idx(i+1));
end


if flag==0
    Z_hat = ones(size(X));
    P_Z = eye(n,n);
    if ~isempty(X_pre)
        Z_pre = ones(size(X_pre));
    end
elseif flag==1
    Z_hat = [ones(size(X)), X];
    P_Z = eye(n,n) - Z_hat*((Z_hat'*Z_hat)\Z_hat');
    if ~isempty(X_pre)
        Z_pre = [ones(size(X_pre)), X_pre];
    end
elseif flag==2
    Z_hat = [ones(size(X)), X, X.^2];
    P_Z = eye(n,n) - Z_hat*((Z_hat'*Z_hat)\Z_hat');
    if ~isempty(X_pre)
        Z_pre = [ones(size(X_pre)), X_pre, X_pre.^2];
    end
elseif flag==3
    Z_hat = [ones(size(X)), X, X.^2, X.^3];
    P_Z = eye(n,n) - Z_hat*((Z_hat'*Z_hat)\Z_hat');
    if ~isempty(X_pre)
        Z_pre = [ones(size(X_pre)), X_pre, X_pre.^2, X_pre.^3];
    end
elseif flag==4
    Z_hat = [ones(size(X)), X, X.^2, X.^3, X.^4];
    P_Z = eye(n,n) - Z_hat*((Z_hat'*Z_hat)\Z_hat');
    if ~isempty(X_pre)
        Z_pre = [ones(size(X_pre)), X_pre, X_pre.^2, X_pre.^3, X_pre.^4];
    end
end

YY = P_Z*Y;
XX = P_Z*X_hat;

%for variance
if ~isempty(ss0)
    MSE_hat = ss0;
else
    [~, fitinfo] = lasso(XX,YY);
    MSE_hat = fitinfo.MSE(end);
end

[history, breaks_hat, Lambda,~] = mylarsp(YY, XX, MSE_hat, alpha, numb);

if ~isempty(X_pre)
    n_pre     = length(X_pre);
    X_hat_pre   = zeros(n_pre,p);
    for i = 1:p
        X_hat_pre(:,i) = (X_pre>dummy_idx(i+1));
    end
    if ~isempty(breaks_hat)
        KK = dummy_idx(breaks_hat+1);
        [b,~] = lasso(XX, YY,'Lambda',Lambda/sqrt(n), 'RelTol',1e-8);
        
        beta = ((Z_hat'*Z_hat)\Z_hat')*(Y-X_hat*b);
        mse = (Y_pre-Z_pre*beta-X_hat_pre*b)'*(Y_pre-Z_pre*beta-X_hat_pre*b);
    else
        KK = [];
        b = 0;
        beta = ((Z_hat'*Z_hat)\Z_hat')*(Y);
        mse = (Y_pre-Z_pre*beta)'*(Y_pre-Z_pre*beta);
    end
else
    if ~isempty(breaks_hat)
        KK = dummy_idx(breaks_hat+1);
        [b,~] = lasso(XX, YY,'Lambda',Lambda/sqrt(n), 'RelTol',1e-8);
        mse = 0;
    else
        KK = [];
        b = 0;
        mse = 0;
    end
end










