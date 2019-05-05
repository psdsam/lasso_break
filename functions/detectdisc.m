function KK = detectdisc(X, Y, alpha)

%implement crossvalidation for continous function g
n = length(Y);
indics = crossvalind('Kfold', n, 10);
mse = zeros(1,5);
for kk = 1:10
    test = indics ==kk; train = ~test;
    X_pre = X(test);
    Y_pre = Y(test);
    X_train = X(train);
    Y_train = Y(train);
    for jj = 0:4
        [~, ~, ~, ~, mmse] = lassocovtest(X_train,Y_train, [], alpha, jj, 5, 0.05, X_pre, Y_pre);
        mse(jj+1) = mse(jj+1)+mmse;
    end
end

[~, idx] = min(mse);


[KK, ~, ~, ~, ~] = lassocovtest(X,Y, [], alpha, idx-1, 5, 0.05, [], []);

if ~isempty(KK)
    KK1 = [min(X); sort(KK); max(X)];
    for ii=1:(length(KK1)-1)
        X1 = X(X>=KK1(ii) & X<KK1(ii+1));
        Y1 = Y(X>=KK1(ii) & X<KK1(ii+1));
        fit1.(sprintf('struct%d', ii)) = polyfit(X1, Y1, idx-1);
    end
else
    fit = polyfit(m, DW, idx-1);
end

%plot
grid = linspace(min(X), max(X), 101);
M = zeros(1,100);
for i  = 1:100
    M(i) = mean(Y((X>=grid(i)&X<grid(i+1))));
end

figure
plot(grid(3:end), M(2:end), 'o')
hold on
if ~isempty(KK)
    for ii=1:(length(KK1)-1)
        X1 = X(X>=KK1(ii) & X<KK1(ii+1));
        p1 = plot(sort(X1), polyval(fit1.(sprintf('struct%d', ii)),sort(X1)), 'Color',[0.85,0.325,0.098]);
        %plot(sort(X2), polyval(fit2,sort(X2)), 'Color',[0.85,0.325,0.098]);
    end
else
    p1 = plot(sort(X), polyval(fit,sort(X)), 'Color',[0.85,0.325,0.098]);
end

yyl = get(gca, 'ylim');
for i = 1:length(KK)
    p2 = line([KK(i) KK(i)], yyl,'Color','k','LineStyle',':','Marker','x');
end

title('Detection of Jump Discontinuity Plot');

xlim([min(X),max(X)])
yyaxis right
[f, xi] = ksdensity(X, grid);
p3 = plot(xi, f,'--', 'Color',[0.466,0.674,0.188]);
if ~isempty(KK)
    legend([p1,p3,p2],{'Polynomial fit','Density plot','Discont detect'},'Location','northwest')
else
    legend([p1,p3],{'Polynomial fit','Density plot'})
end

set(findall(gca,'Type','Line'), 'LineWidth', 2)