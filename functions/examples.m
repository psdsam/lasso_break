%this script implements two examples. The first example is a real data
%example from Lee (2008) U.S. house election data. The second example is a
%simulation example with exp continous function + two jump discontinuities.

clear all
clc

alpha = 0.05;

%%example1
fid      = fopen('lee_2008_data.txt');
formatSpec = '%f%f';
paramIds = textscan(fid,formatSpec,'HeaderLines',1,'Delimiter',',','EmptyValue',-Inf);
fclose(fid);

difdemshare = paramIds{1};
demsharenext = paramIds{2};

meshgrid = linspace(-0.25,0.25,101);

voteshare = zeros(length(meshgrid)-1,1);
for i = 1:length(meshgrid)-1
    voteshare(i) = mean(demsharenext(difdemshare>meshgrid(i) & difdemshare<meshgrid(i+1)));
end
X = difdemshare(difdemshare<0.25 & difdemshare>-0.25);
Y = demsharenext(difdemshare<0.25 & difdemshare>-0.25);
KK = detectdisc(X, Y, alpha);

%%example2
n=2000;
X = randn(n,1);
breaks_J = [-0.5; 1];
Y = -exp(X)+randn(n,1);
for i = 1:length(breaks_J)
    Y = Y+(X>breaks_J(i)).*5;
end
KK = detectdisc(X, Y, alpha)  ;  
    
    
    
    
    
    
    
    
    
    
    
    
    
