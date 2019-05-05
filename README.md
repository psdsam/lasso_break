# lasso_break
Detect discontinuities using Lasso based covariance test as proposed in Backus and Peng (2019). 
If you use this package, please cite the article by Backus and Peng below:

[Matthew Backus and Sida Peng (2019). "On Testing Continuity and the Detection of Failures"](https://61cabf8f-a-ab15346e-s-sites.googlegroups.com/a/umich.edu/mrb/files/backuspeng_2018_wp_lassobreaks.pdf?attachauth=ANoY7cqUviWLqbH1vElYUmGy3YTBEGfmPUgZm1J-q5BJqH8BSWf-zCj6AOTkYcnt5--HepT-HrlNDmWJ5QHpjXzUbx2LQO3mu3Sy2PkDe5wYsSvE1tCIgZJoiS1UbNO50p6MvAEqjB2n-vQ1xfY5bo9ih3wCO1P_VmarpHCNmKGcb53DfhM56IP3YC4PPgb8S5QrXwHT2ms9l8E30EHIOC_4wYauuDbG-d8lDYYCUYcu-JuCgHk48Fo%3D&attredirects=0)

The functions in this package are:

detectdisc: This is the main function for the test using polynomial basis and 10-folded crossvalidation. The input is independent variable X, dependent variable Y and False Discovery Rate (FDR) control alpha. The out put is the detected location of discontinuities and plot.

lassocovtest: This is the main function for the test. Users can customize on the basis functions, cross-validation, maxmium number of discontinuities allowed and support trimming within this function.

mylarsp: This is the function implementing the algorithm proposed in the paper including a moddified combination of 
1) LARs algorithm by [Bradley Efron, Trevor Hastie, Iain Johnstone and Robert Tibshirani 
"Least Angle Regression"](http://statweb.stanford.edu/~tibs/ftp/lars.pdf) and code modified from [Sung Soo Kim's work](https://www.mathworks.com/matlabcentral/fileexchange/23186-lars-algorithm), 
2) Lasso covariance test by [Richard Lockhart, Jonathan Taylor, Ryan J. Tibshirani, 
%Robert Tibshirani "A significance test for the lasso"](http://statweb.stanford.edu/~tibs/ftp/covtest.pdf)
3) False Discovery Control by [Max Grazier G'Sell, Stefan Wager, Alexandra Chouldechova and Robert Tibshirani in "Sequential selection procedures and false discovery rate control"](http://www.stat.cmu.edu/~ryantibs/journalclub/gsell_2014.pdf). 
 
example: This scripts implment our method with one simulation example and one real data example. The plots for both examples are displayed below.

![alt text](https://github.com/psdsam/lasso_break/blob/master/plots/plot2.png)

