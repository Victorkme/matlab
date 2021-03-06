% Assignment #1 - Question 2 & 3
% Author: Victor Me

% import data from CSV file
filename1 = 'AMZN.csv';
filename2 = 'TSLA.csv';
filename3 = 'RY.csv';
delimiterIn = ',';
headerlinesIn = 1;
amzn = importdata(filename1,delimiterIn,headerlinesIn);
tsla = importdata(filename2,delimiterIn,headerlinesIn);
ry = importdata(filename3,delimiterIn,headerlinesIn);

% Closing prices
amznc = amzn.data(:,5);
tslac = tsla.data(:,5);
ryc = ry.data(:,5);

% Calculate length of price vectors
m1 = length(amznc);
m2 = length(tslac);
m3 = length(ryc);

% Calculate log returns of the three assets
logamzn = log(amznc); %Log of prices
logretamzn = log(amznc(2:m1)./amznc(1:m1-1)); %log returns of prices
logtsla = log(tslac); %Log of prices
logrettsla = log(tslac(2:m2)./tslac(1:m2-1)); %log returns of prices
logry = log(ryc); %Log of prices
logretry = log(ryc(2:m3)./ryc(1:m3-1)); %log returns of prices

% Weights of assets in portfolio
w1 = [1/3,1/3,1/3];
w2 = [1/2,1/4,1/4];
V = 10^7;

% expected log returns, volatility, covariance, correlation matrix
mu = [mean(logretamzn),mean(logrettsla),mean(logretry)];
sigmax = [std(logretamzn),std(logrettsla),std(logretry)];
covar = cov([logretamzn, logrettsla, logretry]);
corr = corrcoef([logretamzn, logrettsla, logretry]);
nu = 4;

% Estimated log-returns of portfolio
mup = w1*mu';
mup2 = w2*mu';
sigmap = sqrt(w1*covar*w1');
sigmap2 = sqrt(w2*covar*w2');

% Calculate VaR of portfolio
varnorm95 = -V*mup+V*sigmap*icdf('Normal',0.95,0,1); %at 0.95
varnorm99 = -V*mup+V*sigmap*icdf('Normal',0.99,0,1); %at 0.99
varw2norm95 = -V*mup2+V*sigmap2*icdf('Normal',0.95,0,1); %at 0.95
varw2norm99 = -V*mup2+V*sigmap2*icdf('Normal',0.99,0,1); %at 0.99

% Calculate ES of portfolio
esnorm95 = -V*mup+V*sigmap*pdf('Normal',icdf('Normal',0.95,0,1),0,1)/0.05; % ES at 95%
esnorm99 = -V*mup+V*sigmap*pdf('Normal',icdf('Normal',0.99,0,1),0,1)/0.01; % ES at 99%
esw2norm95 = -V*mup2+V*sigmap2*pdf('Normal',icdf('Normal',0.95,0,1),0,1)/0.05; % ES at 95%
esw2norm99 = -V*mup2+V*sigmap2*pdf('Normal',icdf('Normal',0.99,0,1),0,1)/0.01; % ES at 99%

% Generate the first 4 moments of log-returns
muamzn = mean(logretamzn);
varamzn = var(logretamzn);
skewamzn = skewness(logretamzn);
kurtamzn = kurtosis(logretamzn);

% Estimate Empirical pdf and cdf
[p,x1] = hist(logretamzn,50);
[f,x] = ecdf(logretamzn);

% Normal pdf and cdf with data mu and sigma
ncdf = normcdf(x,muamzn,std(logretamzn));
npdf = normpdf(x1,muamzn,std(logretamzn));

% pdf graph
f3 = figure;
plot(x,f)
hold on
plot(x,ncdf,'r')
legend('ECDF','Normal CDF')
title('Empirical CDF vs Normal CDF')
xlabel('Log-Returns')
ylabel('Cumulative Density Function')
hold off

% cdf graph
f4 = figure;
hold on
plot(x1,p/sum(p));
plot(x1,npdf/sum(npdf),'r');
legend('EPDF','Normal PDF')
title('Empirical PDF vs Normal PDF')
xlabel('Log-Returns')
ylabel('Probability Density Function')
hold off

%kstest statistics at 5%
test_ncdf = makedist('Normal','mu',muamzn,'sigma',std(logretamzn));
test_tcdf = makedist('tlocationscale','mu',muamzn,'sigma',(sigmax(1)*sqrt((nu-2)/nu)),'nu',nu);
tpdf = pdf('tlocationscale',x1,muamzn,std(logretamzn),4);
[pdfhist1,x2] = hist(logretamzn,50);
[h,p,kstat,cv] = kstest(logretamzn,'CDF',test_ncdf);
[h2,p2,kstat2,cv2] = kstest(logretamzn,'CDF',test_tcdf);

%Graph of distributions
f5 = figure;
hold on
plot(x2,pdfhist1/sum(pdfhist1),'b');
plot(x2,npdf/sum(npdf),'r');
plot(x2,tpdf/sum(tpdf),'g');
legend('EPDF','Normal PDF','t-student pdf')
title('Empirical PDF vs Normal PDF vs t-student')
xlabel('Log-Returns')
ylabel('Probability Density Function')