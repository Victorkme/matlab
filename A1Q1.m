% Assignment #1 - Question 1
% Author: Victor Me

% import data from CSV file
filename = 'AMZN.csv';
delimiterIn = ',';
headerlinesIn = 1;
amzn = importdata(filename,delimiterIn,headerlinesIn);

% Closing Prices
amznc = amzn.data(:,5);

% Calculate length of price vector
m = length(amznc);

% Plot asset prices
f1 = figure;
plot(t,amznc);
title('Amazon Historical Stock Prices (5-Yr)')
xlabel('Date (Days)')
ylabel('Stock Price ($)')

% Calculate log returns of asset
logamzn = log(amznc); %Log of prices
logretamzn = log(amznc(2:m)./amznc(1:m-1)); %log returns of prices
t1 = t(2:m); %loss of one observation

% Plot log returns
f2 = figure;
plot(t1,logretamzn);
title('Amazon Log-Returns (5-Yr)')
xlabel('Date (Days)')
ylabel('Log-Return')

% Calculate VaR at 0.95 and 0.99
V = 10^7; % Initial value of portfolio
loss = -V*logretamzn; % approximation based on taylor expansion

% Historic Approach
varhist = quantile(loss,[0.95 0.99]); %computes historic VaR

% Parametric Gaussian
mu = mean(logretamzn); %mean of log-returns
sigma = std(logretamzn); %std dev of log returns
varnorm95 = -V*mu+V*sigma*icdf('Normal',0.95,0,1); %at 0.95
varnorm99 = -V*mu+V*sigma*icdf('Normal',0.99,0,1); %at 0.99

% Parametric t-student (scaled and located)
nu = 4; %degrees of freedom
varstud95 = -V*mu+V*sigma*sqrt((nu-2)/nu)*icdf('T',0.95,nu); %at 0.95
varstud99 = -V*mu+V*sigma*sqrt((nu-2)/nu)*icdf('T',0.99,nu); %at 0.99

% Monte-Carlo VaR  - Gaussian
ns = 10^6; %number of simulations
mclossnorm = normrnd(-V*mu, V*sigma, 1, ns); %generating ns losses
varmcnorm = quantile(mclossnorm, [0.95 0.99]);

% Monte-Carlo VaR - t-student
mclossstud = -V*mu+V*sigma*sqrt((nu-2)/nu)*random('T',4,1,ns); %generating ns losses
varmcstud = quantile(mclossstud,[0.95 0.99]);

%Calculate ES at 0.95 and 0.99

%Historic Approach
condloss95 = (loss>varhist(1)); %losses greater than VaR at 95%
condloss99 = (loss>varhist(2)); %losses greater than VaR at 99%
eshist95 = sum(loss.*condloss95)/sum(condloss95); % ES at 95%
eshist99 = sum(loss.*condloss99)/sum(condloss99); % ES at 99%

%Parametric Gaussian
esnorm95 = -V*mu+V*sigma*pdf('Normal',icdf('Normal',0.95,0,1),0,1)/0.05; % ES at 95%
esnorm99 = -V*mu+V*sigma*pdf('Normal',icdf('Normal',0.99,0,1),0,1)/0.01; % ES at 99%

%Parametric t-student
esstud95 = -V*mu+(1/(1-0.95))*(nu+icdf('T',0.95,nu)^2)/(nu-1)*V*sigma*tpdf(icdf('T',0.95,nu),nu); % ES at 95%
esstud99 = -V*mu+(1/(1-0.99))*(nu+icdf('T',0.99,nu)^2)/(nu-1)*V*sigma*tpdf(icdf('T',0.99,nu),nu); % ES at 99%

%Monte-Carlo Gaussian
mccondloss95 = (mclossnorm>varmcnorm(1)); %losses greater than VaR 95%
mccondloss99 = (mclossnorm>varmcnorm(2)); %losses greater than VaR 99%
mceshist95 = sum(mclossnorm.*mccondloss95)/sum(mccondloss95); % ES at 95%
mceshist99 = sum(mclossnorm.*mccondloss99)/sum(mccondloss99); % ES at 95%

%Monte-Carlo t-student
mcstudloss95 = (mclossstud>varmcstud(1)); %losses greater than VaR 95%
mcstudloss99 = (mclossstud>varmcstud(2)); %losses greater than VaR 99%
mcstudeshist95 = sum(mclossstud.*mcstudloss95)/sum(mcstudloss95); % ES at 95%
mcstudeshist99 = sum(mclossstud.*mcstudloss99)/sum(mcstudloss99); % ES at 95%

%Backtesting Models

% Historic method
numhist95 = sum(condloss95);
numhist99 = sum(condloss99);

% Gaussian method
numnorm95 = sum((loss>varnorm95));
numnorm99 = sum((loss>varnorm99));

% t-student method
numstud95 = sum((loss>varstud95));
numstud99 = sum((loss>varstud99));

% Monte-Carlo Gaussian
nummcnorm95 = sum(loss>mceshist95);
nummcnorm99 = sum(loss>mceshist99);

% Monte-Carlo t-student
nummcstud95 = sum(loss>mcstudeshist95);
nummcstud99 = sum(loss>mcstudeshist99);

%Backtesting statistic

% Historic method
bthist95 = (numhist95-(m-1)*(1-0.95))/sqrt((m-1)*0.95*(1-0.95));
bthist99 = (numhist99-(m-1)*(1-0.99))/sqrt((m-1)*0.99*(1-0.99));

% Gaussian method
btnorm95 = (numnorm95-(m-1)*(1-0.95))/sqrt((m-1)*0.95*(1-0.95));
btnorm99 = (numnorm99-(m-1)*(1-0.99))/sqrt((m-1)*0.99*(1-0.99));

% t-student method
btstud95 = (numstud95-(m-1)*(1-0.95))/sqrt((m-1)*0.95*(1-0.95));
btstud99 = (numstud99-(m-1)*(1-0.99))/sqrt((m-1)*0.99*(1-0.99));

% Monte-Carlo Gaussian
btmcnorm95 = (nummcnorm95-(m-1)*(1-0.95))/sqrt((m-1)*0.95*(1-0.95));
btmcnorm99 = (nummcnorm99-(m-1)*(1-0.99))/sqrt((m-1)*0.99*(1-0.99));

% Monte-Carlo t-student
btmcstud95 = (nummcstud95-(m-1)*(1-0.95))/sqrt((m-1)*0.95*(1-0.95));
btmcstud99 = (nummcstud99-(m-1)*(1-0.99))/sqrt((m-1)*0.99*(1-0.99));


