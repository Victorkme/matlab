% Assignment #2 - Question 1
% Author: Victor Me

% import data from CSV file
filename = 'AMZN.csv';
delimiterIn = ',';
headerlinesIn = 1;
amzn = importdata(filename,delimiterIn,headerlinesIn);

% Closing Prices
amznc = amzn.data(:,5);

% Calculate length of price vector & time
m = length(amznc);
t = datetime(amzn.textdata(2:m+1,1));

% Calculate log returns of asset
logamzn = log(amznc); %Log of prices
logretamzn = log(amznc(2:m)./amznc(1:m-1)); %log returns of prices
t1 = datetime(t(2:m)); %loss of one observation

% Plot log returns
f1 = figure;
plot(t1,logretamzn);
title('Amazon Log-Returns (5-Yr)')
xlabel('Date (Days)')
ylabel('Log-Return')

% Compute autocorrelation function
f2 = figure;
autocorr(logretamzn);
title('Amazon Autocorrelation Log-Returns')
f3 = figure;
autocorr(logretamzn.^2);
title('Amazon Autocorrelation Squared Log-Returns')
f4 = figure;
autocorr(abs(logretamzn));
title('Amazon Autocorrelation Absolute Value Log-Returns')

% Fit GARCH(1,1) with Gaussian noise

model = garch('Offset',NaN,'GARCHLags',1,'ARCHLags',1)
estMdl = estimate(model,logretamzn);
[Vn,En] = simulate(estMdl,m,'NumPaths',5);
f5 = figure
subplot(2,1,1)
plot(Vn)
xlim([0,100])
title('Conditional Variances (Gaussian)')
subplot(2,1,2)
plot(En)
xlim([0,100])
title('Innovations (whitenoise - Gaussian)')

% Fit GARCH(1,1) with t-student noise
tdist = struct('Name','t','DoF',NaN)
modelt = garch('Offset',NaN,'GARCHLags',1,'ARCHLags',1,'Distribution',tdist)
estMdlt = estimate(modelt,logretamzn);
[Vnt,Ent] = simulate(estMdlt,m,'NumPaths',5);
f6 = figure
subplot(2,1,1)
plot(Vnt)
xlim([0,100])
title('Conditional Variances (t-dist)')
subplot(2,1,2)
plot(Ent)
xlim([0,100])
title('Innovations (whitenoise - t-dist)')

% Calculate VaR & ES using GARCH(1,1) Gaussian white noise
V = 10^9;
sigmat = 0.0005; %given std dev of log returns
returnt = logretamzn(end); % last log return
sigmat1 = sqrt(estMdl.Constant + estMdl.ARCH{1,1}*returnt^2+estMdl.GARCH{1,1}*sigmat^2); % estimate sigma t+1
varnorm95 = V*sigmat1*icdf('Normal',0.95,0,1); % VaR at 0.95
esnorm95 = V*sigmat1*pdf('Normal',icdf('Normal',0.95,0,1),0,1)/0.05; % ES at 0.95

% Calculate VaR using GARCH(1,1) t-student white noise
nu = estMdlt.Distribution.DoF; %degrees of freedom
V = 10^9;
sigmatstud = 0.0005; %given std dev of log returns
returntstud = logretamzn(end); % last log return
sigmatstud1 = sqrt(estMdlt.Constant + estMdlt.ARCH{1,1}*returntstud^2+estMdlt.GARCH{1,1}*sigmatstud^2); % estimate sigma t+1
varstud95 = V*sigmatstud1*icdf('T',0.95,nu); % VaR at 0.95
esstud95 = (1/(1-0.95))*(nu+icdf('T',0.95,nu)^2)/(nu-1)*V*sigmatstud1*tpdf(icdf('T',0.95,nu),nu); % ES at 0.95