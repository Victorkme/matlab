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
f3 = figure;
autocorr(logretamzn.^2);
f4 = figure;
autocorr(abs(logretamzn));