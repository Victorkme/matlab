% Question1
% graph the 2d gumbel copula with theta1=1 theta2=2
clear all
close all
u = linspace(0,1,50);
[u1,u2]=meshgrid(u,u);
theta1 = 1;
theta2 = 3;
y1 = copulacdf('Gumbel',[u1(:),u2(:)],theta1);
y2 = copulacdf('Gumbel',[u1(:),u2(:)],theta2);

figure(1)
subplot(2,1,1)
surf(u1,u2,reshape(y1,50,50))
title('Gumbel, \theta=1')
xlabel('u1')
ylabel('u2')
subplot(2,1,2)
surf(u1,u2,reshape(y2,50,50))
title('Gumbel, \theta=2')
xlabel('u1')
ylabel('u2')

% Question 1 b)
% Graph the density of a Gumbel copula with theta 1,2

y3 = copulapdf('Gumbel',[u1(:),u2(:)],theta1);
y4 = copulapdf('Gumbel',[u1(:),u2(:)],theta2);
figure(2)
subplot(2,1,1)
surf(u1,u2,reshape(y3,50,50))
title('Gumbel, \theta=1')
xlabel('u1')
ylabel('u2')
subplot(2,1,2)
surf(u1,u2,reshape(y4,50,50))
title('Gumbel, \theta=2')
xlabel('u1')
ylabel('u2')

% Question 1 c)
% Generate 500 pairs of numbers from the Gumbel Copula

n = 500;
U3 = copularnd('Gumbel',theta1,n)
U4 = copularnd('Gumbel',theta2,n)
figure(3)
subplot(2,1,1)
plot(U3(:,1),U3(:,2),'.')
title('Gumbel, \theta=1')
xlabel('u1')
ylabel('u2')
subplot(2,1,2)
plot(U4(:,1),U4(:,2),'.')
title('Gumbel, \theta=2')
xlabel('u1')
ylabel('u2')

%Q2a
%Generate 500 points from a Gaussian copula
%with marginal uniform cdf.
%rho = 0.7
rho = 0.7;
rng('default')
U1 = copularnd('Gaussian',[1,rho;rho,1],n);
figure(4)
subplot(2,1,1)
plot(U1(:,1),U1(:,2),'.')
rho = -0.7;
rng('default')
U2 = copularnd('Gaussian',[1,rho;rho,1],n);
figure(4)
subplot(2,1,2)
plot(U2(:,1),U2(:,2),'.')

%Q2b
%generate 500 points from gaussian copula
%with exponential marginals

rho1=0.7;
Z = mvnrnd([0,0],[1,rho1;rho1,1],n);
U = normcdf(Z);
X = [expinv(U(:,1),1),expinv(U(:,2),1)];
figure(5)
scatterhist(X(:,1),X(:,2),'Direction','Out')

%Q3
%load stockreturns.mat
%a) plot scatterplot
load stockreturns
x = stocks(:,1);
y = stocks(:,2);
figure(6)
scatterhist(x,y)
%b) rescale the data to have uniformly distributed marginals from [0,1]
u = ksdensity(x,x,'function','cdf');
v = ksdensity(y,y,'function','cdf');
figure(7)
scatterhist(u,v)

%c) fit a t-student copula to the data
%maximum likelihood method
[rho,nu]=copulafit('t',[u,v],'Method','ApproximateML')
