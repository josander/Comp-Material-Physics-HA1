%% Home assignment 1
% Task 1

clc
clear all

% Declaration of alpha
alpha = {0.297104, 1.236745, 5.749982, 38.216677};

r = linspace(1,100);

chi = @(r) exp(-alpha(1).*r.^2);


%% Task 2
clc
clear all


% Plot the Hartree potential
r = linspace(0,100,10000);
V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
plot(r, V(r));
xlabel('Radial distance r');
ylabel('The Hartree potential V');

%%
tau = 10^(-12);
o = 10^9;

ratio = (1+tau^2*o^2)- tau*o^2*(1+tau^2*o^2)/(1-tau*o^2*(1+tau^2*o^2))^(1/4)