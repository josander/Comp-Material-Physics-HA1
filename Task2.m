%% Task 2: Get single Hartree potential

clc
clf
clear all

% Cutoff radius
rMax = 50;

% Number of points in the grid
N = 1001; 

% Radial, discetizised points 
x = linspace(0,rMax, N);

% Initialise an array with zeros
Y = zeros(N,1);

% The step length between two points
h = rMax/(N-1);

% Single orbital density for the hydrogen atom
a0 = 1; % Bohr radius
densConst = a0^(-3)/(pi);
densFunc = @(r) densConst*exp(-2*r/a0);

diffConst = 2*pi*h^2;
eDens = densFunc(x(2:end));
eDens = eDens.*x(2:end)*diffConst;

nIterations = 5*10^5;

for m = 1:nIterations

    % Loop through the coordinates and calculate new solution
    % Y(0) = Y(N) = 0
    for i = 2:N-1
        Y(i) = eDens(i-1) + 0.5*Y(i+1) + 0.5*Y(i-1);
    end

end

save Task2.mat

