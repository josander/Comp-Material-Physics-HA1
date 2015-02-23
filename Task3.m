%% Task 3: Solve Khon-Sham by eigenvalues

clc
clf
clear all

% Borh radius
a0 = 1;

% Cutoff radius
rMax = 15;

% Number of points
N = 3001; 

% Radial, discetizised points 
x = linspace(10^(-9),rMax, N);

%p length between two points
h = rMax/(N-1);

% Define the potential
pot = -1./x;

% Solve the Khon-Sham equation to get the eigenvalues and the eigenvectors 
[A B] = solveKS(pot, x);

% Get the eigenvalues
e = (diag(B));

% Find index of the minimal eigenvalue
index = find(e == min(e));

% Get the minimal eigenvalue in Hartree energy
minEig = e(index)

% Get energy in eV
Energy = 27.211396132*minEig

% Analytic wave function for ground state hydrogen
Psi = @(r) 2/(a0^(3/2))*exp(-r/a0);

waveFuncTask3 = abs(A(:,index)/sqrt(h))./x';

save Task3.mat

