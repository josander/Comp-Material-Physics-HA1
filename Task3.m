%% Task 3: Solve Khon-Sham by eigenvalues

clc
clf
clear all

% Borh radius
a0 = 1;

% Cutoff radius
rMax = 50;

% Number of points
N = 2001; 

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

waveFuncTask3 = abs(A(:,index)/sqrt(h))./x;

save Task3.mat

%% Plot the gound state wave function

clf
set(gcf,'renderer','painters','PaperPosition',[0 0 6 3]);
plot(x, Psi(x), 'LineWidth', 2)
hold on
plot( x, waveFunkTask3,'.', 'color', 'red')
axis([0 10 0 0.8]);

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave functions [-]','Interpreter','latex', 'fontsize', 12);    
title('Ground state wave function in hydrogen','Interpreter','latex', 'fontsize', 14);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

l = legend('Analytic wave function $\Psi(r)/\Psi(0)$','Numerical wave function $\Psi(r)/\Psi(0)$');
set(l,'Interpreter','latex')

plotTickLatex2D

print(gcf,'-depsc2','task3.eps')

