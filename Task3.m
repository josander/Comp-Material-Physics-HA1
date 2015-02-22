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

%% Plot the gound state wave function

clf
clc
set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);
plot(x, Psi(x)/Psi(x(1)), 'LineWidth', 1)
hold on
plot( x(2:10:end), waveFuncTask3(2:10:end)/waveFuncTask3(2),'--', 'color', 'red', 'LineWidth', 1)
axis([0 10 0 1]);

plotTickLatex2D

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave functions [-]','Interpreter','latex', 'fontsize', 12);    
title('Ground state wave function in hydrogen','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0]);

l = legend('Analytic wave function $\Psi(r)/\Psi(0)$','Wave function, Task 3 $\Psi_3(r)/\Psi_3(0)$');
set(l,'Interpreter','latex')

print(gcf,'-depsc2','task3.eps')

