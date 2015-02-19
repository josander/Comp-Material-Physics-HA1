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

% Initialise a matrix with zeros
Y = zeros(N,N);

% Construct a, b and c
for i = 1:N
    a(i) = 1/h^2-1/x(i);
end
b = - 1/(2*h^2);
c = - 1/(2*h^2);

% Construct the Y solution
for i = 1:N-1
       Y(i,i) = a(i);
       Y(i,i+1) = b;
       Y(i+1,i) = c;
end
% Implement the boundary conditions
Y(1,1) = 1;
Y(1,2) = 0;
Y(end,end-1) = 0;
Y(end,end) = 1;

% Solve the eigenvalue problem
[A B] = eig(Y);

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


%% Plot the gound state wave function

clf
set(gcf,'renderer','painters','PaperPosition',[0 0 6 3]);
plot(x, Psi(x).*x, 'LineWidth', 2)
hold on
plot( x, abs(A(:,index)/sqrt(h)),'.', 'color', 'red')
axis([0 10 0 0.8]);

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
%y = ylabel('PDF [1/$a_0$]','Interpreter','latex', 'fontsize', 12);    
title('Ground state radial wave function in hydrogen','Interpreter','latex', 'fontsize', 14);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
%set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
l = legend('Analytic wave function','Numerical wave function');
set(l,'Interpreter','latex')

plotTickLatex2D
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

print(gcf,'-depsc2','task3.eps')

