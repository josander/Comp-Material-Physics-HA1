%% Home assignment 1.
% Task 1

clc
clear all
format long

% Initialise new matrices with zeros
h = zeros(4, 4);
S = zeros(4, 4);
Q = zeros(4, 4, 4, 4);
F = zeros(4, 4);

% Declaration of alpha
alpha = [0.297104, 1.236745, 5.749982, 38.216677];

% Initialisation of C
C = [1, 1, 1, 1]';

% Construct the matrices h, S and C
h = getH(alpha);
S = getS(alpha);
Q = getQ(alpha);

% Variable to keep track of the energy difference
energyDiff = 1;
Eold = 0;

% Normalize C via overlap maxtrix
C = normC(C, S);

% Iterate until the convergence condition; the maximal difference in the
% solution is smaller or equal to 10^-6
while energyDiff > 10^(-6) % [eV]

    % Construct the matrix F
    F = getF(h, C, Q);

    % Solve the generalised eigenvalue problem
    [A B]= eig(F, S);

    % Find the lowest real eigenvalue
    [x y] = find(B == min(diag(B)));

    % New C-array
    C = A(:,y); 

    % Normalize C via overlap maxtrix
    C = normC(C, S);

    % Get the energy of the state
    E  = getEG(h, C, Q);

    % Calculate the new energy difference
    energyDiff = abs(Eold - E);
    
    % Save the solution
    Eold = E;

end

% Print the results
disp('The ground state energy:');
E
disp('Coefficients in wave func:')
C

%% Plot the found wave function
x = linspace(-5,5,1000);
phi = @(r) exp(-alpha(1)*r.^2).*C(1) + exp(-alpha(2)*r.^2).*C(2) + ...
    exp(-alpha(3)*r.^2).*C(3)+ exp(-alpha(4)*r.^2).*C(4);

plot(x,abs(phi(x)))

xlabel('Radial distance r');
ylabel('The wave function');


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
Ynew = zeros(N,1);

% The step length between two points
h = rMax/(N-1);

% Single orbital density for the hydrogen atom
a0 = 1; % Bohr radius
Psi = @(r) 2*exp(-r/a0)/a0^(3/2); % Enligt Thijssen eq (3.23)
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


%% Plot the Hartree potentials
clf
clc

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
Vsh = Y(2:end)'./x(2:end) + 1/rMax;
plot(x(2:end), V(x(2:end)),'.');
hold on
plot( x(2:end), Vsh,'--', 'MarkerSize', 12, 'Color', 'red')


set(gcf,'renderer','painters','PaperPosition',[0 0 6 3]);
X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
%y = ylabel('PDF [1/$a_0$]','Interpreter','latex', 'fontsize', 12);    

title('Electron potential in hydrogen','Interpreter','latex', 'fontsize', 14);

%set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

l = legend('Analytic Hartree potential $V_H$','Numerical single Hartree potential $V_{sH}$');
plotTickLatex2D
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
set(l,'Interpreter','latex')
plotTickLatex2D
print(gcf,'-depsc2','task2.eps')

%% Task 3
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


%% Task 4
clc

% Cutoff radius
rMax = 10;

% Number of points
N = 1001; 

% Radial, discetizised points 
x = linspace(10^(-9),rMax, N);

% Initialise an array with zeros
Psi = zeros(1,N);

% Length between two points
h = rMax/(N-1);

% Initialise a matrix with zeros
Y = zeros(N,N);

% Number of relaxations for the single Hartree potential
nRelax = 50000;

% Variable to keep track of the energy difference
energyDiff = 1;
Eold = 0;

% Iterate until the convergence condition; the maximal difference in the
% solution is smaller or equal to 10^-6
while energyDiff > 10^(-5) % [eV]

    % Get the single Hartree potential
    V = getVSH(N, rMax, nRelax, Psi);

    % Construct a, b and c
    for i = 1:N
        a(i) = 1/h^2-2/x(i)+V(i);
    end
    b = - 1/(2*h^2);
    c = - 1/(2*h^2);

    % Construct the solution
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
    
    % 
    Psi = A(:,index)';

    % Get the minimal eigenvalue in Hartree energy
    minEig = e(index);

    % Get energy in eV
    E = 27.211396132*minEig;

    % Calculate the new energy difference
    energyDiff = abs(Eold - E)
    
    % Save the solution
    Eold = E;

end

Energy = E
