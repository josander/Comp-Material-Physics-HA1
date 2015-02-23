%% Task 4: Get ground state energy 
clear all
clc

rMaxInit = 1;
rMaxFinal = 20;
dr = 0.5;

h = 0.005;


% FIND rMax-CONVERGENCE 
for rMax = rMaxInit:dr:rMaxFinal
    
    % Number of grid spaces
    N = 1 + rMax/h; 

    % Radial, discetizised points 
    x = linspace(10^(-9),rMax, N);
    
    % Coefficients of wave function from task 1
    C = [-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];
    
    % Declaration of alpha
    alpha = [0.297104, 1.236745, 5.749982, 38.216677];

    % Initialise an array with zeros
    psi_r = (exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
        exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4)).*x;

    % Normalise psi_r : 4pi int(r^2U0^2) = 1  
    psi_r = psi_r/sqrt(trapz(4*pi.*x.^2.*psi_r.^2));

    % Length between two points
    h = rMax/(N-1);

    % Variable to keep track of the energy difference
    energyDiff = 1;
    Eold = 0;

    % Iterate until the convergence condition; the maximal denergy difference
    % 10^-5
    while energyDiff > 10^(-5) % [eV]
        
        % Get the single Hartree potential
        Vsh = solveVSH(x, psi_r);

        % Define the potential
        pot = -2./x+Vsh;

        % Solve the Khon-Sham equation and get the eigenvalues and the
        % eigenvectors
        [A B] = solveKS(pot, x);

        % Get the eigenvalues
        e = (diag(B));

        % Find index of the minimal eigenvalue
        index = min(find(e == min(e)));

        % The new radial wave function
        psi_r = A(:,index)';
        
        % Normalise psi_r : 4pi int(r^2U0^2) = 1  
        psi_r = psi_r/sqrt(trapz(4*pi.*x.^2.*psi_r.^2));

        % Get the minimal eigenvalue in Hartree energy
        minEig = e(index);

        % Get energy in eV
        E = 27.211396132*minEig;

        % Calculate the new energy difference
        energyDiff = abs(Eold - E);

        % Save the solution
        Eold = E;
        
    end

    % Save energy and rMax
    Eigen((rMax-rMaxInit)/dr+1) = minEig;
    RMax((rMax-rMaxInit)/dr+1) = rMax; 
    
    rMax
    
end

save Task4rMax2.mat

%% Plot the energies with respect to rMax
clf
clc



set(gcf,'renderer','painters','PaperUnits','centimeters','PaperPosition',[0 0 12 8]);

plot(RMax,Energy/27.211396132,'r-', 'LineWidth', 1);
axis([0 20 -2.3 .5])
X = xlabel('Cutoff radius $r_{max}$ [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Eigenvalue [$E_h$]','Interpreter','latex', 'fontsize', 12);    
title('Convergence with respect to cutoff radius','Interpreter','latex', 'fontsize', 14);
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);


print(gcf,'-depsc2','convRMax.eps') 

%%

clc
clear all

nPointsInit = 101;
nPointsFinal = 4501;
dn = 50;


% FIND GRIDPOINT-CONVERGENCE 
for N = nPointsInit:dn:nPointsFinal
    
    % Cutoff radius
    rMax = 15;
    
    % Grid density
    h = rMax/(N-1);

    % Radial, discetizised points 
    x = linspace(10^(-9),rMax, N);
    
    % Coefficients of wave function from task 1
    C = -1*[-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];

    % Declaration of alpha
    alpha = [0.297104, 1.236745, 5.749982, 38.216677];

    % Initialise psi_r
    psi_r =( exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
        exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4)).*x;
    

    % Normalise psi_r : 4pi int(r^2U0^2) = 1  
    psi_r = psi_r/sqrt(trapz(4*pi.*x.^2.*psi_r.^2));

    % Length between two points
    h = rMax/(N-1);

    % Variable to keep track of the energy difference
    energyDiff = 1;
    Eold = 0;

    % Iterate until the convergence condition; the maximal denergy difference
    % 10^-5
    while energyDiff > 10^(-5) % [eV]

        % Get the single Hartree potential
        Vsh = solveVSH(x, psi_r);

        % Define the potential
        pot = -2./x+Vsh;

        % Solve the Khon-Sham equation and get the eigenvalues and the
        % eigenvectors
        [A B] = solveKS(pot, x);

        % Get the eigenvalues
        e = (diag(B));

        % Find index of the minimal eigenvalue
        index = min(find(e == min(e)));

        % The new radial wave function
        psi_r = A(:,index)';

        % Normalise psi_r : 4pi int(r^2U0^2) = 1  
        psi_r = psi_r/sqrt(trapz(4*pi.*x.^2.*psi_r.^2));

        % Get the minimal eigenvalue in Hartree energy
        minEig = e(index);

        % Get energy in eV
        E = 27.211396132*minEig;

        % Calculate the new energy difference
        energyDiff = abs(Eold - E);

        % Save the solution
        Eold = E;

    end
   

    Eigen((N-nPointsInit)/dn+1) = minEig;
    gridSize((N-nPointsInit)/dn+1) = N;
    
    N
    
end

save Task4gridPoints.mat

%% Plot the different energies with respect to the number of gridpoints

clf
clc

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

plot(gridSize,Eigen,'-', 'Color', 'red', 'LineWidth', 1);
%axis([0 4000 -55 -15]);
plotTickLatex2D

X = xlabel('Grid points  N [-]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Eigenvalue [$E_h$]','Interpreter','latex', 'fontsize', 12);    

title('Convergence with respect to number of grid points','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0]);


print(gcf,'-depsc2','convGrid.eps')

%%

clc
clear all

% Cutoff radius
rMax = 15;

% Number of points
N = 1001; 

% Radial, discetizised points 
x = linspace(10^(-9),rMax, N);

% Coefficients of wave function from task 1
C = -1*[-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];

% Declaration of alpha
alpha = [0.297104, 1.236745, 5.749982, 38.216677];

% Initialise the wave function with solution from task 1
psi_r = (exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
    exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4)).*x;

% Get u from psi_r
u = sqrt(4*pi)*x.*psi_r;

% Normalise u
u = u/sqrt(trapz(u.^2));

% Length between two points
h = rMax/(N-1);

% Variables to keep track of the energy difference
energyDiff = 1;
Eold = 0;

% Iterate until the convergence condition; the maximal energy difference
% 10^-5
while energyDiff > 10^(-5) % [eV]

    % Get the single Hartree potential
    Vsh = solveVSH(x, u);
    
    % Define the potential
    pot = -2./x+Vsh;
    
    % Solve the Khon-Sham equation and get the eigenvalues and the
    % eigenvectors
    [A B] = solveKS(pot, x);

    % Get the eigenvalues
    e = (diag(B));

    % Find index of the minimal eigenvalue
    index = min(find(e == min(e)));
    
    % The new radial wave function
    u = A(:,index)';
    
    % Normalise
    u = u/sqrt(trapz(u.^2));

    % Get the minimal eigenvalue in Hartree energy
    minEig = e(index);

    % Get energy in eV
    E = 27.211396132*minEig;

    % Calculate the new energy difference
    energyDiff = abs(Eold - E)
    
    % Save the solution
    Eold = E;

end

% Print value of lowest eigenvalue
minEig = minEig

% Get ground state energy in Hartree
Energy0 = 2*minEig - 2*trapz(u.^2.*Vsh/2)

% Energy in eV
EnergyEV = E

% Get wave function
waveFuncTask4 = A(:,index)'./x;

% Save all variables
save Task4.mat

%% Plot the wave functions

clf
clc

psi = psi_r./x;

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

plot(x(2:end),psi(2:end)/psi(2), 'LineWidth', 1)
hold on
plot(x(2:end), waveFuncTask4(2:end)./waveFuncTask4(2),'--', 'LineWidth', 1, 'Color', 'red');
axis([0 5 0 1]);

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave function [-]','Interpreter','latex', 'fontsize', 12);    

plotTickLatex2D
title('Wave function for Helium','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0]);

l = legend('Wave function, Task 1 $\Psi_1(r)/\Psi_1(0)$','Wave function, Task 4 $\Psi_4(r)/\Psi_4(0)$');
set(l,'Interpreter','latex')


print(gcf,'-depsc2','task4.eps')
