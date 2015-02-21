%% Task 4: Get ground state energy 

clc

rMaxInit = 5;
rMaxFinal = 50;
dr = 2;


% FIND rMax-CONVERGENCE 
for rMax = rMaxInit:dr:rMaxFinal
    
    % Number of grid spaces
    N = 2001; 

    % Radial, discetizised points 
    x = linspace(10^(-9),rMax, N);
    
    % Coefficients of wave function from task 1
    C = [-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];
    
    % Declaration of alpha
    alpha = [0.297104, 1.236745, 5.749982, 38.216677];

    % Initialise an array with zeros
    U0 = exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
        exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4);

    % Length between two points
    h = rMax/(N-1);

    % Variable to keep track of the energy difference
    energyDiff = 1;
    Eold = 0;

    % Iterate until the convergence condition; the maximal denergy difference
    % 10^-5
    while energyDiff > 10^(-5) % [eV]
        
        % Get the single Hartree potential
        Vsh = solveVSH(x, U0);

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
        U0 = A(:,index)';

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
    Energy((rMax-rMaxInit)/dr+1) = E;
    RMax((rMax-rMaxInit)/dr+1) = rMax; 
    
    rMax
    
end


%% Plot the energies with respect to rMax
clf
clc

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

plot(RMax,Energy,'-', 'MarkerSize', 12, 'Color', 'red');

X = xlabel('Cut off radius $r_max$ [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Ground state energy [eV]','Interpreter','latex', 'fontsize', 12);    
title('Ground state energy for Helium','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

plotTickLatex2D
print(gcf,'-depsc2','convRMax.eps')

%%

clc
clear all

nPointsInit = 101;
nPointsFinal = 4001;
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

    % Initialise an array with zeros
    U0 = exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
        exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4);

    % Length between two points
    h = rMax/(N-1);

    % Variable to keep track of the energy difference
    energyDiff = 1;
    Eold = 0;

    % Iterate until the convergence condition; the maximal denergy difference
    % 10^-5
    while energyDiff > 10^(-5) % [eV]

        % Get the single Hartree potential
        Vsh = solveVSH(x, U0);
set(gcf,'renderer','painters','PaperPosition',[0 0 6 3]);
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
        U0 = A(:,index)';

        % Get the minimal eigenvalue in Hartree energy
        minEig = e(index);

        % Get energy in eV
        E = 27.211396132*minEig;

        % Calculate the new energy difference
        energyDiff = abs(Eold - E);

        % Save the solution
        Eold = E;

    end

    Energy((N-nPointsInit)/dn+1) = E;
    gridSize((N-nPointsInit)/dn+1) = N;
    
    N
    
end

%% Plot the different energies with respect to the number of gridpoints

clf
clc

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

plot(gridSize,Energy,'-', 'MarkerSize', 12, 'Color', 'red');

X = xlabel('Grid size  N [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Energy [eV]','Interpreter','latex', 'fontsize', 12);    
title('Ground state energy for Helium','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

plotTickLatex2D
print(gcf,'-depsc2','convGrid.eps')

%%

clc

% Cutoff radius
rMax = 15;

% Number of points
N = 3001; 

% Radial, discetizised points 
x = linspace(10^(-9),rMax, N);

% Coefficients of wave function from task 1
C = -1*[-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];

% Declaration of alpha
alpha = [0.297104, 1.236745, 5.749982, 38.216677];

% Initialise the wave function with solution from task 1
U0 = (exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
    exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4)).*x;

% Normalise U0 4pi int(r^2U0^2) = 1  
U0 = U0/sqrt(trapz(4*pi.*x.^2.*U0.^2));

% Length between two points
h = rMax/(N-1);

% Variables to keep track of the energy difference
energyDiff = 1;
Eold = 0;

% Iterate until the convergence condition; the maximal denergy difference
% 10^-5
while energyDiff > 10^(-5) % [eV]

    % Get the single Hartree potential
    Vsh = solveVSH(x, U0);
    
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
    U0 = A(:,index)';

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
waveFuncTask4 = A(:,index)'./x;

%% Plot the wave functions

clf
clc

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

plot(x(2:end),psi_r(2:end)./(x.*psi_r(2)))
hold on
plot(x(2:end), waveFuncTask4(2:end)./waveFuncTask4(2),'--', 'MarkerSize', 12, 'Color', 'red');

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave function [-]','Interpreter','latex', 'fontsize', 12);    

title('Electron potential in hydrogen','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

l = legend('Normalised wave function from Task 1 $\Psi_1(r)/\Psi_1(0)$','Normalised wave function from Task 4 $\Psi_4(r)/\Psi_4(0)$');
set(l,'Interpreter','latex')
plotTickLatex2D

print(gcf,'-depsc2','task4.eps')
