%% Task 4: Get ground state energy 
% ----- Find rMax convergence -----
clear all
clc

rMaxInit = 2;
rMaxFinal = 15;
dr = 0.5;

h = 0.005;

% FIND rMax-CONVERGENCE 
for rMax = rMaxInit:dr:rMaxFinal
    
    % Number of grid spaces
    N = 1 + rMax/h; 

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
    u = u/sqrt(trapz(x,u.^2));

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
        u = u/sqrt(trapz(x,u.^2));

        % Get the minimal eigenvalue in Hartree energy
        minEig = e(index);

        % Get energy in eV
        E = 27.211396132*minEig;

        % Calculate the new energy difference
        energyDiff = abs(Eold - E);

        % Save the solution
        Eold = E;

    end
    
    % Get ground state energy in Hartree
    Energy0 = 2*minEig - 2*trapz(x,u.^2.*Vsh/2);

    % Energy in eV
    EnergyEV = Energy0*27.211396132;
    
    % Save energy and rMax vectors to plot
    Energy((rMax-rMaxInit)/dr+1) = EnergyEV;
    RMax((rMax-rMaxInit)/dr+1) = rMax;
    
    rMax
    
end

save Task4rMax2.mat

%% ------ Find grid point convergence ------

clc
clear all

nPointsInit = 101;
nPointsFinal = 4501;
dn = 50;


% FIND GRIDPOINT-CONVERGENCE 
for N = nPointsInit:dn:nPointsFinal
    
    % Cutoff radius
    rMax = 15;

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
    u = u/sqrt(trapz(x,u.^2));

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
        u = u/sqrt(trapz(x,u.^2));

        % Get the minimal eigenvalue in Hartree energy
        minEig = e(index);

        % Get energy in eV
        E = 27.211396132*minEig;

        % Calculate the new energy difference
        energyDiff = abs(Eold - E);

        % Save the solution
        Eold = E;

    end

    % Get ground state energy in Hartree
    Energy0 = 2*minEig - 2*trapz(x,u.^2.*Vsh/2);

    % Energy in eV
    EnergyEV = Energy0*27.211396132;

    % Save energy and grid size vectors to plot them
    Energy((N-nPointsInit)/dn+1) = EnergyEV;
    gridSize((N-nPointsInit)/dn+1) = N;
    
    N
    
end

save Task4gridPoints.mat

%% ----- Get ground state energy Helium -----

clc
clear all

% Cutoff radius
rMax = 15;

% Number of points
N = 4001; 

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
u = u/sqrt(trapz(x,u.^2));

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
    u = u/sqrt(trapz(x,u.^2));

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
minEig

% Get ground state energy in Hartree
Energy0 = 2*minEig - 2*trapz(x,u.^2.*Vsh/2)

% Energy in eV
EnergyEV = Energy0*27.211396132

% Get wave function
waveFuncTask4 = A(:,index)'./x;

% Save all variables
save Task4.mat
