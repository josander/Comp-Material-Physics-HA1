%% Task 4: Get ground state energy 

clc

rMaxInit = 1;
rMaxFinal = 20;
dr = 1;

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
    U0 = exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
        exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4);

    % Length between two points
    h = rMax/(N-1);

    % Initialise a matrix with zeros
    Y = zeros(N,N);

    % Number of relaxations for the single Hartree potential
    nRelax = 5000;

    % Variable to keep track of the energy difference
    energyDiff = 1;
    Eold = 0;

    % Iterate until the convergence condition; the maximal difference in the
    % solution is smaller or equal to 10^-6
    while energyDiff > 10^(-3) % [eV]
        
        % Get the single Hartree potential
        V = solveVSH(x, U0);

        % Define the potential
        pot = -2./x+V;

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
plot(RMax,Energy);

%%

clc
clear all

nPointsInit = 501;
nPointsFinal = 3001;
dn = 100;

% FIND GRIDPOINT-CONVERGENCE 
for N = nPointsInit:dn:nPointsFinal
    
    % Cutoff radius
    rMax = 10;

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

    % Initialise a matrix with zeros
    Y = zeros(N,N);

    % Number of relaxations for the single Hartree potential
    nRelax = 5000;

    % Variable to keep track of the energy difference
    energyDiff = 1;
    Eold = 0;

    % Iterate until the convergence condition; the maximal difference in the
    % solution is smaller or equal to 10^-6
    while energyDiff > 10^(-3) % [eV]

        % Get the single Hartree potential
        V = solveVSH(x, U0);

        % Define the potential
        pot = -2./x+V;

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
plot(gridSize,Energy,'.')

%%

clc

% Cutoff radius
rMax = 30;

% Number of points
N = 2001; 

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

% Initialise a matrix with zeros
Y = zeros(N,N);

% Number of relaxations for the single Hartree potential
nRelax = 50000;

% Variables to keep track of the energy difference
energyDiff = 1;
Eold = 0;

% Iterate until the convergence condition; the maximal difference in the
% solution is smaller or equal to 10^-6
while energyDiff > 10^(-8) % [eV]

    % Get the single Hartree potential
    V = solveVSH(x, U0);
    
    % Define the potential
    pot = -2./x+V;
    
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


%%
clf
plot(A(:,index))

