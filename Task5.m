%% Task 5

clc
clear all

% Cutoff radius
rMax = 6;

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
    
    % Get the exchange potential
    Vx = solveVEx(x,u);
    
    % Define the potential
    pot = -2./x+2*Vsh+Vx;
    
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

%Print value of lowest eigenvalue
minEig = minEig

% Get eigenvalue of exchange function
epsilonX = getEp(x, u, 0);

% Get the single Hartree potential
Vsh = solveVSH(x, u);
    
% Get the exchange potential
Vx = solveVEx(x,u);

% Get ground state energy in Hartree
Energy0 = 2 * minEig - 2* trapz(x, u.^2.*(Vsh + Vx - epsilonX))

% Energy in eV
EnergyEV = Energy0*27.211396132

% Get wave function
waveFuncTask5 = A(:,index)'./x;

save Task5.mat

