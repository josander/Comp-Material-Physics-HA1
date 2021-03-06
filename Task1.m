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
    [A, B]= eig(F, S);

    % Find the lowest real eigenvalue
    [x y] = find(B == min(diag(B)));

    % New C-array
    C = A(:,y); 

    % Normalize C via overlap maxtrix
    C = normC(C, S);

    % Get the energy of the state
    E  = getEG(h, C, Q);
    
    % Get energy in eV
    E = 27.211396132*E;

    % Calculate the new energy difference
    energyDiff = abs(Eold - E);
    
    % Save the solution
    Eold = E;

end


% Print the results
disp('The ground state energy in eV:');
E

% Print the results
disp('The ground state energy in Hartree:');
E = E/27.211396132

disp('Coefficients in wave func:')
C

save Task1.mat
