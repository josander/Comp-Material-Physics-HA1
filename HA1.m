%% Home assignment 1.
% Task 1

clc
clear all

% Initialise new matrices with zeros
h = zeros(4, 4);
S = zeros(4, 4);
Q = zeros(4, 4, 4, 4);
F = zeros(4, 4);

% Declaration of alpha
alpha = [0.297104, 1.236745, 5.749982, 38.216677];

% Initialisation of C
C = [1, .5, .5, 1]';

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
% solution is smaller or equal to 10^-5
while energyDiff > 10^(-5) % [eV]

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
disp('C:')
C
disp('alpha:')
alpha


%% Task 2
clc
clear all

% Cutoff radius
rMax = 50;

% Number of points
N = 1001; 

% Radial, discetizised points 
x = linspace(0,rMax, N);

% Initialise two arrays with zeros
Y = zeros(N,1);
Ynew = zeros(N,1);

% Maximal difference in the solution, initially put to 1
maxDiff = 1;
% The step length between two points
h = rMax/N;

% Electron density for the hydrogen atom
a_0 = 1; % Bohr radius

eDens = @(r) 2*r.^2.*exp(-2*r./a_0)/(a_0^4);


% Iterate until the convergence condition; the maximal difference in the
% solution is smaller or equal to 10^-3
while maxDiff > 10^-3
    

    % Loop through the coordinates and calculate new solution
    for i = 2:N-1
        Ynew(i) = eDens(x(i))*x(i)*h^2 + 0.5*Y(i+1) + 0.5*Y(i-1);
    end
   
    % Maximal change in the solution compared to the last iteration
    maxDiff = max(abs(Ynew - Y));
    % Save new solution
    Y = Ynew;
    
end


% Plot the Hartree potential


V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
Y = Y %+ x'/rMax;
plot(x, V(x), x,Y);
xlabel('Radial distance r');
ylabel('The Hartree potential V');

%% Task 3
