%% Task 4: Get ground state energy 

clc

% FIND rMax-CONVERGENCE 
for r = 1:20
    
    % Number of points
    N = 1001; 

    % Radial, discetizised points 
    x = linspace(10^(-9),rMax, N);
    
    % Coefficients of wave function from task 1
    C = [-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];
    
    % Declaration of alpha
    alpha = [0.297104, 1.236745, 5.749982, 38.216677];

    % Initialise an array with zeros
    Psi = exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
        exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4);

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
    while energyDiff > 10^(-3) % [eV]

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
        energyDiff = abs(Eold - E);

        % Save the solution
        Eold = E;

    end

    Energy(r) = E
    
end

%%

clc
clear all

nPointsInit = 501;
nPointsFinal = 601;

% FIND GRIDPOINT-CONVERGENCE 
for N = nPointsInit:1:nPointsFinal
    
    % Cutoff radius
    rMax = 10;

    % Radial, discetizised points 
    x = linspace(10^(-9),rMax, N);

    % Coefficients of wave function from task 1
    C = [-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];
    
    % Declaration of alpha
    alpha = [0.297104, 1.236745, 5.749982, 38.216677];

    % Initialise an array with zeros
    Psi = exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
        exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4);

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
    while energyDiff > 10^(-3) % [eV]

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
        index = min(find(e == min(e)));
        
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

    Energy((N-nPointsInit)) = E
    
end

%% Plot the different energies with respect to the number of gridpoints

clf

plot(1:20,Energy)

%%

clc

% Cutoff radius
rMax = 10;

% Number of points
N = 1001; 

% Radial, discetizised points 
x = linspace(10^(-9),rMax, N);

% Coefficients of wave function from task 1
C = [-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];

% Declaration of alpha
alpha = [0.297104, 1.236745, 5.749982, 38.216677];

% Initialise wave function with solution from task 1
Psi = (exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
    exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4)).*x;

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


%%
clf
plot(A(:,index))