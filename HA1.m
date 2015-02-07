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
    Eold = E;


    pause(1);
end


disp('The ground state energy:');
E
disp('C:')
C
disp('alpha:')
alpha


%% Task 2
clc
clear all

rMax = 50;
N = 1000 + 2; %Number of points
x = linspace(0,rMax, N - 1);
Y = zeros(N,1);
Ynew = zeros(N,1);
maxError = 1;
h = rMax/N;


eDens = @(r) r*(r - rMax);

while maxError > 10^-3
    
    Y(3) = Y(1);
    Y(N-2) = Y(N);
    for i = 2:N-1
        Ynew(i) = 2*pi*eDens(x(i))*x(i)*h^2 + 0.5*Y(i+1) + 0.5*Y(i-1);

    end
    maxError = max(abs(Ynew - Y));
        
    Y = Ynew;
        
    
end

% Plot the Hartree potential
r = linspace(0,rMax,10000);
V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
plot(r, V(r));
xlabel('Radial distance r');
ylabel('The Hartree potential V');

%%
