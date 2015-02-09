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
% solution is smaller or equal to 10^-5
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

%% Plot task 1
x = linspace(-5,5,1000);
phi = @(r) exp(-alpha(1)*r.^2).*C(1) + exp(-alpha(2)*r.^2).*C(2) + ...
    exp(-alpha(3)*r.^2).*C(3)+ exp(-alpha(4)*r.^2).*C(4);

plot(x,abs(phi(x)))

xlabel('Radial distance r');
ylabel('The wave function');

%% Task 2
clc
clf
clear all

% Cutoff radius
rMax = 50000;

% Number of points in the grid
N = 1001; 


% Radial, discetizised points 
x = linspace(0,rMax, N);

% Initialise two arrays with zeros
Y = zeros(N,1);
Ynew = zeros(N,1);

% Maximal difference in the solution, initially put to 1
maxDiff = 1;

% The step length between two points
h = rMax/(N-1);

% Single orbital density for the hydrogen atom
a0 = 1; % Bohr radius
Psi = @(r) 2*exp(-r/a0)/a0^(3/2); % Enligt Thijssen eq (3.23)
densConst = a0^(-3)/(pi);
densFunc = @(r) densConst*exp(-2*r/a0);


eDens = densFunc(x(2:end));
m = 0;

% Iterate until the convergence condition; the maximal difference in the
% solution is smaller or equal to 10^-3

%while maxDiff > 10^-10
while m < 5*10^5
    
    % Loop through the coordinates and calculate new solution
    for i = 2:N-1
        Y(i) = 2*pi*eDens(i)*x(i)*h^2 + 0.5*Y(i+1) + 0.5*Y(i-1);
    end

    % Maximal change in the solution compared to the last iteration
    %maxDiff = max(abs(Ynew - Y));

    % Save new solution
    %Y = Ynew;
    
    m= m+1;

end


%% Plot the Hartree potential

V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
Vsh = Y(2:end)'./x(2:end) + 1/rMax;
plot(x(2:end), V(x(2:end)), x(2:end), Vsh);
xlabel('Radial distance r');
ylabel('The Hartree potential V');

%% Task 3

clc
clf
clear all

% Cutoff radius
rMax = 50;

% Number of points
N = 101; 

% Radial, discetizised points 
x = linspace(10^(-9),rMax, N);

% Length between two points in the grid
h = rMax/N;

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

% The last element in the matrix
Y(N,N) = a(N);

% Solve the eigenvalue problem
[A B] = eig(Y);

% Get the eigenvalues
diag(B)

%% Task 4

