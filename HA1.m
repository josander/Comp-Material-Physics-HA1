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
C = [1, 2, 2, 1]';
C = C/norm(C);

% Construct the matrices h, S and C
h = getH(alpha);
S = getS(alpha);
Q = getQ(alpha);

% Construct the matrix F
F = getF(h, C, Q);

% Solve the generalised eigenvalue problem
Eigen = (F*C)\(S*C); % [4 x 4]x[4 x 1]\[4 x 4]x[4 x 1] = [4 x 1]\[4 x 1]



%% Task 2
clc
clear all


% Plot the Hartree potential
r = linspace(0,100,10000);
V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
plot(r, V(r));
xlabel('Radial distance r');
ylabel('The Hartree potential V');
