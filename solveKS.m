function [A B] = solveKS(pot, x)
% Function that solves the Khon-Sham equation. 'pot' is the potentials in the Hamiltonian except for -nabla/2. 

% Get number of gridpoints
[M N] = size(x);

% Initialise a matrix with zeros
Y = zeros(N,N);

% Get the length between two points
h = x(3) - x(2);

% Construct a, b and c
for i = 1:N
    a(i) = 1/h^2+pot(i);
end
b = - 1/(2*h^2);
c = - 1/(2*h^2);

% Construct the Y solution
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

end