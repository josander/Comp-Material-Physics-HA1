function [ VSH ] = solveVSH( r, psi_r )
%SOLVEVHS Summary of this function goes here
%   r : radial grid with N ponits from [0,rMax]
%   psi_r : radial wavefunction over r, normalized
%           according to 4*pi*int(r^2psi_r^2 dr) = 1

[M N] = size(r);
A = zeros(N,N);
rMax = (r(end));
h = r(3) - r(2);

% Construct a, b and c
a = -2/(h^2);
b = 1/(h^2);
c = 1/(h^2);

% Construct the A solution
for i = 1:N-1
       A(i,i) = a;
       A(i,i+1) = b;
       A(i+1,i) = c;
end

% Implement the boudary conditions
A(1,1) = 1;
A(1,2) = 0;
A(end,end) = 1;
A(end, end -1) = 0;

% Get B matrix
B = -4*pi*r.*psi_r.^2;

% Solve system Au = B
u = A\B';

% Translate into the single Hartree potential
VSH = u'./r + 1/rMax;

plot(VSH,'r');
hold off


end

