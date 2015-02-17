function [ VSH ] = getVSH( N, rMax, nIter, psi)
%GETVSH Summary of this function goes here
%   


x = linspace(rMax/N,rMax, N);

Y = zeros(1,N+2);
u_sq = zeros(1,N); 

h = rMax/(N);

u_sq = 4*pi.*x.*psi.^2 *2*pi*h^2; % u^2/r * 2 *pi *h^2 

for m = 1:nIter

    % Loop through the coordinates and calculate new solution
    % Y(0) = Y(N) = 0
    for i = 2:N+1
        Y(i) = u_sq(i-1) + 0.5*Y(i+1) + 0.5*Y(i-1);
    end

end
VSH = Y(2:end-1)./x + 1/rMax;


end
