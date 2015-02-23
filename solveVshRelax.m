function [Vsh] = solveVshRelax(r, psi_r)

[M N] = size(r);
A = zeros(N,N);
rMax = (r(end));
h = r(3) - r(2);

nIterations = 100000;

Vsh = zeros(N,1);

for m = 1:nIterations

    % Loop through the coordinates and calculate new solution
    % Y(0) = Y(N) = 0
    for i = 2:N-1
        Vsh(i) = 2*pi*r(i)*psi_r(i)^2*h^2 + 0.5*Vsh(i+1) + 0.5*Vsh(i-1);
    end

end

Vsh = Vsh';
Vsh = Vsh(1:end)./r(1:end) + 1/rMax;

plot(Vsh,'b');
hold on
end
