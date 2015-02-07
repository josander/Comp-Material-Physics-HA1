function Q = getQelement(p, r, q, s)
% The matrix elements of the kinetic and the Coulomb energy in a 
% Hartee-Fock approximation, with four bases. The indices p,
% r, q and s has to be numbers in the range 1 to 4.


alpha = [0.297104, 1.236745, 5.749982, 38.216677];

Q = 2*pi^(5/2)/((alpha(p)+alpha(q))*(alpha(r)-alpha(s))*sqrt(alpha(p)+alpha(q)+alpha(r)+alpha(s)));

end