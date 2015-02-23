function [ VEx ] = solveVEx( u )
%SOLVEEX Summary of this function goes here
%   psi_r : Normalised wavefunction for heilum
%   VEx = ex +n *ex'
%   ex = - 3/4 * (3n/pi)^(1/3)

VEx = - (6*abs(u).^2/pi).^(1/3);

end

