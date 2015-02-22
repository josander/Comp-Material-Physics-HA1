function [ VEx ] = solveVEx( psi_r )
%SOLVEEX Summary of this function goes here
%   psi_r : Normalised wavefunction for heilum
%   VEx = ex +n *ex'
%   ex = 3/4 * (3n/pi)^(1/3)

const = (-(3/pi)^(-2/3) - 3*(3/pi)^(1/3));

VEx = (2.*abs(psi_r).^2).^(1/3)/4;

end

