function [ VEx ] = solveVEx(r, u  )
%SOLVEEX Summary of this function goes here
%   u : Normalised eigenfunction for heilum
%   VEx = ex +n *ex'
%   ex = - 3/4 * (3n/pi)^(1/3)

phi_r = u./(r.*sqrt(4*pi));

VEx = - (3*abs(phi_r).^2/pi).^(1/3);

end

