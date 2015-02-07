function Cnorm = normC( C, S )
%NORMC Normalizes the C-vector 
%   Detailed explanation goes here

Cnorm = C/sqrt(C'*S*C);

end

