function Cnorm = normTest( C, S )

%NORMC Normalizes the C-vector 

Cnorm = C/sqrt(C'*S*C);

end
