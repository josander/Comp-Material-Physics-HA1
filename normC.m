function Cnorm = normC( C, S )

%NORMC Normalizes the C-vector 

Cnorm = C./sqrt(C'*S*C);

end

