function S = getS(alpha)

% Function that generates the matrix S(p,q) = <chi(p)|chi(q)> 

S = zeros(4, 4);

    for p = 1:4
        for q = 1:4
            S(p, q) = (pi/(alpha(p)+alpha(q)))^(3/2);
        end
    end

end