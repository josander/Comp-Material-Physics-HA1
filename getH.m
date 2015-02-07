function h = getH(alpha)

% Function that generates the matrix h(p,q) = <chi(p)|-nabla^2/2 - 2/r|chi(q)> 

h = zeros(4, 4);

    for p = 1:4
        for q = 1:4
            h(p, q) = alpha(q)*sqrt(pi)/(4*(alpha(p)+alpha(q))^(3/2))...
                    -alpha(q)*3*sqrt(pi)/(4*(alpha(p)+alpha(q))^(5/2))...
                    -1/(alpha(p)+alpha(q));
        end
    end

end