function h = getH(alpha)

% Function that generates the matrix h(p,q) = <chi(p)|-nabla^2/2 - 2/r|chi(q)> 

h = zeros(4, 4);

    for p = 1:4
        for q = 1:4
            h(p, q) = 3*alpha(q)*pi^(3/2)/((alpha(p)+alpha(q))^(3/2))...
                    -alpha(q)^2*3*pi^(3/2)/((alpha(p)+alpha(q))^(5/2))...
                    -4*pi/(alpha(p)+alpha(q));
        end
    end

end