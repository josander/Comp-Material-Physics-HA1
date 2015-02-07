function h = getH(alpha)

% Function that generates the matrix h(p,q) = <chi(p)|-nabla^2/2 - 2/r|chi(q)> 

h = zeros(4, 4);

    for p = 1:4
        for q = 1:4
            h(p, q) = 1; % Whaaat?
        end
    end

end