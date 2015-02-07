function F = getF(Q)

F = zeros(4, 4);
sum = 0;

C = [1, 2, 2, 1];
h = zeros(4, 4);

for p = 1:4
    for q = 1:4
        
        for r = 1:4
            for s = 1:4
                sum = Q(p, r, q, s) * C(r) * C(s);
            end
        end
        
        F(p, q) = h(p, q) + sum;
        sum = 0;
        
    end % End q-loop
end % End p-loop

end