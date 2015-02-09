function F = getF(h, C, Q)

F = zeros(4, 4);
sum = 0;

for p = 1:4
    for q = 1:4
        
        for r = 1:4
            for s = 1:4
                sum = sum + Q(p, r, q, s) * C(r) * C(s);
            end
        end
        
        F(p, q) = h(p, q) + sum;
        sum = 0;
        
    end % End q-loop
end % End p-loop

end