function  E  = getEG(h, C,Q )
%GETEG Summary of this function goes here
%   Detailed explanation goes here

E = 0;

for p = 1:4
    for q = 1:4
        
        E = E + 2*C(p)*C(q)*h(p,q);
        
        for r = 1:4
            for s = 1:4
                sum = Q(p, r, q, s) * C(p) * C(q)* C(r) * C(s);
            end
        end
        
        E = E + sum;
        sum = 0;
        
    end % End q-loop
end % End p-loop



end

