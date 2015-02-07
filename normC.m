function C = normC( C, S )
%NORMC Summary of this function goes here
%   Detailed explanation goes here

sum = 0;
for p = 1:4
    for q = 1:4
        
    sum = sum + C(p)*S(p,q)*C(q);
        
    end % End q-loop
end % End p-loop

C = C/sum;

end

