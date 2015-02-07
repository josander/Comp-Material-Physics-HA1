function Q = getQ(alpha)

% The matrix of the kinetic and the Coulomb energy in a 
% Hartee-Fock approximation, with four bases. The indices p,
% r, q and s run from 1 to the number of basis functions.

Q = zeros(4, 4, 4, 4);

for p = 1:4
   
    for r = 1:4
       
        for q = 1:4
           
            for s = 1:4
                
                Q(p, r, q, s) = 2*pi^(5/2)/((alpha(p)+alpha(q))*...
                    (alpha(r)+alpha(s))*sqrt(alpha(p)+alpha(q)+...
                    alpha(r)+alpha(s)));
                
            end % End s-loop
            
        end % End q-loop
        
    end % End r-loop
    
end % End p-loop

end