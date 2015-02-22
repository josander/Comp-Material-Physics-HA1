function [ Ep ] = getEp( phi_r, c )
%GETEP Get epsilon
%   psi_r : Normalised wavefunction for heilum
%   c : option, 0 = just e_x, 1 = e_xc

[n m] = size(psi_r);
e_c = zeros(1,m); 

if c == 1
    A = -0.0311;
    B = -0.048;
    C = 0.002;
    D = -0.0116;

    gamma = -0.1423;
    beta1 = 1.0529;
    beta2 = 0.3334;
    r_s = (3./(4*pi*2*abs(psi_r).^2)).^(1/3);

    for i=1:m
       if r_s(i) < 1
           e_c(i) = A.*log(r_s(i))+B+C.*r_s(i).*log(r_s(i)) +D*r_s(i);
       else
           e_c(i) = gamma./(1 + beta1*sqrt(r_s(i)) + beta2*r_s(i)); 
       end
    end
end

Ep = -(3/4)*(3*2*abs(psi_r).^2/pi)^(1/3) + e_c;

end
