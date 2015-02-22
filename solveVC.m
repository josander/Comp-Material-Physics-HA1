function [ VC ] = solveVC( psi_r )
%SOLVEVC Summary of this function goes here
%   Detailed explanation goes here

A = -0.0311;
B = -0.048;
C = 0.002;
D = -0.0116;

gamma = -0.1423;
beta1 = 1.0529;
beta2 = 0.3334;


[n m] = size(psi_r);

r_s = (3./(4*pi*2*abs(psi_r).^2)).^(1/3);
e_c = zeros(1,m); 
de_c = zeros(1,m);

for i=1:m
   if r_s(i) < 1
       e_c(i) = A.*log(r_s(i))+B+C.*r_s(i).*log(r_s(i)) +D*r_s(i);
       de_c(i) = (A + C*r_s(i) + D*r_s(i))/r_s(i) + C*log(r_s(i));
   else
       e_c(i) = gamma./(1 + beta1*sqrt(r_s(i)) + beta2*r_s(i));
       de_c(i) = -gamma.*(beta1.*(1./2*sqrt(r_s(i))) + beta2)/(1 + beta1*sqrt(r_s(i)) + beta2*r_s(i)).^2;
   end

       
       
    
end
e_c(1:index-1) = A.*log(r_s(1:index-1))+B+C.*r_s(1:index-1).*log(r_s(1:index-1)) +D*r_s(1:index-1);
e_c(index:end) = gamma./(1 + beta1*sqrt(r_s(index:end)) + beta2*r_s(index:end));

de_c(1:index-1) = (A + C*r_s(1:index-1) + D*r_s(1:index-1))/r_s(1:index-1) + C*log(r_s(1:index-1));
de_c(index:end) = -gamma.*(beta1.*(1./2*sqrt(r_s(index:end))) + beta2)/(1 + beta1*sqrt(r_s(index:end)) + beta2*r_s(index:end)).^2;


VC = e_c + (-1.*(3./4*pi.*r_s).^(-1/3)/(4*pi.*r_s.^4))*de_c; % e_c + de_c*(n*dr_s/dn) 


end

