%% Task 6

clc
clear all

% Cutoff radius
rMax = 15;

% Number of points
N = 3001; 

% Radial, discetizised points 
x = linspace(10^(-9),rMax, N);

% Coefficients of wave function from task 1
C = -1*[-0.146916049461378, -0.393060020070374, -0.411115799349951, -0.261968242091914];

% Declaration of alpha
alpha = [0.297104, 1.236745, 5.749982, 38.216677];

% Initialise the wave function with solution from task 1
psi_r = (exp(-alpha(1)*x.^2).*C(1) + exp(-alpha(2)*x.^2).*C(2) + ...
    exp(-alpha(3)*x.^2).*C(3)+ exp(-alpha(4)*x.^2).*C(4)).*x;

% Normalise psi_r : 4pi int(r^2U0^2) = 1  
psi_r = psi_r/sqrt(trapz(4*pi.*x.^2.*psi_r.^2));

% Length between two points
h = rMax/(N-1);

% Variables to keep track of the energy difference
energyDiff = 1;
Eold = 0;

% Iterate until the convergence condition; the maximal energy difference
% 10^-5
while energyDiff > 10^(-5) % [eV]

    % Get the single Hartree potential
    Vsh = solveVSH(x, psi_r);
    
    % Get the exchange potential
    Vx = solveVEx(psi_r);
    
    % Get the correlation potential
    Vc = solveVC(psi_r);
    
    % Define the potential
    pot = -2./x+2*Vsh+Vx+Vc;
    
    % Solve the Khon-Sham equation and get the eigenvalues and the
    % eigenvectors
    [A B] = solveKS(pot, x);

    % Get the eigenvalues
    e = (diag(B));

    % Find index of the minimal eigenvalue
    index = min(find(e == min(e)));
    
    % The new radial wave function
    psi_r = A(:,index)';
    
    % Normalise psi_r : 4pi int(r^2U0^2) = 1  
    psi_r = psi_r/sqrt(trapz(4*pi.*x.^2.*psi_r.^2));

    % Get the minimal eigenvalue in Hartree energy
    minEig = e(index);

    % Get energy in eV
    E = 27.211396132*minEig;

    % Calculate the new energy difference
    energyDiff = abs(Eold - E)
    
    % Save the solution
    Eold = E;

end

% Energy in eV
Energy = E

% Energy in Eh
EnergyHartree = minEig;

waveFuncTask6 = A(:,index)'./x;

save Task6.mat

%% Plot the wave functions

clf
clc

load Task4.mat
load Task5.mat
load Task6.mat

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

psi = psi_r./x;

plot(x(2:end),psi(2:end)./psi(2),'-', 'LineWidth', 1)
hold on
plot(x(2:end), waveFuncTask4(2:end)./waveFuncTask4(2),'r--', 'LineWidth', 1);
hold on
plot(x(2:end), waveFuncTask5(2:end)./waveFuncTask5(2),'g-.', 'LineWidth', 1);
hold on
plot(x(2:end), waveFuncTask6(2:end)./waveFuncTask6(2),'b', 'LineWidth', 1);
axis([0 5 0 1]);

plotTickLatex2D

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave function [-]','Interpreter','latex', 'fontsize', 12);    

title('Electron potential in hydrogen','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0]);

l = legend('Wave function, Task 1 $\Psi_1(r)/\Psi_1(0)$','Wave function, Task 4 $\Psi_4(r)/\Psi_4(0)$','Wave function, Task 5 $\Psi_5(r)/\Psi_5(0)$','Wave function, Task 6 $\Psi_6(r)/\Psi_6(0)$');
set(l,'Interpreter','latex')

print(gcf,'-depsc2','task6.eps')
