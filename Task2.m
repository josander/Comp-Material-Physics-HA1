%% Task 2: Get single Hartree potential

clc
clf
clear all

% Cutoff radius
rMax = 50;

% Number of points in the grid
N = 1001; 

% Radial, discetizised points 
x = linspace(0,rMax, N);

% Initialise an array with zeros
Y = zeros(N,1);
Ynew = zeros(N,1);

% The step length between two points
h = rMax/(N-1);

% Single orbital density for the hydrogen atom
a0 = 1; % Bohr radius
densConst = a0^(-3)/(pi);
densFunc = @(r) densConst*exp(-2*r/a0);

diffConst = 2*pi*h^2;
eDens = densFunc(x(2:end));
eDens = eDens.*x(2:end)*diffConst;

nIterations = 5*10^5;

for m = 1:nIterations

    % Loop through the coordinates and calculate new solution
    % Y(0) = Y(N) = 0
    for i = 2:N-1
        Y(i) = eDens(i-1) + 0.5*Y(i+1) + 0.5*Y(i-1);
    end

end

save Task2.mat


%% Plot the Hartree potentials

clf
clc

set(gcf,'renderer','painters','PaperPosition',[0 0 12 8]);

V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
Vsh = Y(2:end)'./x(2:end) + 1/rMax;
plot(x(2:end), V(x(2:end)),'.');
hold on
plot( x(2:end), Vsh,'--', 'MarkerSize', 12, 'Color', 'red')

set(gcf,'renderer','painters','PaperPosition',[0 0 6 3]);
X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Electron potential [Hartree/$a_0$]','Interpreter','latex', 'fontsize', 12);    

title('Electron potential in hydrogen','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

l = legend('Analytic Hartree potential $V_H$','Numerical single Hartree potential $V_{sH}$');
set(l,'Interpreter','latex')

plotTickLatex2D
print(gcf,'-depsc2','task2.eps')

