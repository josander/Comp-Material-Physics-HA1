%% Nice plots
% Task 1: Plot the found wave function

clf
clc

load Task1.mat

x = linspace(0,5,1000);
Psi = @(r) -1*exp(-alpha(1)*r.^2).*C(1) + exp(-alpha(2)*r.^2).*C(2) + ...
    exp(-alpha(3)*r.^2).*C(3)+ exp(-alpha(4)*r.^2).*C(4);

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);
plot(x, Psi(x)/Psi(x(1)), 'LineWidth', 1)
hold on

plotTickLatex2D

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave function [-]','Interpreter','latex', 'fontsize', 12);    
title('Ground state wave function in helium','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

print(gcf,'-depsc2','task1.eps')

%% Task 2: Plot the Hartree potentials

clf
clc

load Task2.mat

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);

V = @(r) 1./r - (1 + 1./r) .* exp(-2.*r);
Vsh = Y(2:end)'./x(2:end) + 1/rMax;

plot(x(2:end), V(x(2:end)),'LineWidth', 1);
hold on
plot( x(2:end), Vsh,'--', 'LineWidth', 1, 'Color', 'red')
plotTickLatex2D

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Electron potential [$E_h$]','Interpreter','latex', 'fontsize', 12);    

title('Electron potential in hydrogen','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

l = legend('Analytic Hartree potential $V_H$','Numerical single Hartree potential $V_{sH}$');
set(l,'Interpreter','latex')

print(gcf,'-depsc2','task2.eps')

%% Task 3: Plot the gound state wave function

clf
clc

load Task3.mat 

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);
plot(x, Psi(x)/Psi(x(1)), 'LineWidth', 1)
hold on
plot( x(2:10:end), waveFuncTask3(2:10:end)/waveFuncTask3(2),'--', 'color', 'red', 'LineWidth', 1)
axis([0 10 0 1]);

plotTickLatex2D

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave functions [-]','Interpreter','latex', 'fontsize', 12);    
title('Ground state wave function in hydrogen','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

l = legend('Analytic wave function $\Psi(r)/\Psi(0)$','Wave function, Task 3 $\Psi_3(r)/\Psi_3(0)$');
set(l,'Interpreter','latex')

print(gcf,'-depsc2','task3.eps')

%% Task 4: Plot the energies with respect to rMax
clf
clc

load Task4rMax2.mat

set(gcf,'renderer','painters','PaperUnits','centimeters','PaperPosition',[0 0 12 7]);

plot(RMax,Energy,'r-', 'LineWidth', 1);
axis([2 15 -79 -70])
X = xlabel('Cutoff radius $r_{max}$ [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Ground state energy [$eV$]','Interpreter','latex', 'fontsize', 12);    
title('Convergence with respect to cutoff radius','Interpreter','latex', 'fontsize', 14);
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);


print(gcf,'-depsc2','convRMax.eps') 

%% Task 4: Plot the different energies with respect to the number of gridpoints

clf
clc

load Task4gridPoints.mat

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);

plot(gridSize,Energy,'r-', 'LineWidth', 1);
axis([0 3000 -78 -77])
plotTickLatex2D

X = xlabel('Grid points  N [-]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Ground state energy [$eV$]','Interpreter','latex', 'fontsize', 12);    

title('Convergence with respect to number of grid points','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);


print(gcf,'-depsc2','convGrid.eps')

%% Task 4: Plot the wave functions

clf


load Task4.mat

psi = psi_r./x;

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);

plot(x(2:end),psi(2:end)/psi(2), 'LineWidth', 1)
hold on
plot(x(2:end), waveFuncTask4(2:end)./waveFuncTask4(2),'--', 'LineWidth', 1, 'Color', 'red');
axis([0 5 0 1]);

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave function [-]','Interpreter','latex', 'fontsize', 12);    

plotTickLatex2D
title('Wave function for Helium','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

l = legend('Wave function, Task 1 $\Psi_1(r)/\Psi_1(0)$','Wave function, Task 4 $\Psi_4(r)/\Psi_4(0)$');
set(l,'Interpreter','latex')


print(gcf,'-depsc2','task4.eps')

%% Task 5: Plot the wave functions

clf
clc

load Task4.mat
load Task5.mat

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);

psi = psi_r./x;

plot(x(2:end),psi(2:end)./psi(2), 'LineWidth', 1)
hold on
plot(x(2:10:end), waveFuncTask4(2:10:end)./waveFuncTask4(2),'r--', 'LineWidth', 1);
hold on
plot(x(2:50:end), waveFuncTask5(2:50:end)./waveFuncTask5(2),'g-.', 'LineWidth', 1);
axis([0 5 0 1]);

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave function [-]','Interpreter','latex', 'fontsize', 12);    

plotTickLatex2D
title('Wave function for Helium','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

l = legend('Wave function, Task 1 $\Psi_1(r)/\Psi_1(0)$','Wave function, Task 4 $\Psi_4(r)/\Psi_4(0)$','Wave function, Task 5 $\Psi_5(r)/\Psi_5(0)$');
set(l,'Interpreter','latex')

print(gcf,'-depsc2','task5.eps')

%% Task 6: Plot the wave functions

clf
clc

load Task4.mat
load Task5.mat
load Task6.mat

set(gcf,'renderer','painters','PaperPosition',[0 0 12 7]);

psi = psi_r./x;

plot(x(2:end),psi(2:end)./psi(2),'-', 'LineWidth', 1)
hold on
plot(x(2:end), waveFuncTask4(2:end)./waveFuncTask4(2),'r--', 'LineWidth', 1);
hold on
plot(x(2:end), waveFuncTask5(2:end)./waveFuncTask5(2),'g-.', 'LineWidth', 1);
hold on
plot(x(2:end), waveFuncTask6(2:end)./waveFuncTask6(2),'m:', 'LineWidth', 1);
axis([0 5 0 1]);

plotTickLatex2D

X = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Normalised wave function [-]','Interpreter','latex', 'fontsize', 12);    

title('Electron potential in hydrogen','Interpreter','latex', 'fontsize', 14);
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
set(X, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

l = legend('Wave function, Task 1 $\Psi_1(r)/\Psi_1(0)$','Wave function, Task 4 $\Psi_4(r)/\Psi_4(0)$','Wave function, Task 5 $\Psi_5(r)/\Psi_5(0)$','Wave function, Task 6 $\Psi_6(r)/\Psi_6(0)$');
set(l,'Interpreter','latex')

print(gcf,'-depsc2','task6.eps')
