clear;
close all;
clc;

% 1. import a geometry from pdeModeler or
load('model4.mat');

% 2. define fegeometry and plot it
model = createpde();
gm = decsg(gd,sf,ns);
geometryFromEdges(model, gm);

figure 
pdegplot(model,"EdgeLabels","on"); 
axis equal

%%

mesh = generateMesh(model,Hedge={[33 34 35 36 45 47 49 51 37 39 41 43],0.02}, ...
                           Hvertex={[28 32 11 13 23 25 20 16 19 17 26 22 14 10 29 31],0.01});

figure;
pdemesh(model)
axis equal

%%

% Specify PDE coefficients
applyBoundaryCondition(model, 'dirichlet', 'Edge', 33:36, 'u', 50);
applyBoundaryCondition(model, 'dirichlet', 'Edge', [37 39 41 43], 'u', 50);
%edge 4: leave the default homo Neumann

specifyCoefficients(model, 'm', 0, 'd', 1, 'c', 1, 'a', 0, 'f', 1);

% Define initial condition: u(x, y, 0) = sin(pi*x)*sin(pi*y)
setInitialConditions(model, @(location) 100*exp(-5*((10)*(location.x + 0.894).^2 + (location.y).^2)) + ...
    100*exp(-5*(10*(location.x - 0.799).^2 + (location.y).^2)) + 100*exp(-5*((location.x).^2 + 5*(location.y - 0.923).^2)) ...
    + 100*exp(-5*((location.x).^2 + 5*(location.y+0.836).^2)) + 200*exp(-200*((location.x + 0.023).^2 + (location.y-0.067).^2)));

% Time span
Tmax=0.02;
tlist = linspace(0, Tmax, 50);

% Solve PDE
result = solvepde(model, tlist);

% Extract solution at the time step K
K=1;
u = result.NodalSolution(:, K);

% Plot solution at the final time step
figure;
pdeplot(model, 'XYData', u, 'ColorMap', 'jet', 'Mesh', 'off');
xlabel('x'); ylabel('y'); zlabel('Temperature');
title(['Heat Equation Solution at t = ', num2str(tlist(end))]);
colorbar;
axis equal tight;

%%
close all;
figure;
for i = 1:length(tlist)
    % Plot the solution for the current time step
    pdeplot(model, 'XYData', result.NodalSolution(:, i), 'ColorMap', 'jet', 'Mesh', 'off');
    %%plot mesh as well    
    hold on;
    pdemesh(model);
    hold off;
    xlabel('x'); ylabel('y');
    title(['Time: t = ', num2str(tlist(i))]);
    colorbar;
    clim([min(u(:)), max(u(:))]); % Dynamic color range
    axis equal tight;
    pause(0.1); % Adjust speed of animation
end
