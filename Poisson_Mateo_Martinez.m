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

figure(4);
generateMesh(model,Hedge={[33 34 35 36 45 47 49 51 37 39 41 43],0.02}, ...
                           Hvertex={[28 32 11 13 23 25 20 16 19 17 26 22 14 10 29 31],0.01});
pdemesh(model)
axis equal


%% Poisson equation with Dirichlet BCs

applyBoundaryCondition(model, 'dirichlet', 'Edge', 33:36, 'u', 0);
applyBoundaryCondition(model, 'dirichlet', 'Edge', [37 39 41 43], 'u', 0.2);
applyBoundaryCondition(model, 'dirichlet', 'Edge', [22 25 19 16 31 28 13 10], 'u', 0);
applyBoundaryCondition(model, 'dirichlet', 'Edge', [45 47 51 49], 'u', 0);

f = @(location,state) 100*exp(-5*((10)*(location.x + 0.894).^2 + (location.y).^2)) + ...
    100*exp(-5*(10*(location.x - 0.799).^2 + (location.y).^2)) + 100*exp(-5*((location.x).^2 + 5*(location.y - 0.923).^2)) ...
    + 100*exp(-5*((location.x).^2 + 5*(location.y+0.836).^2)) + 200*exp(-200*((location.x + 0.023).^2 + (location.y-0.067).^2));
specifyCoefficients(model,"m",0,"d",0,"c",1,"a",0,"f",f); %see help


%%

results = solvepde(model);
u = results.NodalSolution;
pdeplot(model,"XYData",u, "ZData",u);
hold on;
pdemesh(model);
title("Numerical Solution");
xlabel("x")
ylabel("y")
hold off;


