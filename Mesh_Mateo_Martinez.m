clear;
close all;

% 1. import a geometry from pdeModeler or
load('MyModel.mat');

% 2. define fegeometry and plot it
model = fegeometry(decsg(gd,sf,ns));
figure(1);
pdegplot(model, EdgeLabels="on", VertexLabels="on");

%%
model1 = generateMesh(model);
figure(2);
pdemesh(model1)

%% thinner mesh
figure(3);
model3 = generateMesh(model,Hmax=0.01);
pdemesh(model3)

%% mesh refinement
figure(4);
model4 = generateMesh(model,Hedge={[33 34 35 36 45 47 49 51 37 39 41 43],0.02}, ...
                           Hvertex={[28 32 11 13 23 25 20 16 19 17 26 22 14 10 29 31],0.01});

pdemesh(model4)

%%
save("model4.mat")



