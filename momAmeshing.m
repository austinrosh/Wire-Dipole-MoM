%% Method of Moments: Meshing
% Austin Rothschild

%

clc
clear all
close all

%% Parameterize Curve

%this is to make a straight line
Nnodes = 81; %Number of Nodes
Nsegments = Nnodes-1;
tag = 'dipole'

%0.48*lambda = 0.0576
L = 0.0576; % resonant length of a lambda/2 dipole antenna at 2.5 GHz, a little longer than calculated to account for finite wire radius
Z = linspace(-L/2,L/2,Nnodes); %Z-axis orientation
X = 0*Z;
Y = 0*Z;
            


%% Make 3xNnodes array p where the rows are the x,y,z coordinates of each node


%This is valid for a wire dipole antenna
P = zeros(3,Nnodes);
for i = 1:Nnodes
    P(3,i) = Z(i); 
end

%% Make the 2xNsegment array l containing the nodes on either side of the nth segment

l = zeros(2, Nsegments);
for i=1:20
    l(1,i) = i;
    l(2,i) = i+1;
end


%% Make an array "S" 1xNsegments containing the length of each segment

S = zeros(1,Nsegments);

for i=1:Nsegments
    S(i) = norm(P(3,i+1)-P(3,i));
end

%% Make structure "Mesh" 

Mesh.P = P;
Mesh.l = l;
Mesh.S = S;
Mesh.Nnodes = Nnodes;
Mesh.Nsegments = Nnodes-1;
Mesh.Line_L = L;

%% Make a function "plotMesh(Mesh)" that visualizes the mesh

save(['Mesh_',tag, '.mat'],'Mesh')
%plotMesh(Mesh)


