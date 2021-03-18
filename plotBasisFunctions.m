%% Method of Moments: Plotting Basis Functions
% Austin Rothschild



%% 5 
%Write a function plotBF.m that accepts the structures Mesh and BF and
%produces a visualization of the mesh?s nodes and segments. Color code the
%plus and minus unit vectors for verification that they are being
%calculated properly.

function [] = plotBasisFunctions(Mesh, BF)


plotMesh(Mesh)
L = BF.Line_L; % length (m)
Z = linspace(-L/2,L/2,Mesh.Nnodes); %Z-axis orientation
X = linspace(0,1,Mesh.Nnodes);
Y = 0*Z;

%% Plot minus basis functions
%These basis functions have a positive linear slope of 1
%figure(1)
%{
for i = 1:(floor(BF.Nbf/2))
    hold on
    plot3(X, Y,  X./(Mesh.P(3,i+1)-Mesh.P(3,i))+ Mesh.P(3,i), 'r')
    hold off
end

for i = (floor(BF.Nbf/2)):BF.Nbf
    hold on 
    plot3(X, Y,  X./(Mesh.P(3,i+1)-Mesh.P(3,i))+ Mesh.P(3,i), 'r') 
    hold off
end

%% Plot plus basis functions
for i = 2:floor(BF.Nbf/2)
    hold on
    plot3(X, Y,  -X./(Mesh.P(3,i+1)-Mesh.P(3,i))+ Mesh.P(3,i+1), 'b')
    hold off
end

for i = floor(BF.Nbf/2):BF.Nbf+1
    hold on
    plot3(X, Y,  -X./(Mesh.P(3,i+1)-Mesh.P(3,i))+ Mesh.P(3,i+1), 'b')
    hold off
end

xlim([-2 2])

axis square
drawnow()
%}
end