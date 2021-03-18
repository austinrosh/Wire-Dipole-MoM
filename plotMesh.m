function [] = plotMesh(Mesh)


figure(1)
Nnodes = Mesh.Nnodes;
Nsegments = Mesh.Nsegments;
P = Mesh.P;
l = Mesh.l;
S = Mesh.S;


L = Mesh.Line_L; % % resonant length of a lambda/2 dipole antenna at 2.5 GHz, a little longer than calculated to account for finite wire radius
Z = linspace(-L/2,L/2,Nnodes); %Z-axis orientation
X = 0*Z;
Y = 0*Z;

figure(1)
plot3(X,Y,Z, '-o')
title("Dipole Antenna Meshing")
xlabel('x')
ylabel('y')
zlabel('z')
axis square
drawnow()










end

