%% Method of Moments: Basis Functions
% Austin Rothschild

clc
clear all
close all

load 'Mesh_dipole.mat'
tag = 'dipole';

%% 1
%Create a list of the indices of the Nbf nodes that are connected to
%exactly two other nodes (i.e., they?re not at either end of the wire).
%Store this list in the 1 × Nbf vector c
c = [];
Nbf = Mesh.Nnodes-2;
for i = 1:Nbf
    c(i) = i + 1;
end
%% 2
%For each basis function described by a center node in c, store the index
%of the ?+? and ?-? segments attached to it. Store these indices in the 1 ×
%Nbf vectors lp and lm.
lp = zeros(1,Nbf);
lm = zeros(1,Nbf);
for i = 1:Nbf
    lp(i) = i+1;
end
for i = 1:Nbf
    lm(i) = i;
end
%% 3
% For each basis function, calculate the unit vectors a+ and a-. 
%Store the x, y, and z components of each unit vector in the rows of the
%3 × Nbf arrays ap and am.
ap = zeros(3, Nbf);
am = zeros(3, Nbf);

% Minus basis functions
for i = 1: Nbf
    am(:,i) = (Mesh.P(:,c(i)) - Mesh.P(:,c(i)-1))/norm((Mesh.P(:,c(i)) - Mesh.P(:,c(i)-1)));
end

% PLus basis functions
for i = 1: Nbf
    ap(:,i) = (Mesh.P(:,c(i)+1) - Mesh.P(:,c(i)))/norm((Mesh.P(:,c(i)+1) - Mesh.P(:,c(i))));
end



%% 4
%. Create a structure BF with the fields c, lp, lm ap, am, and Nbf. Save
%the structure to a file named [?bf ?,tag,?.mat?], where tag is a string
%variable with a short description of the curve shape (e.g., tag =
%?dipole?).


BF.c = c;
BF.lp = lp;
BF.lm = lm;
BF.ap = ap;
BF.am = am;
BF.Nbf = Nbf;
BF.Line_L = Mesh.Line_L;

save(['bf_', tag,'.mat'], 'BF')
plotBasisFunctions(Mesh,BF)

