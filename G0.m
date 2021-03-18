function G = G0(k,r1,r2)
%This function computes the scalar Green's Function given two input vectors
%r1 and r2
    G = (exp(-1j.*k.*norm(r2-r1)))./(4.*pi.*norm(r2-r1));
end

