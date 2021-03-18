f = linspace(10e3, 1e9, 1000);
w = 2.*pi.*f;
e0 = (10e-9)/(36*pi);
u0 = 4*pi*10^-7;
sigma = 5;
gamma_0 = 1j.*w.*sqrt(e0*u0);
gamma = gamma_0.*sqrt(u0.*((1j*sigma)./w));
alpha = real(gamma);
alpha_dB = real(20.*log10(alpha));
semilogx(f,alpha_dB)