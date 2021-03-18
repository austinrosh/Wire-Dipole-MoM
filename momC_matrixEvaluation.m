h%% Method of Moments: Green's Function Evaluations
% Austin Rothschild


clc
clear
close all

%Run Scripts
run('momAmeshing.m')
run('momB_basisfunctions.m')

%Load in Basis Function/Mesh data
load 'bf_dipole.mat'
load 'Mesh_dipole.mat'

%% Computing Z Matrix and Unknown Currents

Z_ = zeros(BF.Nbf, BF.Nbf);% Initialize storage of impedance matrix
Zin = []; % Initialize storage of input impedance vector
I_m = []; %Feed point current 

a = 6.0000e-10; %Wire radius
Z0 = 377; %Impedance of free space, Ohms

V_ = zeros(BF.Nbf,1); %Inititialize voltage
V_(ceil(BF.Nbf/2)) = 1; %Gap voltage source

% Preprocess bf structure single point integration rule
for i = 1:BF.Nbf
        rp_(i,:) = (Mesh.P(:,BF.lp(i)+1)+Mesh.P(:,BF.lp(i)))/2; %midpoint of plus segement of the i^th basis function
        xp_(i,:) = BF.ap(:,i)/2; %one half of the unit vector in plus segment of i^th basis function
        rm_(i,:) = (Mesh.P(:,BF.lm(i))+Mesh.P(:,BF.lm(i)+1))/2; %midpoint of minus segement of the i^th basis function
        xm_(i,:) = BF.am(:,i)/2; %one half of the unit vector in minus segment of i^th basis function
        dp_(i) = -1/(Mesh.S(i)); %divergence of plus segment of the i^th basis function
        dm_(i) = 1/(Mesh.S(i)); %divergense of plus segment of i^th basis function 
end

rp_ = rp_.';
xp_ = xp_.';

rm_ = rm_.';
xm_ = xm_.';

% cycle through all combinations of basis functions

f = linspace(1.5e9,7.0e9, 250);

for i = 1:length(f)
    
    lambda = 3e8./f(i); %wavelength
    k = 2.*pi./lambda; %spatial wavenumber
    
    % cycle through all combinations of basis functions
    for m = 1:BF.Nbf
        for n = 1:BF.Nbf

            %Greens function from lm+/- to lp+/i
            Gpp = G0(k, rp_(:,m), rp_(:,n));
            Gpm = G0(k, rm_(:,m), rp_(:,n));
            Gmp = G0(k, rp_(:,m), rm_(:,n));
            Gmm = G0(k, rm_(:,m), rm_(:,n));

            %deal with singularities
            %if plus of mth bf is same as plus of nth bf
            if abs((rp_(:,m)- rp_(:,n))) <= eps
               S = Mesh.S(n);
               Gpp = ((S/(2*pi))*(log(S/a +sqrt(1+(S^2)/(a^2))) - sqrt(1+(a^2)/(S^2)) + a/S) - (1j.*k.*S^2)/(4*pi))*1/S^2;
            end

            if abs((rm_(:,m)- rp_(:,n))) <= eps
                S = Mesh.S(n);
                Gpm = ((S/(2*pi))*(log(S/a +sqrt(1+(S^2)/(a^2))) - sqrt(1+(a^2)/(S^2)) + a/S) - (1j.*k*S^2)/(4*pi))*1/S^2;
            end

            if abs((rp_(:,m)- rm_(:,n))) <= eps
                S = Mesh.S(n);
                Gmp = ((S/(2*pi))*(log(S/a +sqrt(1+(S^2)/(a^2))) - sqrt(1+(a^2)/(S^2)) + a/S) - (1j.*k*S^2)/(4*pi))*1/S^2;
            end

            if abs((rm_(:,m)- rm_(:,n))) <= eps
                S = Mesh.S(n);
                Gmm = ((S/(2*pi))*(log(S/a +sqrt(1+(S^2)/(a^2))) - sqrt(1+(a^2)/(S^2)) + a/S) - (1j.*k*S^2)/(4*pi))*1/S^2;
            end

            % Compute scalar potential fields 
            Ppp = Gpp;
            Ppm = -Gpm;
            Pmp = -Gmp;
            Pmm = Gmm;
            
            %Scalar Potential Matrix
            P = Ppp + Ppm + Pmp + Pmm;

            %Compute vector potential terms
            App = dot(BF.ap(:,m),BF.ap(:,n))*Mesh.S(BF.lp(m))*Mesh.S(BF.lp(n))*Gpp/4;
            Apm = dot(BF.am(:,m),BF.ap(:,n))*Mesh.S(BF.lm(m))*Mesh.S(BF.lp(n))*Gpm/4;
            Amp = dot(BF.ap(:,m),BF.am(:,n))*Mesh.S(BF.lp(m))*Mesh.S(BF.lm(n))*Gmp/4;
            Amm = dot(BF.am(:,m),BF.am(:,n))*Mesh.S(BF.lm(m))*Mesh.S(BF.lm(n))*Gmm/4;
            
            %Vector Potential Matrix
            A = App + Apm + Amp + Amm;

            
            %Impedance Matrix
            Z_(m,n) = 1j.*k.*Z0.*(A - (1./k.^2).*P);
        end
    end
    
   
    %Compute Current coefficients from V = ZI Relation
    I_ = Z_\V_;
    
    %Save current distribution @ 2.5 GHz
    if i == 47
        I_m = I_;   
    end
     
    %Compute input impedance over each frequency
    Zin(i) = V_(ceil(BF.Nbf/2))/I_(ceil(BF.Nbf/2)); 
    
    
    %characteristic modes calculation
    X = imag(Z_);
    R = real(Z_);
    [v, u] = eigs(X,R,6,'sm');
    num = diag(u);
    lam1(i) = num(1);
    lam2(i) = num(2);
    lam3(i) = num(3);
    lam4(i) = num(4);
    lam5(i) = num(5);
    lam6(i) = num(6);

end
        
    
    



%% Post Processing

%Current Distribution Plots
figure(2)
plot(real(I_m))
title('Curent Distribution of a Half-Wave Dipole @ 2.5 GHz')
xlabel('Node Number')
ylabel('I(z), Amps')
xlim([1 BF.Nbf])

%Input Impedance Plots
figure(3)
hold on
plot(f,real(Zin))
plot(f, imag(Zin), '-')
Zcheck = f.*0;
Rr = Zcheck;
Rr(1,:) = 73; %set radiation resistance to be 73 ohms at target frequency
plot(f, Zcheck)
plot(f, Rr)
title('Input Impedance for a \lambda/2 dipole @ 2.5 GHz')
legend('Real{Z_{in}}', 'Imag{Z_{in}}', 'Z_{in} = 0', 'Z_{in} = 73', 'Location', 'Best')
xlabel('Frequency, Hz')
ylabel('Zin, \Omega')
xline(2.56e9)
hold off
xlim([1.5e9 7e9])

%index 46 or 47 ~= 2.5 GHz ->check Zin values here

%%Zin(45) 153 BFs = 73.7352 +16.5363i
%Zin(45) 121 BFs = 73.6787 +15.7963i
%Zin(45) 99 BFs = 73.6556 +15.4988i
%Zin(45) 79 BFs = 73.6215 +15.0630i
%Zin(45) 73 BFs = 72.0153 +14.8331i
%Zin(45) 63 BFs = 71.7259 -11.8088i
%Zin(47) 47 BFs = 70.9660 -67.1394i
%Zin(47) 35 BFs = 69.544 -68.6936i
%Zin(47) 19 BFs = 67.3780 -75.1384i

%Convergence of Radiation Resistance @ 2.5 GHz vs. # of Basis Functions
Zin_conv = [67.3780 69.544 70.9660 71.7259 72.0153 73.6215 73.6556 73.6787 73.7352]
NBfs = [19 35 47 63 73 79 99 121 153]

%Plotting Convergence
figure(4)
plot(NBfs, Zin_conv, 'o')
title('Accuracy of Solver vs. Number of Basis Functions')
ylabel('Z_{in} \Omega')
xlabel('# of Basis Functions')
yline(73)
%xline(2.5e9)


%% Characteristic Mode Analysis
% GEP X*In = Lambda*R*In
% GEP is capable of diagnolizng both R and X operators
figure(5)
semilogy(f, abs(lam1(1,:)),f, abs(lam2(1,:)),f, abs(lam3(1,:)),f, abs(lam4(1,:)), f, abs(lam5(1,:)),f, abs(lam6(1,:)),'LineWidth',1.5)
xlim([1.5e9 7e9])
%yticks([10^-2 10^0 10^2 10^4 10^6 10^8 10^10])
%yticklabels({'-20','0','20','40','60','80','100'})
xlabel('Frequency (Hz)')
ylabel('|\lambda| (dB) for select Z-Matrix Modes')
legend('\lambda_1(\it f)','\lambda_2(\it f)','\lambda_3(\it f)','\lambda_4(\it f)','\lambda_5(\it f)','\lambda_6(\it f)','location','best')
title('Characteristic Modes for an 81-Basis Function Z matrix')
grid on