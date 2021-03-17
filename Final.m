clc; clear all; close all;
%% 1
alt = 13500; 
Minf=2.2; 
half_angle = 10; %degrees
Lc=8.65;
Rf=1.5;
Lf=40;
K= 0.405*10^-5;
gamma=1.4;
r=287;

Af_cross=pi*(Rf^2);
Af_surface = 2*pi*Rf*Lf;
Ac_surface = pi*Rf*Lc;
Atotal = Af_surface+Ac_surface;
Aplanform = (Lf*2*Rf);%+(2*Rf*Lc/2); 

pc_p1 =  1.34274065; %from CF calc @M=2.2 and 10 degrees
[T, a, palt, rho] = atmosisa(alt);
qinf = 0.5*rho*((Minf*a)^2);
pc=pc_p1*palt;

%from derivation in HW4 #3
Cd_wave=(pc-palt)/qinf %1.a

rcut = 44.62*((Lf/K)^1.053)*(Minf^1.6);
muinf=(1.458*10^-6*(T^(3/2)))/(T+110.4);
re= (rho*Minf*a*Lf)/muinf;
if rcut<re
    r_use=rcut;
else
    r_use=re;
end

Mc =  2.01121997; %from CF calc @M=2.2 and 10 degrees
M3 = 2.41894; %plugging Mc and half cone angle for turn angle into prandyl myer calc

Cd_parasitic=0.455/((log10(r_use)^2.58)*((1+0.144*M3^2)^0.65)) %1.b
Cd0 = Cd_parasitic*(Af_surface/Af_cross)

Cd_total = Cd_wave+Cd0
Drag = Cd_total*qinf*Af_cross; 
Thrust=Drag/1000 %Kn, concorde has uses 121kn of thrust for ref at cruise m=2.05
