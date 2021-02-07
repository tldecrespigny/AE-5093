Mexit = 2;
Aexit = 0.3; %meters
T0 = 3800; %kelvin
%m0=mp1/0.05; %from a payload mass fraction of 5%
gamma = 1.4;
g = 9.81;
r = 287;
Hmax = 100000;
Pexit = 101500; %101.5 kpa
Mpay = 300;

[mach, ToverT0, PoverP0, RhooverRho0, AoverAstar] = flowisentropic(1.4, Mexit);
combustion_stag_pressure = Pexit/PoverP0;
Texit = ToverT0*T0;
m0=Mpay/0.05; %from a payload mass fraction of 5%
Ve = Mexit*sqrt(gamma*r*Texit);
deltaV = sqrt(2*g*Hmax); %initial estimate
Vexit = Mexit*sqrt(gamma*r*Texit);
MpoverM0 = 1-exp(-deltaV/Vexit);
Mp = MpoverM0*m0;
Mstructural = m0-Mp-Mpay;
Mf = Mstructural+Mpay;
At = Aexit/AoverAstar;

