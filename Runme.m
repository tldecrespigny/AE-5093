clc; clear all; close all;

%Initial Conditions
mp = 200; %propellant mass kg
mpl = 300; 
m0 = mp + mpl;
tspan=0:tStep:tend;
tStep = .01;
tend = 100000;
gamma = 1.4;
g = 9.81;
r = 287;
alt = zeros(1,tbend/tbStep);
t_bo = 30;
V0 = 0;

Pexit = 101500; %101.5 kpa sea level

initial_conditions = [V0,m0,0]
%Loop
for tb = 1:tbStep:tbend
   options=odeset('Events',@bevents);
   [t1,r1] = ode45(@f_x,tspan,initial_conditions,options);
    
end


function rho = rho(h)
    Ts = 288.16;
    Rhos = 1.225;   
    g = 9.81;
    a1 = -6.5*(10^-3);
    R = 287;
    temp = Ts + a1*h;
    rho = Rhos*((temp/Ts)^(-1-(g/(a1*R))));
end

%ODE
function x_t = f_x(t,f)
dv = f(1);
dm = f(2);
y = 1.4;
rhoh = rho(f(3)); 
mp = mp + dm;
m = mp+mpl;
V = V + dv;

V_dot = -g-(.5*rhoh*C_d*A*(V^2))/m + (T/m);
if t >= t_bo
    m_dot = 0;
else
    m_dot = mdot;
end
x_dot = V;
x_t = [V_dot,m_dot,x_dot];
end

function m_dot = m_dot(mp1,r,gamma)
Mexit = 2;
Aexit = 0.3; %meters
T0 = 3800; %kelvin
%m0=mp1/0.05; %from a payload mass fraction of 5%
gamma = 1.4;
g = 9.81;
r = 287;
Hmax = 100000; %desired final height of 100km
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
end

%part inout Me Ae To Hdesign (sea level) 

%take in part 4 functions and calc part 5
