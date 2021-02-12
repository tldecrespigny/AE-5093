clc; clear all; close all;
global  g R a1 expo Tsl T11 rho0 Vflow Pexit Aexit Me Te m0 mp m_dot;
%% Initial Conditions
g = 9.80665;    % m/s^2
R = 287;        % J/(kgK)
a1 = -.0065;    % K/m
expo = -g/(a1*R);  % n/a
Tsl = 288.16;   % K
T11 = 216.66;   % K
rho0 = 1.225;   % kg/m^3
Hmax = 100000; %desired final height of 100km
deltaV = sqrt(2*g*Hmax);
psl = 101320; %pa
a = -.00065; %k/m
tStep = 1;
tend = 100;
tspan=0:tStep:tend;
tbStep = 1;
tbend = 25;
gamma = 1.4;
alt = zeros(1,tbend/tbStep);
V0 = 0;
%% Bread
Me = 4;
Aexit = (.375^2)* pi; %meters
T0 = 3800; %kelvin
gamma = 1.4;
Pexit = 101500; %101.5 kpa
Mpay = 300;
[mach, ToverT0, PoverP0, RhooverRho0, AoverAstar] = flowisentropic(gamma, Me);
combustion_stag_pressure = Pexit/PoverP0;
Te = ToverT0*T0;
m0=Mpay/0.05; %from a payload mass fraction of 5%
Ve = Me*sqrt(gamma*R*Te); 
Vexit = Me*sqrt(gamma*R*Te);
MpoverM0 = 1-exp(-deltaV/Vexit);
mp = MpoverM0*m0;
Mstructural = m0-mp-Mpay;
Mf = Mstructural+Mpay;
At = Aexit/AoverAstar;
m_dot = (Me*sqrt(gamma*R*Te))*Aexit;
Pexit = 101500; %101.5 kpa sea level
m0 = Mstructural+Mpay;
initial_conditions = [V0,0,m_dot];
Vflow = Me*sqrt(1.4*R*Te);
%% Butter
[t1,r1] = ode45(@f_x,tspan,initial_conditions);
plot(t1,r1)
legend(["V_dot" "x_dot" "mdot"])

%% ODE

function x_t = f_x(t,f)
global g Vflow Aexit Pexit m0 m_dot mp;
    V = f(1);
    rho = dens(f(2));
    dm = f(3);
    C_d = .3;
    mp = mp - dm;
    m = mp+m0;
    T = m_dot*Vflow + (Pexit-press(f(2)))*Aexit;
    V_dot = -g-(.5*rho*C_d*Aexit*(V^2))/m + (T/m);
    if mp <= 0
        dm = 0;
        m_dot = 0;
        mp = 0;
        Pexit = 0;
    end
    x_dot = V;
    x_t = [V_dot; x_dot; m_dot;];
end

%% Atmospheric Parameters
function T = temp(h)    % K
global a1 Tsl 
    T = Tsl-a1.*h;
end
function P = press(h)   % Pa
global g R expo Tsl T11
    if h <= 11000
        P = 101320.*(temp(h)./Tsl).^expo;
    else
        P = 22346.*exp(-(g.*(h-11000))./(R.*T11));
    end
end
function rho = dens(h)  % kg/m^3
global g R a1
    Ts = 288.16;
    Rhos = 1.225;   
    temp = Ts + a1*h;
    rho = Rhos*((temp/Ts)^(-1-(g/(a1*R))));
end


%part inout Me Ae To Hdesign (sea level) 

%take in part 4 functions and calc part 5
