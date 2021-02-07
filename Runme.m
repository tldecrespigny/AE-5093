clc 
clear all 
close all

%Initial Conditions
mp = 200;
mpl = 300;
m0 = mp + mpl;
tspan=0:tStep:tend;
tStep = .01;
tend = 100000;
gamma = 1.4;
g = 9.8;
alt = zeros(1,tbend/tbStep);
t_bo = 30;
V0 = 0;
initial_conditions = [V0,m0,0]
%Loop
for tb = 1:tbStep:tbend
   options=odeset('Events',@bevents);
   [t1,r1] = ode45(@f_x,tspan,initial_conditions,options);
    
end


function rho = rho(h)
    Ts = 288.16;
    Rhos = 1.225;   
    g = 9.8;
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


