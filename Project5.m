clc; clear all ; close all;
global m_inf gamma p_r deltan tw_tinf
gamma = 1.4;
p_r = 0.71;
m_inf = 3;
deltan = 0.05;
tw_tinf = 0;
guess = 0.5;

array = code(guess, gamma, p_r,m_inf, deltan, tw_tinf)

function  array = code(guess, gamma, p_r, m_inf, deltan, tw_tinf, case_num)
guess = 0.5;
initial = [0,0,tw_tinf,0,guess];
array1 = initial;
while m~=u/v_inf
    array1 = rk4_step(f(array1))
    
    n=n+deltan;
end
array2 = intial;
while g~=T/T_inf
    array2 = rk4_step(f(array2))
    
    n=n+deltan;
end
end

function z = f(x)
x(1) = f;
x(2) = m;
x(3) = g;
x(4) = n;
x(5) = p;
fprime = m;
mprime = n;
gprime = p;
nprime = n*(((((1/3)*g^(-1))*p) - f*(g^(1/3))));
pprime = p*((1/3)*(g^(-1))*p - (g^(-1/3))*p_r*f) - p_r*(y-1)*(M_inf^2)*(n^2);
z = [fprime,mprime,gprime,nprime,pprime];
end
  

function x_tplus_dt = rk4_step(x_t)
global dt
k1 = dt* f_x(x_t);
k2 = dt* f_x(x_t + (1/2)*k1);
k3 = dt* f_x(x_t + (1/2)*k2);
k4 = dt* f_x(x_t + k3);

x_tplus_dt = x_t + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
end
