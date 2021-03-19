clc; clear all ; close all;
%% part 2
global dt
y = 1.4;
p_r = 0.71;
m_inf = 3;
deltan = 0.005;
eta = 0:deltan:20;
tw_tinf = 1;
dt = deltan;
guessP = 0.677;
guessN = 0.477;



CaseaA = code(guessP, guessN, y, p_r,m_inf, deltan, tw_tinf, 1);
% tw_tinf = 1;
% CaseaB = code(guess, y, p_r,m_inf, deltan, tw_tinf, 2);
% tw_tinf = 1/4;
% CaseaC = code(guess, y, p_r,m_inf, deltan, tw_tinf, 2);
% m_inf = 2;
% CaseaD = code(guess, y, p_r,m_inf, deltan, tw_tinf, 2);

figure(1)
plot(CaseaA(:,2),eta');
hold on;
plot(CaseaA(:,3),eta');
%% function
function  solution_array = code(guessP, guessN, y, p_r, m_inf, deltan, tw_tinf, case_num)
global m_inf y p_r 
initial = [0;0;tw_tinf;guessN;guessP];
array = initial;
m=array(2);
g=array(3);

eta = 0:deltan:20;
for i=1:length(eta)
    array = rk4_step(array);
    solution_array(i,:) = array';
end
end



function z = f_x(x)
global m_inf y p_r 
f= x(1);
m= x(2);
g = x(3);
n = x(4);
p = x(5);
fprime = m;
mprime = n;
gprime = p;
nprime = n*(((((1/3)*g^(-1))*p) - f*(g^(1/3))));
pprime = p*((1/3)*(g^(-1))*p - (g^(-1/3))*p_r*f) - p_r*(y-1)*(m_inf^2)*(n^2);
z = [fprime;mprime;gprime;nprime;pprime];
end
  

function x_tplus_dt = rk4_step(x_t)
global dt
k1 = dt.* f_x(x_t);
k2 = dt.* f_x(x_t + (1/2).*k1);
k3 = dt.* f_x(x_t + (1/2).*k2);
k4 = dt.* f_x(x_t + k3);

x_tplus_dt = x_t + (1/6).*k1 + (1/3).*k2 + (1/3).*k3 + (1/6).*k4;
end
