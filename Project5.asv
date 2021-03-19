clc; clear all ; close all;
%% part 2
global dt
y = 1.4;
p_r = 0.71;
deltan = 0.005;
eta = 0:deltan:20;
dt = deltan;

%% Case A
guessPA = 0;
guessNA = 0.57;
tw_tinf = 2.7;
m_inf = 3;
CaseaA = code(guessPA, guessNA, y, p_r,m_inf, deltan, tw_tinf);
figure(1)
plot(CaseaA(:,2),eta');
hold on;
plot(CaseaA(:,3),eta');
ylabel('eta');title('Case A: m and g vs eta');legend('m','g');

figure(2)
plot(CaseaA(:,2),(sqrt(2).*eta'));
hold on;
plot(CaseaA(:,3),(sqrt(2).*eta'));
ylabel('eta');title('Case A: m and g vs sqrt(2)eta');legend('m','g');

%% Case B
tw_tinf = 1;
m_inf = 3;
guessPB = 0.625;
guessNB = 0.455;
CaseaB = code(guessPB, guessNB, y, p_r,m_inf, deltan, tw_tinf);
figure(3)
plot(CaseaB(:,2),eta');
hold on;
plot(CaseaB(:,3),eta');
ylabel('eta');title('Case B: m and g vs eta');legend('m','g');

figure(4)
plot(CaseaB(:,2),(sqrt(2).*eta'));
hold on;
plot(CaseaB(:,3),(sqrt(2).*eta'));
ylabel('eta');title('Case B: m and g vs sqrt(2)eta');legend('m','g');

%% Case C
tw_tinf = 1/4;
m_inf = 3;
guessPC = 0.613;
guessNC = 0.304;
CaseaC = code(guessPC, guessNC, y, p_r,m_inf, deltan, tw_tinf);
figure(5)
plot(CaseaC(:,2),eta');
hold on;
plot(CaseaC(:,3),eta');
ylabel('eta');title('Case C: m and g vs eta');legend('m','g');

figure(6)
plot(CaseaC(:,2),(sqrt(2).*eta'));
hold on;
plot(CaseaC(:,3),(sqrt(2).*eta'));
ylabel('eta');title('Case C: m and g vs sqrt(2)eta');legend('m','g');

%% Case D
tw_tinf = 0.25;
m_infd = 2;
guessPD = 0.412;
guessND = 0.319;
CaseaD = code(guessPD, guessND, y, p_r,m_infd, deltan, tw_tinf);
figure(7)
plot(CaseaD(:,2),eta');
hold on;
plot(CaseaD(:,3),eta');
ylabel('eta');title('Case D: m and g vs eta');legend('m','g');

figure(8)
plot(CaseaD(:,2),(sqrt(2).*eta'));
hold on;
plot(CaseaD(:,3),(sqrt(2).*eta'));
ylabel('eta');title('Case D: m and g vs sqrt(2)eta');legend('m','g');


%% function
function  solution_array = code(guessP, guessN, y, p_r, minf, deltan, tw_tinf)
global m_inf y p_r 
m_inf=minf;
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
