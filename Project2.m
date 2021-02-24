clc; clear all; close all;
global y theta_3 theta_2 P2_P1 P3_P1 a_2 a_3 b_2 b_3 c
M1 = 2.5;
alt = 5000; %meters
P1 = 54000; %54kpa
T1 = 225; %k
y = 1.3;

theta_2 = 12*(pi/180); %input degrees
theta_3 = -8*(pi/180); %inputdegrees

beta2 = collar(M1,abs(theta_2));
beta3 = collar(M1,abs(theta_3));

[P2_P1,p2_p1,T2_T1,P02_P01,M2] = obliqueshock(M1, beta2);
[P3_P1,p3_p1,T3_T1,P03_P01,M3] = obliqueshock(M1, beta3);

a_2 = 1+(y*(M2^2));
a_3 = 1+(y*(M3^2));
b_2 = (((2*y)/(y+1))*(M2^2)) - ((y-1)/(y+1));
b_3 = (((2*y)/(y+1))*(M3^2)) - ((y-1)/(y+1));
c = (y-1)/(y+1);

%% Newton Rap Son
P4_P1 = fsolve(@f,0);

%% Pt4
Y= P4_P1;
theta_4 = atan((((Y/P3_P1)-1)/(a_3-(Y/P3_P1)))*sqrt((b_3-(Y/P3_P1))/(c+(Y/P3_P1)))); 
theta_4prime = atan((((Y/P2_P1)-1)/(a_2-(Y/P2_P1)))*sqrt((b_2-(Y/P2_P1))/(c+(Y/P2_P1))));
beta4 = collar(M3,abs(theta_4));
beta4prime = collar(M2,abs(theta_4prime));
[P4prime_P2,p4prime_p2,T4prime_T2,P04prime_P02,M4prime] = obliqueshock(M2, beta4prime);
[P4_P3,p4_p3,T4_T3,P04_P03,M4] = obliqueshock(M3, beta4);

P04_P01 = P04_P03*P03_P01;
P04prime_P01 = P04prime_P02*P02_P01;

T4_T1=T4_T3*T3_T1;
T4=T4_T1*T1;
T4prime_T1=T4prime_T2*T2_T1;
T4prime=T4prime_T1*T1;

V4=M4*sqrt(y*287*T4);
V4prime=M4prime*sqrt(y*287*T4prime);
[mach_nshock, T_nshock, P_nshock, rho_nshock, downstream_mach_nshock, P0nshock, P1nshock] = flownormalshock(y, M1);

diff_eff=100*(P04_P01)+(P04prime_P01)/2*(P0nshock)
deltaV4_V4prime = (2*abs(V4-V4prime))/(V4+V4prime)



%% functions
function Beta = collar(M,theta)
global y
A = M^2-1;
B = ((y+1)/2)*M^4*tan(theta);
C = (1+((y+1)/2)*M^2)*tan(theta);
x = zeros(1,20);
x(1) = sqrt(M^2-1);
n=2;
accuracy = 1;

    while accuracy > 0.000001
        x(n)= sqrt(A-(B/(x(n-1)+C)));
        accuracy = abs(x(n-1)- x(n));
        n=n+1;
    end
Beta = acot(x(n-1));
end

function [P2_P1,p2_p1,T2_T1,P02_P01,M2] = obliqueshock(M1,B)
global y
P2_P1 = ((2*y*(M1^2)*(sin(B)^2))/(y+1)) - ((y-1)/(y+1));
p2_p1 = ((y+1)*(M1^2)*(sin(B)^2))/((y-1)*(M1^2)*(sin(B)^2) +2);
T2_T1 = ((1+((y-1)/2)*(M1^2)*(sin(B)^2)*(((2*y)/y-1)*(M1^2)*(sin(B)^2) -1)))/((((y+1)^2)/(2*(y-1)))*(M1^2)*(sin(B)^2));
P02_P01 = (((((y+1)/2)*(M1^2)*(sin(B)^2))/(1+((y-1)/2)*(M1^2)*(sin(B)^2)))^(y*(y-1)))*((1/((2*y/(y+1))*(M1^2)*(sin(B)^2) - ((y-1)/(y+1))))^(1/(y-1)));
M2 = sqrt(((1+((y-1)/2)*M1^2)/((y*M1^2*sin(B)^2)-((y-1)/2)))+((M1^2*cos(B)^2)/(1+((y-1)/2)*M1^2*sin(B)^2)));
end

function F = f(Y)
global theta_3 theta_2 P3_P1 P2_P1 a_2 a_3 b_2 b_3 c
theta_4 = atan((((Y/P3_P1)-1)/(a_3-(Y/P3_P1)))*sqrt((b_3-(Y/P3_P1))/(c+(Y/P3_P1)))); 
theta_4prime = atan((((Y/P2_P1)-1)/(a_2-(Y/P2_P1)))*sqrt((b_2-(Y/P2_P1))/(c+(Y/P2_P1))));
phi_4 = theta_3 + theta_4;
phi_4prime = theta_2 - theta_4prime;
F = phi_4 - phi_4prime;
end