clc; clear all; close all;
global y
M1 = 2;
alt = 5000; %meters
p_s = 54000; %54kpa
t_s = 225; %k
y = 1.3;

theta2 = 9*(pi/180); %input degrees
theta3 = -3*(pi/180); %inputdegrees

beta1 = collar(M1,theta2)*(180/pi)

[mach, T, P, rho, downstream_mach, P0, P1] = flownormalshock(y, M1)
diff_eff=(P04/P01)+(P04prime/P01)/2*(p02/p01)

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

function [P2_P1,p2_p1,T2_T1,P02_P01,M2] = obliqueshock(M1)
global y
P2_P1 = ((2*y*(M1^2)*(sin(B)^2))/(y+1)) - ((y-1)/(y+1));
p2_p1 = ((y+1)*(M1^2)*(sin(B)^2))/((y-1)*(M1^2)*(sin(B)^2) +2);
T2_T1 = ((1+((y-1)/2)*(M1^2)*(sin(B)^2)*(((2*y)/y-1)*(M1^2)*(sin(B)^2) -1)))/((((y+1)^2)/(2*(y-1)))*(M1^2)*(sin(B)^2));
P02_P01 = (((((y+1)/2)*(M1^2)*(sin(B)^2))/(1+((y-1)/2)*(M1^2)*(sin(B)^2)))^(y*(y-1)))*((1/((2*y/(y+1))*(M1^2)*(sin(B)^2) - ((y-1)/(y+1))))^(1/(y-1)));
M2 = sqrt(((1+((y-1)/2)*M1^2)/((y*M1^2*sin(B)^2)-((y-1)/2)))+((M1^2*cos(B)^2)/(1+((y-1)/2)*M1^2*sin(B)^2)));
end
