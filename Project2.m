clc; clear all; close all;
M1 = 2.5;
alt = 5000; %meters
p = 54000; %54kpa
t = 225; %k
gamma = 1.4;

theta2 = 12*(pi/180); %input degrees
theta3 = -3*(pi/180); %inputdegrees

beta1 = collar(M1,theta2, gamma)*(180/pi)

function shock_angle = collar(M,theta,gamma)
A = M^2-1;
B = ((gamma+1)/2)*M^4*tan(theta);
C = (1+((gamma+1)/2)*M^2)*tan(theta);
x(1) = sqrt(M^2-1);
n=2;
while x(n-1)~= x(n)
x = sqrt(A-(B/(x+C)));
shock_angle = acot(x);
n=n+1;
end
end