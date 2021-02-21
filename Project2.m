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



function Beta = collar(M,theta)
global y
A = M^2-1;
B = ((y+1)/2)*M^4*tan(theta);
C = (1+((y+1)/2)*M^2)*tan(theta);
x = zeros(1,20);
x(1) = sqrt(M^2-1);
n=2;
    while x(n-1)~= x(n)
        x(n)= sqrt(A-(B/(x(n-1)+C)));
        n=n+1;
    end
Beta = acot(x);
end