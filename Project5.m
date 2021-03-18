clc
clear all 
close all



function z = f(x)
x(1) = f;
x(2) = m;
x(3) = g;
x(4) = n;
x(5) = p;
fprime = m;
mprime = n;
gprime = p;
nprime = n*((((1/2)*(g^(-1/3))*p) - f*(g^(-1/3))));
pprime = p*((1/3)*(g^(-1/3))*p - (g^(-1/3))*p_r*f) - p_r*(y-1)*(M_inf^2)*(n^2);
z = [fprime,mprime,gprime,nprime,pprime];
end
    