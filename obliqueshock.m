
function [] = obliqueshock(M1)
global y
P2_P1 = ((2*y*(M1^2)*(sin(B)^2))/(y+1)) - ((y-1)/(y+1));
p2_p1 = ((y+1)*(M1^2)*(sin(B)^2))/((y-1)*(M1^2)*(sin(B)^2) +2);
T2_T1 = ((1+((y-1)/2)*(M1^2)*(sin(B)^2)*(((2*y)/y-1)*(M1^2)*(sin(B)^2) -1)))/((((y+1)^2)/(2*(y-1)))*(M1^2)*(sin(B)^2));
P02_P01 = (((((y+1)/2)*(M1^2)*(sin(B)^2))/(1+((y-1)/2)*(M1^2)*(sin(B)^2)))^(y*(y-1)))*((1/((2*y/(y+1))*(M1^2)*(sin(B)^2) - ((y-1)/(y+1))))^(1/(y-1)))
end