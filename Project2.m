clc; clear all; close all;

M1_design = 2.5;
alt_design = 5000; %meters
P1_design = 54000; %54kpa
T1_design = 225; %k
theta_2design = 12*(pi/180); %input degrees
theta_3design = -8*(pi/180); %inputdegrees
[P4_P1design, diff_eff_design, deltaV4_V4prime_design] = diffuser(theta_2design, theta_3design, M1_design, alt_design, P1_design, T1_design);



M1_off = 2.1;
alt_off= 5000;
P1_off = 54000;
T1_off = 225;
theta_2off = 16*(pi/180); %input degrees
theta_3off = -6.5*(pi/180); %inputdegrees
[P4_P1off, diff_eff_off, deltaV4_V4off] = diffuser(theta_2off, theta_3off, M1_off, alt_off, P1_off, T1_off);
error=100*(P4_P1design-P4_P1off)/(P4_P1design) %under 0.5% error


%% Trade Study

P4_P1ts = zeros(1,10);
diff_effts = zeros(1,10);
deltaV4_V4ts = zeros(1,10);
i = 1;
theta2=4;
while theta2 <= 16
  theta3 = -4;
  while theta3 >= -16
      [P4_P1ts(i), diff_effts(i), deltaV4_V4ts(i)] = diffuser((theta2*(pi/180)), (theta3*(pi/180)), M1_design, alt_design, P1_design, T1_design);
      theta3 = theta3 -4;
      i = i+1;
  end
  theta2= theta2 + 4;
end
theta2p=[4:4:16];
theta3p=[-4:-4:-16];
difftheta4=[diff_effts(1),diff_effts(2),diff_effts(3),diff_effts(4)];
difftheta8=[diff_effts(5),diff_effts(6),diff_effts(7),diff_effts(8)]; 
difftheta12=[diff_effts(9),diff_effts(10),diff_effts(11),diff_effts(12)];
difftheta16=[diff_effts(13),diff_effts(14),diff_effts(15),diff_effts(16)];

P4_P1theta4=[P4_P1ts(1),P4_P1ts(2),P4_P1ts(3),P4_P1ts(4)];
P4_P1theta8=[P4_P1ts(5),P4_P1ts(6),P4_P1ts(7),P4_P1ts(8)]; 
P4_P1theta12=[P4_P1ts(9),P4_P1ts(10),P4_P1ts(11),P4_P1ts(12)];
P4_P1theta16=[P4_P1ts(13),P4_P1ts(14),P4_P1ts(15),P4_P1ts(16)];

deltaVtheta4=[deltaV4_V4ts(1),deltaV4_V4ts(2),deltaV4_V4ts(3),deltaV4_V4ts(4)];
deltaVtheta8=[deltaV4_V4ts(5),deltaV4_V4ts(6),deltaV4_V4ts(7),deltaV4_V4ts(8)]; 
deltaVtheta12=[deltaV4_V4ts(9),deltaV4_V4ts(10),deltaV4_V4ts(11),deltaV4_V4ts(12)];
deltaVtheta16=[deltaV4_V4ts(13),deltaV4_V4ts(14),deltaV4_V4ts(15),deltaV4_V4ts(16)];

figure(1)
hold on;
plot(theta3p,difftheta4);
plot(theta3p,difftheta8);
plot(theta3p,difftheta12);
plot(theta3p,difftheta16);
xlabel('theta3 (degrees)');
ylabel('diffuser efficiency');
legend('theta2 = 4','theta2 = 8','theta2 = 12','theta2 = 16');

figure(2)
hold on;
plot(theta3p,P4_P1theta4);
plot(theta3p,P4_P1theta8);
plot(theta3p,P4_P1theta12);
plot(theta3p,P4_P1theta16);
xlabel('theta3 (degrees)');
ylabel('P_4/P_1');
legend('theta2 = 4','theta2 = 8','theta2 = 12','theta2 = 16');

figure(3)
hold on;
plot(theta3p,deltaVtheta4);
plot(theta3p,deltaVtheta8);
plot(theta3p,deltaVtheta12);
plot(theta3p,deltaVtheta16);
xlabel('theta3 (degrees)');
ylabel('deltaV');
legend('theta2 = 4','theta2 = 8','theta2 = 12','theta2 = 16');

%% Thee Function
function [P4_P1, diff_eff, deltaV4_V4prime] = diffuser(theta_2,theta_3, M1,alt, P1,T1)
global y
y = 1.3;
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

diff_eff=100*(P04_P01)+(P04prime_P01)/2*(P0nshock);
deltaV4_V4prime = 100*(2*abs(V4-V4prime))/(V4+V4prime);

function Beta = collar(M,theta)
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

P2_P1 = ((2*y*(M1^2)*(sin(B)^2))/(y+1)) - ((y-1)/(y+1));
p2_p1 = ((y+1)*(M1^2)*(sin(B)^2))/((y-1)*(M1^2)*(sin(B)^2) +2);
T2_T1 = ((1+((y-1)/2)*(M1^2)*(sin(B)^2)*(((2*y)/y-1)*(M1^2)*(sin(B)^2) -1)))/((((y+1)^2)/(2*(y-1)))*(M1^2)*(sin(B)^2));
P02_P01 = (((((y+1)/2)*(M1^2)*(sin(B)^2))/(1+((y-1)/2)*(M1^2)*(sin(B)^2)))^(y*(y-1)))*((1/((2*y/(y+1))*(M1^2)*(sin(B)^2) - ((y-1)/(y+1))))^(1/(y-1)));
M2 = sqrt(((1+((y-1)/2)*M1^2)/((y*M1^2*sin(B)^2)-((y-1)/2)))+((M1^2*cos(B)^2)/(1+((y-1)/2)*M1^2*sin(B)^2)));
end

function F = f(Y)

theta_4 = atan((((Y/P3_P1)-1)/(a_3-(Y/P3_P1)))*sqrt((b_3-(Y/P3_P1))/(c+(Y/P3_P1)))); 
theta_4prime = atan((((Y/P2_P1)-1)/(a_2-(Y/P2_P1)))*sqrt((b_2-(Y/P2_P1))/(c+(Y/P2_P1))));
phi_4 = theta_3 + theta_4;
phi_4prime = theta_2 - theta_4prime;
F = phi_4 - phi_4prime;
end
end
