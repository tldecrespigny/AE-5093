clc; close all; clear all;
%% variable inputs
g = 9.80665;    % m/s^2
R = 287;        % J/(kgK)
Me = 3;
y = 1.3;
R_t = 0.5389; %throat radius of project 1 nozzle meters
R_e = 0.7; %exit radius of project 1 nozzle meters

%% part1
[M num mu ] = PMF(y,Me,0,0);
theta_w_max = num/2; %degrees

deltatheta_a_minus1 = 0.3404;
N_f=7;
deltatheta=(theta_w_max-deltatheta_a_minus1)/N_f;
%% part2
mu_t=zeros(1,N_f+1);
nu_t=zeros(1,N_f+1);
M_t = zeros(1,N_f+1);
Kminus = zeros(1,N_f+1);
Kplus = zeros(1,N_f+1);
theta = zeros(1,N_f+1);
theta(1)=deltatheta_a_minus1;
N_t = [1:N_f+1];

for N=1:1:N_f
[M_t(N), nu_t(N), mu_t(N)] = PMF(y,0,theta(N),0);

Kminus(N)=theta(N)+nu_t(N);
Kplus(N)=theta(N)-nu_t(N);
theta(N+1)=theta(N)+deltatheta;
end
[M_t(8), nu_t(8), mu_t(8)] = PMF(y,0,theta(8),0);
Kminus(8)=theta(8)+nu_t(8);
Kplus(8)=theta(8)-nu_t(8);

point=N_t'; %for table implimentation
nu=nu_t';
Kplus=Kplus';
Kminus=Kminus';
theta=theta';
M=M_t';
mu=mu_t';
T = table (point,Kminus,Kplus,theta,nu,M,mu)
%% part3
figure(1);
MinLengthNozzle(y,Me,N+1)

figure(2);
MinLengthNozzle(y,Me,20)

Lmin=xwall(21)*R_t;
Yexit=ywall(21)*R_t;

Lmin_over_2yexit = Lmin/2*Yexit

L= Lmin_over_2yexit*R_e*2

%% part4

%yes the nozzle is practicle for the rocket based on the RL-10 having
%similar dimensions. https://www.rocket.com/space/liquid-engines/rl10-engine#features
% RL10 L/D is 1.5 and ours is 1.3, but we belive that the length given on 
%the website includes the engine as well as nozzle so the L/D would be close to ours.
%However it is not common to use a single engine on a rocket of this size for many reasons.  

%not practicle cuz long as fuck

figure(3);
Me4 = 1.5;
MinLengthNozzle(y,Me4,20)

Lmin4=xwall(21)*R_t;
Yexit4=ywall(21)*R_t;

Lmin_over_2yexit = Lmin4/2*Yexit4
r4=0.549;
L4= Lmin_over_2yexit*r4*2



%% functions
function [ M nu mu ] = PMF(G,M,nu,mu)

Gp=G+1;
Gm=G-1;

% for known M
if M~=0;

    nu = sqrt(Gp/Gm).*atand(sqrt(Gm*(M.^2-1)/Gp))-atand(sqrt(M.^2-1));

    mu = asind(1./M);
    

% for known nu
elseif norm(nu)~=0;
    
    % Find M
        
    %Nu = @(Mg)sqrt(Gp/Gm)*atand(sqrt(Gm*(Mg.^2-1)/Gp))-atand(sqrt(Mg.^2-1))-nu;
    
    for i=1:length(nu(1,:))
        for j = 1:length(nu(:,1))
            M(j,i) = fzero(@(Mg)sqrt(Gp/Gm)*atand(sqrt(Gm*(Mg.^2-1)/Gp))...
                -atand(sqrt(Mg.^2-1))-nu(j,i),[1 100]);
        end
    end

    mu = asind(1./M);
    
    
% for known mu
elseif mu~=0;
    
    M=1./sind(mu);
    
    nu=sqrt(Gp/Gm)*atand(sqrt(Gm*(M.^2-1)/Gp))-atand(sqrt(M.^2-1));
    
end
end

function MinLengthNozzle(G,Me,n)

%{
    Defines geometry for a minimum length nozzle based on a design exit
    mach number for a certain gas, given a finite number (n) of mach waves.
    Based on the information described in Anderson, Modern Compressible
    Flow 3rd Edition (Library of Congress CN: 2002067852).

Input parameters
    G is gamma, the ratio of specific heats (Cp/Cv)
    Me is the design exit mach number
    n is the finite number of expansion waves used in approximation
    
%}
    
%% Initialize datapoint matrices
Km = zeros(n,n);    % K- vlaues (Constant along right running characteristic lines)
Kp = zeros(n,n);    % K- vlaues (Constant along left running characteristic lines)
Theta = zeros(n,n); % Flow angles relative to the horizontal
Mu = zeros(n,n);    % Mach angles
M = zeros(n,n);     % Mach Numbers
x = zeros(n,n);     % x-coordinates
y = zeros(n,n);     % y-coordinates

%% Find NuMax (maximum angle of expansion corner)
[~, B, ~] = PMF(G,Me,0,0);
NuMax = B/2;

%% Define flow of first C+ line
y0 = 1;
x0 = 0;

dT = NuMax/n;
Theta(:,1) = (dT:dT:NuMax);

Nu = Theta;
Km = Theta + Nu;
Kp = Theta - Nu;
[M(:,1) Nu(:,1) Mu(:,1)] = PMF(G,0,Nu(:,1),0);

%% Fill in missing datapoint info along first C+ line
y(1,1) = 0;
x(1,1) = x0 - y0/tand(Theta(1,1)-Mu(1,1));
for i=2:n;
    
    s1 = tand(Theta(i,1)-Mu(i,1));
    s2 = tand((Theta(i-1,1)+Mu(i-1,1)+Theta(i,1)+Mu(i,1))/2);
    x(i,1) = ((y(i-1,1)-x(i-1,1)*s2)-(y0-x0*s1))/(s1-s2);
    y(i,1) = y(i-1) + (x(i,1)-x(i-1,1))*s2;
    
end

%% Find flow properties in characteristic web
for j=2:n;
    for i=1:1+n-j;
        
        Km(i,j) = Km(i+1,j-1);
        
        if i==1;
            
            Theta(i,j) = 0;
            Kp(i,j) = -Km(i,j);
            Nu(i,j) = Km(i,j);
            [M(i,j) Nu(i,j) Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            x(i,j) = x(i+1,j-1) - y(i+1,j-1)/s1;
            y(i,j) = 0;
            
        else
            
            Kp(i,j) = Kp(i-1,j);
            Theta(i,j) = (Km(i,j)+Kp(i,j))/2;
            Nu(i,j) = (Km(i,j)-Kp(i,j))/2;
            [M(i,j) Nu(i,j) Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            s2 = tand((Theta(i-1,j)+Mu(i-1,j)+Theta(i,j)+Mu(i,j))/2);
            x(i,j) = ((y(i-1,j)-x(i-1,j)*s2)-(y(i+1,j-1)-x(i+1,j-1)*s1))/(s1-s2);
            y(i,j) = y(i-1,j) + (x(i,j)-x(i-1,j))*s2;
            
        end
        
    end
end

%% Find wall datapoint info
xwall = zeros(1,n+1);
ywall = zeros(1,n+1);

xwall(1,1) = x0;
ywall(1,1) = y0;

walls = tand(NuMax);
webs = tand(Theta(n,1)+Mu(n,1));

xwall(1,2) = ((y(n,1)-x(n,1)*webs)-(ywall(1,1)-xwall(1,1)*walls))/(walls-webs);
ywall(1,2) = ywall(1,1)+(xwall(1,2)-xwall(1,1))*walls;

for j=3:n+1;
    
    walls = tand((Theta(n-j+3,j-2)+Theta(n-j+2,j-1))/2);
    webs = tand(Theta(n-j+2,j-1)+Mu(n-j+2,j-1));
    xwall(1,j) = ((y(n-j+2,j-1)-x(n-j+2,j-1)*webs)-(ywall(1,j-1)-xwall(1,j-1)*walls))/(walls-webs);
    ywall(1,j) = ywall(1,j-1) + (xwall(1,j)-xwall(1,j-1))*walls;
    
end

%% Provide wall geometry to user and plot
assignin('base','xwall',xwall)
assignin('base','ywall',ywall)

grid=1;
if grid == 1
    
    plot(xwall,ywall,'-')
    axis equal
    axis([0 ceil(xwall(1,length(xwall))) 0 ceil(ywall(1,length(ywall)))])
    hold on
    
    for i=1:n
        plot([0 x(i,1)],[1 y(i,1)])
        plot([x(n+1-i,i) xwall(1,i+1)],[y(n+1-i,i) ywall(1,i+1)])
    end
    
    for i=1:n-1
        plot(x(1:n+1-i,i),y(1:n+1-i,i))
    end

    for c=1:n
        for r=2:n+1-c
            plot([x(c,r) x(c+1,r-1)],[y(c,r) y(c+1,r-1)])
        end
    end

xlabel('Length [x/y0]')
ylabel('Height [y/y0]')

end

end

