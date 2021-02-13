function rho = dens(h)
%H in Meters
% assumes temperature lapse rate is zero
R = 8.3144598;
M = .0289644;
g = 9.80665;
    if(h<11000)
        rho0 = 1.225;
        T0 = 288.15;
        h0=0;
    elseif(h<20000)
        rho0 = .36391;
        T0 = 216.65;
        h0=11000;
    elseif(h<32000)
        rho0 = .08803;
        T0 = 216.65;
        h0=20000;
    elseif(h<47000)
        rho0 = .01322;
        T0 = 228.65;
        h0=32000;
    elseif(h<51000)
        rho0 = .00143;
        T0 = 270.65;
        h0=47000;
    elseif(h<71000)
        rho0 = .00086;
        T0 = 270.65;
        h0=51000;
    else
        rho0 = .000064;
        T0 = 214.65;
        h0=71000; 
    end
rho = rho0*exp((-g*M*(h-h0))/(R*T0));
end


function P = press(h)
%H in Meters
% assumes temperature lapse rate is zero
R = 8.3144598;
M = .0289644;
g = 9.80665;
    if(h<11000)
        P0 = 101325;
        T0 = 288.15;
        h0=0;
    elseif(h<20000)
        P0 = 22632.10;
        T0 = 216.65;
        h0=11000;
    elseif(h<32000)
        P0 = 5474.89;
        T0 = 216.65;
        h0=20000;
    elseif(h<47000)
        P0 = 868.02;
        T0 = 228.65;
        h0=32000;
    elseif(h<51000)
        P0 = 110.91;
        T0 = 270.65;
        h0=47000;
    elseif(h<71000)
        P0 = 66.94;
        T0 = 270.65;
        h0=51000;
    else
        P0 = 3.96;
        T0 = 214.65;
        h0=71000; 
    end
P = P0*exp((-g*M*(h-h0))/(R*T0));
end