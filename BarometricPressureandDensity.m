function rho = dens(h)
%H in Meters
R = 8.3144598;
M = .0289644;
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


