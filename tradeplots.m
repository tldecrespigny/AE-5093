close all; clear all; clc;
mach = [1,1.5,2,3,4,5];
alt = [71065,79903,100160,114450,132950,0];
plot(mach,alt)
xlabel('Exit Mach')
ylabel('Altitude (Meters)')

frac = [0.02,0.03,0.04,0.05,0.06];
alt = [100160,59329,41877,44399,19775];
plot(frac,alt)
xlabel('Payload Mass Fraction')
ylabel('Altitude (Meters)')

ratio = [1.2,1.25,1.3,1.35,1.4];
alt = [103680,97686,100160,89041,89531];
plot(ratio,alt)
xlabel('Specific Heat Ratio')
ylabel('Altitude (Meters)')

temp = [2800,3050,3300,3550,3800];
alt = [92437,96649,91574,86999,100810];
plot(temp,alt)
xlabel('Combustion Chamber Stagnation Temperature')
ylabel('Altitude (Meters)')

Hdesign = [0,5000,10000,20000,40000];
alt = [100810,91001,86255,81495,80541];
plot(Hdesign,alt)
xlabel('Nozzle Design Altitude')
ylabel('Altitude (Meters)')