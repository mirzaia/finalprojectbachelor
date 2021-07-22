function ROP=fobdrilling(x)
%optimized variable
WOB=x(1) %weight on bit (tons)
RPM=x(2) %rotary speed (RPM)
qmud=x(3) %flow rate lumpur (L/s)
pc=x(4) %pressure choke (bar)
rho= 1403.805631 %density lumpur (kg/m3)

%Parameter drilling model
R0=5; %formation drillability (m/hr)
pf=470; %formation pressure at bottom well (bar)
tt0= 10; %trip time hours
wdbmax=178.583; %maximum WOB per diameter (tons/m)
wdbt=63/R0; %threshold WOB per m diameter (tons/m)
N=100; %drill string RPM (min-1)
db=0.254; %drill bit diameter (m)
theta1=900; %annulus friction parameter (kg/m4s)
dstring=0.1; %drill string outter diameter (m)
D=3000; %depth well (m)
dp=0.005; %particle diameter (m)
phi=pi; 
rhos=2700; %formation density (kg/m3)
g=9.8; %gravity (m/s2)


%Rate of Penetration Modelling
kc=5e5/R0; %valve flow konstan (m2)
Aa=(db^2-dstring^2)*phi/4; %cross-sectional area of annulus (m2)
An=3*pi*(0.01^2)/4; %nozzle x-sec area (m2)
F=x(1)/db*kc*db^2*phi/4; %force
T=F*dstring/2000; %torque (kNm)
newrho=rho*(1-T/100)+rhos*T/100; %rho friction (kg/m3)
pbh=x(4)+(theta1*x(3)/1000)+newrho*g*D/1e5   %bottom hole pressure (bar)
P=T*2*phi*x(2)/60; %power (kW)
Fj=x(3)^2*rho/An/1e6; %hydraulic jet impact force (N)
ROP=(R0*((x(1)/db-wdbt)/(71.4-wdbt))*((N/60)^0.7)*((Fj/4482)*0.3)*exp(0.01*(pf-pbh)))

K=15*R0; %formation of abrasiveness constant (hours)
va=qmud/100/Aa; %annulus velocity (m/s)
vslip=sqrt(8/9*g*dp*(rhos-newrho)/rho); %slip velocity (m/s)
vT=va-vslip; %transport velocity (m/s)
qs=ROP*phi/4*db^2*1000/3600; %cutting feed rate (1/s)
xc=qs/1000/Aa/vT*100; %cuttings fraction

%cost of drilling
td0=K*(60/ROP)^1.7*(wdbmax-x(1)/db)/(wdbmax-71.433); %bit life-time
cost=D*((1/ROP)+(tt0/(ROP*td0))) %cost drilling
end