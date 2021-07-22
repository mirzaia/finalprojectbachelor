%Script for simulating performance of an automated pressure control 
%structure using a simple dynamic pressure model. 
%Script needs to use external functions calcf.m.
%Written by: ?NS 2009 
%Modified by: DEH 2010 
  
clear all 
close all 
clc 
  
% set up simulation 
% stepsize 
dT = 0.1; 
% time vector 
time=0:dT:200; 
% system parameters 
Va = 128.45; 
Vd = 17.02; 
VaDot = 0; 
betaa = 20000; 
betad = 20000; 
M = 8384; 
% massa jenis mud
rhoa = 0.0140; 
rhod = 0.0140; 
g = 9.81; 
theta1= 900; 
theta2= 60000; 
hBit = 3500; 
WOB = 35.8

Atot = 0.0485; %Va/hBit + Vd/hBit 
Apipe = 0.00218; %7" OD 3" ID 
  
%intial conditions 
pc0 = 13; 
pp0 = 208; 
qbit0 = 3000/60000; 
zc0 = 0.5; 
pbit0 = 575; 
pcref0 = 13; 
hmud0 = 3000; 
qbck0 = 235/60000; 
qp0 = 3000/60000; 
hpipe0 = 3000; 
qc0 = qp0+qbck0; 
  
% storage arrays 
pc = zeros(length(time),1); 
pp = zeros(length(time),1); 
qbit = zeros(length(time),1); 
qbck = zeros(length(time),1); 
qp = zeros(length(time),1); 
zc = zeros(length(time),1); 
pbit = zeros(length(time),1); 
pcref = zeros(length(time),1); 
hmud = zeros(length(time),1); 
qc = zeros(length(time),1); 
hpipe = zeros(length(time),1); 
  
pc(1) = pc0; 
pp(1) = pp0; 
qbit(1)= qbit0; 
zc(1) = zc0; 
pbit(1) = pbit0; 
pcref(1) = pcref0; 
hmud(1) = hmud0; 
qbck(1) = qbck0; 
qp(1) = qp0; 
hpipe(1) = hpipe0;
qc(1) = qc0; 
  
%PI controller to keep pbit at a set point 
pbitref = 575; 
Kppcref = 4; 
Kipcref = 10; 
eipcref = pcref0/Kipcref; 
  
%PI controller to keep pc at a set point 
Kp = 0.01; 
Ki = 5e-3; 
ei = zc0/Ki; 
  
%PI controller to keep zc at a set point 
zcref = 0.5; 
Kpbck = 0.01; 
Kibck = 0.005; 
eibck = qbck0/Kibck; 
  
%PI controller to keep qp at a set point 
qpref = 3000/60000; 
Kpqp = 0.1; 
Kiqp = 0.1; 
eiqp = qp0/Kiqp; 
  
%PI controller to keep hpipe at a set point 
hpiperef = 3000; 
Kphpipe = 0.1; 
Kihpipe = 0.5; 
eihpipe = hpipe0/Kihpipe; 
  
%Euler integration 
for i=1:length(time)-1 
%Define current state vectors 
x = [pc(i); pp(i);qbit(i)]; 
hBit = hBit + 0.125
pbitref = pbitref + 0.0185
%Set inputs 
%if time(i)<50 
    
%qpref = 3000/60000; 
%hpiperef = 3000; 
%pbit(i) = 470; 
%elseif time(i)<200 
%Ramp down pump 
% qpref = 0; 
% zcref = 0.15; 
%pbit(i) = 621; 
% % elseif time(i)<300 
%Trip out drill string 
% hpiperef = 2973; 
% % pbit(i) = 470; 
% % elseif time(i)<400 
%Trip in drill string 
% hpiperef = 3000; 
% % pbit(i)= 440; 
% %  
% % elseif time(i)>=400 
%Ramp up pump 
% qpref = 3000/60000; 
% zcref = 0.5; 
% % pbit(i) = 430; 
% %  
%end 
  
%Controllers 
if i>=2 
  
%Cascade controller to keep pbit at set point by changing pcref 
epcref = pbitref - pbit(i - 1); 
eipcref = eipcref + dT*epcref; 
pcref(i) = max(Kppcref*epcref + Kipcref*eipcref,0); 

        
qpref(i) = 0.001*(max(Kppcref*epcref + Kipcref*eipcref,0)); 
% qpref(i) = qpref;

%Controller to keep qp at set point, with switch 
% eqp = qpref - qp(i - 1); 
eqp = qpref(i-1) - qp(i - 1); 
eiqp = eiqp + dT*eqp; 
if abs(epcref) < 5 
qp(i) = Kpqp*eqp + Kiqp*eiqp; 
else 
qp(i) = qp(i - 1); 
end 
  
%Controller to keep hpipe at set point, with switch 
ehpipe = hpiperef - hpipe(i - 1); 
eihpipe = eihpipe + dT*ehpipe; 
if abs(epcref) < 5 
hpipe(i) = Kphpipe*ehpipe + Kihpipe*eihpipe; 
else 
hpipe(i) = hpipe(i-1); 
end 
  
%Simple controller to keep pc at pcref by changing zc 
e = pc(i-1) - pcref(i-1); 
ei = ei + dT*e; 
zc(i) = min(max(Kp * e + Ki * ei,0),100); 
  
%Input reset controller to keep zc at set point by changing qbck 
ebck = zcref - zc(i - 1); 
eibck = eibck + dT*ebck; 
qbck(i) = max( Kpbck*ebck + Kibck*eibck,0); 
  
end 
  
%Calculate right?hand side 
[f, qccalc] = calcf(time(i),x,zc(i),qbck(i),qp(i),hBit,Va,...
    VaDot,Vd,betaa,betad,M,theta1,theta2,rhod,rhoa,g);

%Step states 
x = x + dT*f; 
  
if i>2 
vpipe = (hpipe(i)-hpipe(i-1))/dT; 
else 
vpipe = 0; 
end 
  
qc(i+1) = qccalc + Apipe*vpipe; 
hmud(i+1) = hmud(i)+Apipe*vpipe/Atot*dT + (qp(i)+ qbck(i)-... 
qc(i))/Atot*dT; 
  
%Store results 
  
pc(i+1) = x(1); 
pp(i+1) = x(2); 
qbit(i+1) = x(3); 
  
%Calculate pbit 
pbit(i+1) = pc(i+1)+theta1*qbit(i+1)+rhoa*g*hBit+WOB; 
  
end 
  
%Store last results 
zc(length(time)) = zc(length(time)-1); 
qp(length(time)) = qp(length(time)-1); 
qbck(length(time)) = qbck(length(time)-1); 
pcref(length(time)-1:length(time)) = pcref(length(time)-2); 
hmud(length(time)) = hmud(length(time)-1); 
hpipe(length(time)) = hpipe(length(time)-1); 
qc(length(time)) = qc(length(time)-1); 
  
%% Plot results
set(0,'defaultaxesfontsize',14); 
set(0,'defaulttextfontsize',14); 
set(0,'DefaultLineLineWidth',1.5); 
set(0,'DefaultFigureColor','none'); 
legFontSize = 16; 
scrsz = get(0,'ScreenSize'); 
  
figure
plot(time,pc,time,pcref,time,pbit-400) 
h(1)=legend('$p_c$','$p_{c,s}$','$p_{bit}$','Location','Best'); 
xlabel('Time [s]'); 
ylabel('Pressure [bar]') 
set(h,'Interpreter','latex','FontSize',legFontSize) 
set(gcf, 'PaperPositionMode', 'auto') 
% Use screen size 
print -djpeg cascadecontrol 
  
figure
plot(time,zc*100) 
h(1)=legend('$z_c$'); 
%axis([0 200 0 100]); 
ylabel('Valve Position [%]'); 
xlabel('Time [s]'); 
set(h,'Interpreter','latex','FontSize',legFontSize) 
set(gcf, 'PaperPositionMode', 'auto') % Use screen size 
print -djpeg cascadecontrol2 
  
figure 
plot(time,qp*60000,time,qbit*60000) 
h(1)=legend('$q_p$','$q_{bit}$','Location','Best'); 
%axis([0 600 0 2500]); 
ylabel('Volume flow [liter/min]'); 
xlabel('Time [s]'); 
set(h,'Interpreter','latex','FontSize',legFontSize) 
set(gcf, 'PaperPositionMode', 'auto') % Use screen size 
print -djpeg cascadecontrol3 
  
figure
plot(time,qbck*60000,time,qc*60000) 
h(1)=legend('q_{bck}','q_c','Location','Best'); 
