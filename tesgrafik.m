% load(NEW_KWA_data.mat,'-mat',y6);
% x = linspace(0,100,10);
% y = linspace(0,100,10);
x1 = Maxitrwa;
x2 = MaxGenerationda;
x3 = generasiga;
y1 = fitnessmaxrwa;
y2 = fitnessmaxfda;
y3 = max_objective_functionga;
% y1 = y6; %L(x)
% y2 = (-27/97750)*x.^3+(4/5)*x; %g(x)
% y3 = (-81/9775)*x + (1726/1955)*x -(108/391); %q(x)
% y4 = (27/97750)*x.^3 + (-162/1955)*x.^2 +6.6864*x -116.2103; %h(x)
% y5 = -1.6*x + 120; %L2(x)
plot(Maxitrwa,fitnessmaxrwa);

% hold on;
% % plot(x,y);
% 
% plot(x2,y2);
% plot(x3,y3);
xlim([0 200]);
ylim([0 15000]);
% plot(x1,y1);
% plot(x,y1,'-r')
% plot(x,y2,'-g')
% plot(x,y3,'-b')
% plot(x,y4,'-m')
% plot(x,y5,'-y')
% plot(x,y6,'-k')
% legend('L(X)','g(x)','q(x)','h(x)','L2(x)','f(x)','location','northeast')
xlabel('Iterasi')
ylabel('Objective Function')
