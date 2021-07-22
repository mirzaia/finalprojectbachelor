clear all
clc

Dimension = 4;                          % dimensi diganti sesuai dengan jumlah variabel yang dioptimasi
UB = [200 200 50 100 ];                 % Upper Bounds diganti sesuai dengan constraint fungsi objektif
LB = [10 100 0.001 0.001 ];             % Lower Bounds diganti sesuai dengan constraint fungsi objektif

Npop = 200;     %Number of Population                             
Maxit = 100;    % Maximum Iteration
el = 0.95;      % Elitisim
Pc = 0.8;       % Probabilitas CrossOver
Pm = 0.01;      % Probabilitas Mutasi
%Psilang = 0.8;  %
Nbit = 20;      %


eBangkit = [];
Individu = [];
eIndividu = [];
david = [];
Dadatfit = [];
Datfit = [];
summary = [];
eDadatfit = [];
efitnessmax = [];
eIndividuMax = [];


Bangkit = round(rand(Npop,Nbit*Dimension));
popsize = size(Bangkit,1);
for i = 1:Dimension
batas(i) = UB(i)-LB(i);
end
for i =1:Npop
for j = 1:Dimension
Desimal(i,j)= bi2de(Bangkit(i,((j*Nbit)-(Nbit-1)):(j*Nbit)),'left-msb');
Individu(i,j) = floor((Desimal(i,j)*batas(:,j)-batas(:,j)+LB(:,j)*(2^Nbit-1))/(2^Nbit-1));
end
end
Datfit = [];
variabel = [];
for i = 1:size(Individu,1)
fitness = fobdrilling(Individu(i,:));
Datfit = [Datfit;fitness];
[fitemax,nmax]=max(Datfit);
end
Dadatfit = [];
for generasi=1:Maxit
disp('GA processing')
clear command windows
clear command history
clear memory
if generasi > 1
sort_fit = sortrows(sort,Nbit*Dimension+1);
Individu1 = sort_fit(round((1-el)*Npop+1):Npop,:);
remain = sort_fit(round(el*Npop)+1:Npop,:);
X = Individu1;
M = size(X,1);
sumfitness = sum(Datfit);
for i=1:M
Prob(i) = Datfit(i)/sumfitness;
end
for i=2:M
Prob(i) = Prob(i)+Prob(i-1);
end
for i=1:M
n=rand;
k=1;
for j=1:M-1
if (n>Prob(j))
k=j+1;
end
end
Xparents(i,:) = X(k,:);
end
%=============Crossover========
[M,d] = size(Xparents);
Xcrossed = Xparents;
for i=1:2:M-1
c=rand;
if (c<=Pc)
p=ceil((d-1)*rand);
Xcrossed(i,:) = [Xparents(i,1:p) Xparents(i+1,p+1:d)];
Xcrossed(i+1,:) = [Xparents(i+1,1:p) Xparents(i,p+1:d)];
end
end
if (M/2~=floor(M/2))
c=rand;
if (c<=Pc)
p=ceil((d-1)*rand);
str=ceil((M-1)*rand);
Xcrossed(M,:) = [Xparents(M,1:p) Xparents(str,p+1:d)]; %the first child is chosen
end
end
%============Mutation==============
[M,d] = size(Xcrossed);
Xnew=Xcrossed;
for i=1:M
for j=1:d
p=rand;
if (p<=Pm)
Xnew(i,j)=1-Xcrossed(i,j);
end
end
end
disp('New fitness calculation');
Bangkit = [Xnew(:,1:Nbit*Dimension);remain(:,1:Nbit*Dimension)];
end
eBangkit = [eBangkit; Bangkit];
for i =1:Npop
for j = 1:Dimension;
Desimal(i,j) = bi2de(Bangkit(i,((j*Nbit)-(Nbit-1)):(j*Nbit)),'left-msb');
Individu(i,j) = (Desimal(i,j)*batas(:,j)-batas(:,j)+LB(:,j)*(2^Nbit-1))/(2^Nbit-1);
end
end
Datfit = [];
for i = 1:Npop
fitness = fobdrilling(Individu(i,:));
Datfit = [Datfit;fitness];
[fitemax,nmax] = max(Datfit);
end
Dadatfit = Datfit;
eDadatfit = [eDadatfit;Dadatfit];
eIndividu = [eIndividu;Individu];
[fitnessmax,nmax] = max(eDadatfit);
efitnessmax = [efitnessmax;fitnessmax];
h = plot (efitnessmax);
xlabel('Iteration');
ylabel('Best ROP');
hold on
grid on
refreshdata (h, 'caller')
drawnow;
hold off
BangkitMax = eBangkit(nmax,:);
IndividuMax = eIndividu(nmax,:);
eIndividuMax = [eIndividuMax;IndividuMax];
BangkitMaxlast = BangkitMax;
schedmax = BangkitMax;
clear z;
sort = [Bangkit Dadatfit];
summary = [summary; sort];
david = [david; Dadatfit];
end