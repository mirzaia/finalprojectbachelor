function profit=fobjco2egr (x)
%optimized variable
minj=x(1) %laju aliran massa injeksi (kg/s)
pinj=x(2) %tekanan injeksi (psia)
tinj=x(3) %temepratur injeksi (C)

%Parameter model injection well
dwell=0.127; %diameter well (m)
g=9.8; %percepatan gravitasi (m/s2)
gc=1; %faktor gravitasi (kg m/N s2)
depth=1700; %kedalaman sumur (m) % DIGANTI
thick=0.0053848; %ketebalan dinding (m)
rough=0.0000675; %kekasaran dinding (m)
trev=38.889; %temperatur reservoir (C) % DIGANTI
tling=31; %temperatur lingkungan (C)

%INJECTION WELL

mjci= 750.159157983176-15.4538661097692*tinj+0.294070293859911*pinj; %massa jenis co2 injection well (kg/m3)
volci= minj/mjci; %Volume flowrate co2 (m3/s)
viscoci= 0.000074101376759184+2.44308946085363E-08*pinj+-1.22646600912301E-06*tinj; %viscositas co2 injection well (kg/ms)
vlci= volci/(3.14*((dwell^2)/4)); %velocity co2 injection well (m/s)
Gammai=0.07275*(1-0.002*((tinj+273)-291)); % Surface Tension
LVNi=vlci*(mjci/(g*Gammai)^0.25);%Liquid Velocity Number CO2
NFRi=(vlci^2)/(g*dwell); %Froude Number Co2 Injection 
    
    % Mentukan Pola Aliran
    lambdai=1;
    Eksi=log(lambdai);
    L1i=exp(-4.62-(3.757*Eksi)-(0.481*Eksi^2)-(0.0207*Eksi^3));
    L2i=exp(1.061-(4.602*Eksi)-(1.609*Eksi^2)-(0.179*Eksi^3)+(0.635*(10^-3)*Eksi^5));
    if ((NFRi<L1i)) || ((NFRi<L2i))
        pola = 1; %1=segregated
    elseif  ((NFRi>L1i) &&(NFRi>L2i))
            pola = 2;%2=distributed
    else
        pola=3; %3=intermittent
    end
    Hlsci=(0.98*(lambdai^0.4846)/(NFRi^0.0868));
    Hlici=(0.84*(lambdai^0.5351)/(NFRi^0.0173));
    Hldci=(1.06*(lambdai^0.5824)/(NFRi^0.0609));
     %Nilai C Aliran
    Csi=(1-lambdai)*log((4.7*LVNi^0.0868)/((lambdai^0.3692)*(NFRi^0.5056)));
    Cii=(1-lambdai)*log((4.7*LVNi^0.1244)/((lambdai^0.3692)*(NFRi^0.5056)));
    Cdi=(1-lambdai)*log((4.7*LVNi^0.1244)/((lambdai^0.3692)*(NFRi^0.5056)));
    if pola == 1
        HLI = Hlsci;
        Ci = Csi;
    elseif pola == 2
        HLI = Hldci;
        Ci = Cdi;
    else HLI = Hlici;
        Ci = Cii;
    end
    
    PSIi=1+Ci*(sin(1.8*-90)-((1/3)*sin(-90)^3));
    HLItetha = HLI * PSIi; 
    yi = lambdai/(HLItetha^2);
    Si = (log(yi))/(-0.0523+3.812*log(yi)-0.8725*(log(yi)^2)+0.01853*(log(yi)^4));
    
    NREi = mjci*vlci*dwell/viscoci;
    fnsi = (2*log10(NREi/(((4.5223*log10(NREi))-3.8215))))^(-2);
    ftpi = exp(Si)*fnsi;
   
    %Pressure Drop
    pdfi = (ftpi*minj*vlci^2)/(2*gc*dwell)
    pdei = (g*mjci)/gc;
    pdtoti = (pdfi+pdei)*0.0000442075*depth*3.28084;
    
     %heat transfer co2 injection well
    Kci=-8.62245548121198E-06*pinj+(0.00165497333443547*tinj)-0.00120940124857462;   %Konduktivitas Thermal CO2 injection well (W/mK)
    Cpci= (22.9307000677034+-0.100097810122466*tinj+-0.00261700285458468*pinj)*1000;%Heat capacity CO2 injection well (kJ/kgK)
    NPrci=Cpci*viscoci/Kci;%Prandlt Number CO2 injection well
    NNuci=0.023*(NREi^0.8)*(NPrci^0.3);%Nuselt Number CO2 injection well
    hci=NNuci*Kci/dwell;%heat transfer coefficient CO2 injection well (W/m2K)
    Rci=thick/(Kci*(3.14 *(dwell^2)/4));%Resistansi konduktivitas thermal CO2 injection well (K/W) 0.005-wallthicknes
    Uci=hci+(1/(Rci*3.14*dwell^2))%All Heat transfer coefficient CO2 injection well (W/m2K)
    Qci=minj*Cpci*(tinj-trev)*35%energy panas CO2 injection well (watt)
    dtci=Qci/(Uci*3.14*dwell*depth*2) %Delta T CO2 injection well (C)
    
%RESERVOIR

tcr=tinj-dtci %temperature co2+oil pada reservoir C
pcr=(pinj+pdtoti)*1.1 %pressure Co2+oil pada reservoir psi
    
%parameter model
frac=0.5582; %fraksi liquid awal
resthick=12.192; %reservoir thickness (m) % DIGANTI
reslength=100; %reservoir length (m)% DIGANTI
Krock=0.02630; %thermal conductivity (W/mK) % DIGANTI
permeabil= 8.48734E-14; %permeability (m2)% DIGANTI
permeabilv= 8.48734E-15; %permeability vertikal (m2) % DIGANTI
poros= 0.089; %porosity batuan % DIGANTI
Lling = 116.6862; % 3316.625; %Luas Lingkaran % DIGANTI

    %pressure drop
    mjcr= 46.8300924806693-1.3210030053444*tcr+0.0865162316754536*pcr%massa jenis phase liquid co2+Oil reservoir (kg/m3)
    volcr= minj/mjcr;%Volume flowrate co2+Oil reservoir (m3/s)
    vlcr=volcr/(3.14*((resthick^2)/4)); %velocity co2+oil reservoir (m/s)
    viscocr =9.82074864020945E-06+tcr*-5.57764956509715E-08+pcr*6.21051761785609E-09;%viscositas phase liquid co2+Oil reservoir (kg/ms)
    NRecr=mjcr*vlcr*resthick/viscocr; %reynold number co2+oil reservoir
    pdcr= -((viscocr*volcr*reslength)/(permeabil*3.14*resthick*reslength)*0.000145038);%pressure drop co2 reservoir psi. thicknes 35.052m

    %heat transfer
    Cpcr=2.37305418766803+tcr*-0.0109693863865791+0.000240054091308146*pcr;  %Heat capacity CO2+oil di reservoir (J/kgK)
    Kcr=0.0174865883625187+-0.000124017589052688*tcr+0.0000167603402817045*pcr; %Thermal conductivity co2+oil di reservoir (W/mK)
    NPrcr=Cpcr*viscocr/Kcr; %PRandlt Number CO2+oil di reservoir 1.73-thermal conductivity batuan
    NNucr=0.023*(NRecr^0.8)*(NPrcr^0.3);%Nuselt Number CO2+oil di reservoir
    Lpcr=(3.14*0.25*(resthick^2)*reslength)^(1/3); %length characteristic di reservoir
    hcr=NNucr*Krock/Lpcr ;%heat transfer coefficient di reservoir 5.678263 W/mK=Kthermal batuan
    Rkvcr=1/(hcr*(3.14*0.25*(resthick^2)));%tahanan konveksi perpindahan panas di rservoir
    Rkdcr=Lpcr/(Krock*3.14*0.25*(resthick^2));%tahanan konduksi perpindahan panas di reservoir
    Qcr=((tcr-48.33)/(Rkvcr+Rkdcr))*100;%energy panas co2+oil yg terbuang di reservoir (J)
    dtcr=Qcr/((Rkvcr+Rkdcr)*Cpcr); %Delta T Co2+oil di rservoir (C)

%PRODUCTION WELL

tcp=tcr-dtcr %temperature co2+oil production well (C)
pcp=pcr+pdcr%pressure co2+oil production well psi

mjcp= 52.8443199090267-1.35123899611828*tcp+0.0844646786671024*pcp; %massa jenis co2 injection well (kg/m3)
volcp= minj/mjcp; %Volume flowrate co2 (m3/s)
viscocp= 9.93335855782682E-06+tcp*-6.23116671510213E-08+pcp*6.27709117124863E-09; %viscositas co2 injection well (kg/ms)
vlcp= volcp/(3.14*((dwell^2)/4)); %velocity co2 injection well (m/s)
Gammap=0.07275*(1-0.002*((tinj+273)-291)); % Surface Tension
LVNp=vlcp*(mjcp/(g*Gammap)^0.25);%Liquid Velocity Number CO2
NFRp=(vlcp^2)/(g*dwell); %Froude Number Co2 Injection 
    
    % Mentukan Pola Aliran
    lambdap=1;
    Eksp=log(lambdap);
    L1p=exp(-4.62-(3.757*Eksp)-(0.481*Eksp^2)-(0.0207*Eksp^3));
    L2p=exp(1.061-(4.602*Eksp)-(1.609*Eksp^2)-(0.179*Eksp^3)+(0.635*(10^-3)*Eksp^5));
    if ((NFRp<L1p)) || ((NFRp<L2p))
        polap = 1; %1=segregated
    elseif  ((NFRi>L1p) &&(NFRi>L2p))
            polap = 2;%2=distributed
    else
        polap=3; %3=intermittent
    end
    Hlscp=(0.98*(lambdap^0.4846)/(NFRp^0.0868));
    Hlicp=(0.84*(lambdap^0.5351)/(NFRp^0.0173));
    Hldcp=(1.06*(lambdap^0.5824)/(NFRp^0.0609));
     %Nilai C Aliran
    Csp=(1-lambdap)*log((4.7*LVNp^0.0868)/((lambdap^0.3692)*(NFRi^0.5056)));
    Cip=(1-lambdap)*log((4.7*LVNp^0.1244)/((lambdap^0.3692)*(NFRi^0.5056)));
    Cdp=(1-lambdap)*log((4.7*LVNp^0.1244)/((lambdap^0.3692)*(NFRi^0.5056)));
    if polap == 1
        HLP = Hlscp;
        Cp = Csp;
    elseif polap == 2
        HLP = Hldcp;
        Cp = Cdp;
    else HLP = Hlicp;
        Cp = Cip;
    end
    
    PSIp=1+Cp*(sin(1.8*90)-((1/3)*sin(90)^3));
    HLPtetha = HLP * PSIp; 
    yp = lambdap/(HLPtetha^2);
    Sp = (log(yp))/(-0.0523+3.812*log(yp)-0.8725*(log(yp)^2)+0.01853*(log(yp)^4));
    
    NREp = mjcp*vlcp*dwell/viscocp;
    fnsp = (2*log10(NREp/(((4.5223*log10(NREp))-3.8215))))^(-2);
    ftpp = exp(Sp)*fnsp;
   
    %Pressure Drop
    pdfp = (ftpp*minj*vlcp^2)/(2*gc*dwell);
    pdep = (g*mjcp)/gc;
    pdtotp = (pdfp+pdep)*0.0000442075*depth*3.28084;
    pco = pcp+pdtotp;
     %heat transfer co2 injection well
    Kcp=0.0174865883625187+-0.000124017589052688*tcp+0.0000167603402817045*pcp;   %Konduktivitas Thermal CO2 injection well (W/mK)
    Cpcp= (12.44045889320179+tcp*-0.0105800396299696+0.000240054091308146*pcp)*1000;%Heat capacity CO2 injection well (kJ/kgK)
    NPrcp=Cpcp*viscocp/Kcp;%Prandlt Number CO2 injection well
    NNucp=0.023*(NREp^0.8)*(NPrcp^0.3);%Nuselt Number CO2 injection well
    hcp=NNucp*Kcp/dwell;%heat transfer coefficient CO2 injection well (W/m2K)
    Rcp=thick/(Kcp*(3.14 *(dwell^2)/4));%Resistansi konduktivitas thermal CO2 injection well (K/W) 0.005-wallthicknes
    Ucp=hcp+(1/(Rcp*3.14*dwell^2));%All Heat transfer coefficient CO2 injection well (W/m2K)
    Qcp=minj*Cpcp*(tcp-tling);%energy panas CO2 injection well (watt)
    dtcp=Qcp/(Ucp*3.14*dwell*depth*2*minj) %Delta T CO2 injection well (C)

    %GAS RECOVERY
 mjng=35.34614832+-0.672856976292827*tcr+0.049907089*pcr; %Massa jenis natural gas (kg/m3)
 viscong=0.000010055+-1.33802257199713E-08*tcr+3.48620620536712E-09*pcr; %viskositas natural gas (kg/ms)
 Cpng=2.921948336+-0.006217212*tcr+0.00019548*pcr;
 kng=0.028944032+-0.00000363428*tcr+0.0000111577*pcr;
 Lpres=3.14*resthick*reslength;
 z = 0.86;
 ppz=2260.668973;
 pipzi=47867.359;
 Bg=2.8793*z*48.33/1960;
 M=(viscong/viscocr);
 fg=(1/(1+M));
 Sg=-0.0000005*(fg^5) + 0.00004*(fg^4) - 0.001*(fg^3) + 0.0071*(fg^2) + 0.0521*fg+ 0.2623;
 G = (Lling*reslength*poros*Sg)/Bg;
 Gp = G*(1-(ppz/pipzi))
 
 mco2inj=minj*3600*24;
 VCO2inj=mco2inj/mjci;
 tinj=G/VCO2inj;
 CNGR=0.6829*Gp;
 CNGRpd=CNGR/tinj;
 CCGR=(1-0.6829)*Gp;
 CCGRpd=CCGR/tinj;
 Rng=CNGRpd/28.263682*2.7;
 Rcg=CCGRpd*6.289814*0.935*68.096;
 Rtot=Rng+Rcg;
 
 
 co2price= 20.220; % 15.5432; %harga co2/ton
 %biaya pembelian co2 per day % DIGANTI
 bc=mco2inj*co2price/1000; %USD/day 
    
 listrik=0.06; %harga listrik industri di Ohio,USA per Mei 2020 % DIGANTI
 %biaya operasional pompa injeksi co2 per day
 kpci= (pco-pinj)*volci;%kinerja pompa co2 perjam 0.8efficiencypompa watt
 bopc= (kpci)*3600*listrik; %USD/day 0.06 harga listrik CA,USA per Mei 2020 % DIGANTI
    
 %biaya recycling co2
 brc= CCGRpd*20;%USD/day % DIGANTI

 profit=(Rtot-bc-brc-bopc)
end