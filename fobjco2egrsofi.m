function profit=fobjco2egr (x)
%optimized variable
minj=x(1) %laju aliran massa injeksi (kg/s)
pinj=x(2) %tekanan injeksi (psia)
tinj=x(3) %temepratur injeksi (C)

%Parameter model injection well
dwell=0.127; %diameter well (m)
g=9.8; %percepatan gravitasi (m/s2)
gc=1; %faktor gravitasi (kg m/N s2)
depth=2600; %kedalaman sumur (m)
thick=0.0053848; %ketebalan dinding (m)
rough=0.0000675; %kekasaran dinding (m)
trev=80.048; %temperatur reservoir (C)
tling=50; %temperatur lingkungan (C)

%INJECTION WELL

mjci= (861.751959949294+(0.00044477617075414*pinj)+(0.000392620606660566*tinj))/1.175; %massa jenis co2 injection well (kg/m3)
volci= minj/mjci; %Volume flowrate co2 (m3/s)
viscoci= 0.000096257645437626+(2.88475488287845E-09*tinj)+(9.72132894618796E-10*pinj); %viscositas co2 injection well (kg/ms)
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
    pdfi = (ftpi*minj*vlci^2)/(2*gc*dwell)%dp/dz friction
    pdei = (g*mjci)/gc;%dp/dz elevation
    pdtoti = (pdfi+pdei)*0.0000442075*depth*3.28084;%pressure drop(psi)
    
     %heat transfer co2 injection well
    Kci= 0.0258461471030894+(1.31403231349784E-08*tinj)+(3.71098287715882E-09*pinj);   %Konduktivitas Thermal CO2 injection well (W/mK)
    Cpci= 2330.98123930636-(0.241835581197486*tinj)-(0.0679611252333465*pinj);%Heat capacity CO2 injection well (kJ/kgK)
    NPrci=Cpci*viscoci/Kci;%Prandlt Number CO2 injection well
    NNuci=0.023*(NREi^0.8)*(NPrci^0.3);%Nuselt Number CO2 injection well
    hci=NNuci*Kci/dwell;%heat transfer coefficient CO2 injection well (W/m2K)
    Rci=thick/(Kci*(3.14 *(dwell^2)/4));%Resistansi konduktivitas thermal CO2 injection well (K/W) 0.005-wallthicknes
    Uci=hci+(1/(Rci*3.14*dwell^2))%All Heat transfer coefficient CO2 injection well (W/m2K)
    Qci=minj*Cpci*(tinj-trev)*200%energy panas CO2 injection well (watt)
    dtci=Qci/(Uci*3.14*dwell*depth*2) %Delta T CO2 injection well (C)
    
%RESERVOIR

tcr=tinj-dtci %temperature co2+oil pada reservoir C
pcr=(pinj+pdtoti)*1.1 %pressure Co2+oil pada reservoir psi
    
%parameter model
frac=0.5582; %fraksi liquid awal
resthick=120; %reservoir thickness (m)
reslength=200; %reservoir length (m)
Krock=0.000506126543186054; %thermal conductivity (W/mK)
permeabil= 2.26987E-16; %permeability (m2)
permeabilv= 2.26987E-17; %permeability vertikal (m2)
poros= 0.04; %porosity batuan
Lling = 11304; %Luas Lingkaran

    %pressure drop
    mjcr= -18815.7103042052+(238.074712376813*tcr)+(0.0230158247252583*pcr);%massa jenis phase liquid co2+Oil reservoir (kg/m3)
    volcr= minj/mjcr;%Volume flowrate co2+Oil reservoir (m3/s)
    vlcr=volcr/(3.14*((resthick^2)/4)); %velocity co2+oil reservoir (m/s)
    viscocr = -0.00257128806153902+(0.0000324094205486587*tcr)+(2.92008742548968E-09*pcr);%viscositas phase liquid co2+Oil reservoir (kg/ms)
    NRecr=mjcr*vlcr*resthick/viscocr; %reynold number co2+oil reservoir
    pdcr= -((viscocr*volcr*reslength)/(permeabil*3.14*resthick*reslength)*0.000145038);%pressure drop co2 reservoir psi. thicknes 35.052m

    %heat transfer
    Cpcr= -13770.7201705957+(-0.265156936748156*pcr)+(228.752765593037*tcr);  %Heat capacity CO2+oil di reservoir (J/kgK)
    Kcr= -4.74829179191263+(0.0600397249326838*tcr)+(5.80434279999517E-06*pcr); %Thermal conductivity co2+oil di reservoir (W/mK)
    NPrcr=Cpcr*viscocr/Kcr; %PRandlt Number CO2+oil di reservoir 1.73-thermal conductivity batuan
    NNucr=0.023*(NRecr^0.8)*(NPrcr^0.3);%Nuselt Number CO2+oil di reservoir
    Lpcr=(3.14*0.25*(resthick^2)*reslength)^(1/3); %length characteristic di reservoir
    hcr=NNucr*Krock/Lpcr ;%heat transfer coefficient di reservoir 5.678263 W/mK=Kthermal batuan
    Rkvcr=1/(hcr*(3.14*0.25*(resthick^2)));%tahanan konveksi perpindahan panas di rservoir
    Rkdcr=Lpcr/(Krock*3.14*0.25*(resthick^2));%tahanan konduksi perpindahan panas di reservoir
    Qcr=((tcr-80.048)/(Rkvcr+Rkdcr))*55*100;%energy panas co2+oil yg terbuang di reservoir (J)
    dtcr=Qcr/((Rkvcr+Rkdcr)*Cpcr); %Delta T Co2+oil di rservoir (C)

%PRODUCTION WELL

tcp=tcr-dtcr %temperature co2+oil production well (C)
pcp=pcr+pdcr%pressure co2+oil production well psi

mjcp= 547.749635477461-(0.000121094523548624*pcp)-(8.84095259528778E-06*tcp); %massa jenis co2 produksi well (kg/m3)
volcp= minj/mjcp; %Volume flowrate co2 (m3/s)
viscocp= 0.0000403651865689896+(2.76925923159356E-10*pcp)+(2.68440827059536E-11*tcp); %viscositas co2 injection well (kg/ms)
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
    
    NREp = mjcp*vlcp*dwell/viscocp;%massa jenis produksi*velocity*diameter wel/visco
    fnsp = (2*log10(NREp/(((4.5223*log10(NREp))-3.8215))))^(-2);
    ftpp = exp(Sp)*fnsp;
   
    %Pressure Drop
    pdfp = (ftpp*minj*vlcp^2)/(2*gc*dwell);%dp/dz friction
    pdep = (g*mjcp)/gc;%elevation
    pdtotp = (pdfp+pdep)*0.0000442075*depth*3.28084;%pressure drop
    pco = pcp+pdtotp;%Pco2 produksi well dipermukaan
     %heat transfer co2 injection well
    Kcp=-0.0207426758856552+(1.49501435213156E-06*pcp)+(0.000018929302557936*tcp);   %Konduktivitas Thermal CO2 injection well (W/mK)
    Cpcp= 2108.17762754639-(0.00301185424667767*tcp)-(0.0381364000806737*pcp);%Heat capacity CO2 injection well (kJ/kgK)
    NPrcp=Cpcp*viscocp/Kcp;%Prandlt Number CO2 injection well
    NNucp=0.023*(NREp^0.8)*(NPrcp^0.3);%Nuselt Number CO2 injection well
    hcp=NNucp*Kcp/dwell;%heat transfer coefficient CO2 injection well (W/m2K)
    Rcp=thick/(Kcp*(3.14 *(dwell^2)/4));%Resistansi konduktivitas thermal CO2 injection well (K/W) 0.005-wallthicknes
    Ucp=hcp+(1/(Rcp*3.14*dwell^2));%All Heat transfer coefficient CO2 injection well (W/m2K)
    Qcp=minj*Cpcp*(tcp-tling);%energy panas CO2 injection well (watt)
    dtcp=Qcp/(Ucp*3.14*dwell*depth*2*minj) %Delta T CO2 injection well (C)

    %GAS RECOVERY
 mjng=-14237.7271547194+(0.00207578076147587*pcr)+(182.455928567706*tcr); %Massa jenis gas (kg/m3)
 viscong=-0.00205546302813121+(2.95126688293442E-10*pcr)+(0.000961097651600211*tcr); %viskositas gas (kg/ms)
 Cpng= -8122.25347633883-(0.184614656999314*pcr)+(159.268435731113*tcr)/1000;%heat capacity mixture
 kng= 4.9490527096244+(0.0000797081912706261*pcr)-(0.0687648438977208*tcr);%kondutivitas mixture
 Lpres=3.14*resthick*reslength;%length vharacteristic
 z = 0.23;%derivation factor
 ppz=13571.429;
 pipzi=1408662.921;
 Bg=2.8793*z*80.048/3800;
 M=(viscong/viscocr);%mobiliti ratio
 fg=(1/(1+M));
 Sg=-0.0000005*(fg^5) + 0.00004*(fg^4) - 0.001*(fg^3) + 0.0071*(fg^2) + 0.0521*fg+ 0.2623;
 G = (Lling*reslength*poros*Sg)/Bg; %OGIP
 Gp = G*(1-(ppz/pipzi))
 
 mco2inj=minj*3600*24;%massa injeksi c02 perhari*24 jam*3600 detik
 VCO2inj=mco2inj/mjci;%massa co2 injeksi/massa jenis co2 pada injeksi
 tinj=G/(VCO2inj);%waktu injeksi = OGIP/volume injeksi
 CNGR=0.6829*Gp; %cummulative CH4 recovery(m3)=fraksi mol natgas*Gas recovery
 CNGRpd=CNGR/tinj;%cumulatif ch4 recovery/ waktu injeksi= recovery ch4(m3/day)
 CCGR=(1-0.6829)*Gp;%cumultive gas condensate(m3)
 CCGRpd=CCGR/tinj;%cumulative gas condensate/waktu inj injeksi= recovery condensate (m3/day)
 Rng=CNGRpd/28.263682*2.377; %cum gas condensat/(m3/day dikonversi MMBTU/day) *harga kondensate= pendapatan ch4
 Rcg=CCGRpd*6.289814*0.935*68.096;%pendapatan condensate=recovery condensate dkonversi (bbl/day)*harga kondensat per bbl/day
 Rtot=Rng+Rcg;
 
 
 co2price=15.5432; %harga co2/ton
 %biaya pembelian co2 per day
 bc=mco2inj*co2price/1000; %USD/day 
    
 listrik=0.1049; %harga listrik industri di Texas,USA
 %biaya operasional pompa injeksi co2 per day
 kpci= (pco-pinj)*volci;%kinerja pompa co2 perjam =(tekanan produksi well - p injeksi)* Volume flowrate co2 (m3/s)
 bopc= (kpci)*3600*listrik; %USD/day 0.1049 harga listrik texas,USA
    
 %biaya recycling co2
 brc= CCGRpd*15;%recovery gas condensate USD/day

 profit=(Rtot-bc-brc-bopc)
end