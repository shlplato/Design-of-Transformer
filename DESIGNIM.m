clc
close all
clear all
n=input('Enter the synchronous speed in rpm')
ns=n/60 %synchronous speed in rps 
f=50 
p=2*f/ns % no.of poles 
rating=input('enter the kW rating of IM') 
%choosing values of magnetic and electrical loadings from the standards 
%according to the machine specifications required. The machine is designed 
%to have good performance alongwith lesser cost 
Bav=input('value of magnetic loading') %0.3 TO 0.6 Wb/m2
ac=input('value of electric loading') %5000 To 45000 amp.cond./m
Kw=0.955; %taking winding factor as 0.955
Co=11*Kw*Bav*ac*1e-3 %output cofficient 
%given efficiency=0.85 and pf(at full load)=0.83 
Q=rating/(0.85*0.83) % input kVA rating 
D2L=Q/(Co*ns) % as Q=Co*D2L*ns 
% For a cheap design ratio L/t(pole length to pole face length) should be 
%b/w 1.5 to 2. so for cheapest we 
% take L/t=1.5 ...so L/(pi*D/4), 
%where, D=inner Dia of stator 
% L=length of pole or core 
% t=pole face length 
% L/D=1.18; 
D=(D2L/1.18)^(1/3) 
L=1.18*D 
Li=L*0.9 %stacking factor=0.9 
t=L/1.5

%STATOR DESIGN----WINDING 

%Machine is designed to start with star-delta starter and operate as delta 
Es=400 %stator voltage per phase 
phim=Bav*t*L % Flux per pole 
Ts=floor(Es/(4.44*f*phim*Kw))%as Es=4.44*f*phim*Kw*Ts, stator turns/phase 
qs=3 ; %slot/pole/phase 
Ss=qs*p*3 %stator slots=qs*no. of poles*no. of phases 
% as there are many slots, thus slot harmonics and tooth pulsation is reduced 
yss=3.14*D*1e3/Ss % Stator slot pitch = pi*D/(stator slots) 
%stator slot pitch is allowable as we are using semi-closed slots 
Zss=floor(6*Ts/Ss) % Total Stator conductors/slot

% Now we are using single layer MUSH winding. As single layer thus the no. 
% of stator coils =1/2(stator slots). 
% We are not using double layer winding as then the slot area would be 
% large and as the slot pitch is very small thus mechanical strength of the 
% stator tooth will be poor

Cs=Ss/p % Coil span 
% as the coil span is odd (9), thus no need of shorting the coil and we 
% will use it as it is, 
angle_of_chording=0 
Kp=cos(angle_of_chording/2) %pitch factor 
%slot pitch=180/9=20deg 
Kd=sin(pi/6)/(3*sin(10*pi/180)) 
Kws=Kd*Kp % stator winding factor

% CONDUCTOR SIZE 
Is=Q*1e3/(3*Es) %stator current per phase 
Isl=sqrt(3)*Is 
% Taking current density as 4 A/mm2 
asc=Is/4 % area of stator conductor in mm2 
dsc=sqrt(4*asc/pi) 
% nearest standard diameter= 0.95mm 
ascn=pi*((0.95)^2)/4 % new area of stator conductor used 
cds=Is/ascn % current density in stator cond. 
%using medium covering the Dia of enamelled conductor is 1.041 mm

% SLOT DIMENSIONS 
%taking a space factor of 0.4 
As=round(Zss*ascn/0.4) % area of each slot 
% maximun allowable flux density is 1.7 Wb/mm2, so minimum tooth width to 
% keep flux density within limits
Wtsmin=phim*1e3/(1.7*Cs*Li) %min width in mm 
%tooth of width 6mm is chosen. lip=1mm, wedge=3mm 
swts=pi*(D*1e3+2*4)/24-6 % slot width of the portion near to rotor; 2*4 for 
%	wedge and lip for dimetrically opposite slots %slot width top of stator
%slot width at bottom= pi*(D+8+2*h)/24-6= 9.1+(pi*h/12) ; h=height of slot
%	area of conductor portion =.5*(slot width at bottom+slot width at top)*h
%	equating we get
h=11.1
swbs=9.1+(pi*h/12) %slot width at bottom of stator
dss=h+4 %depth of slot
Lmts=2*L+2.3*t+0.24 % length of mean turn

% STATOR TEETH
Bst=phim/(Cs*6*1e-3*Li) % flux density in stator teeth
%	flux density is coming out very small, thus we can increase the magnetic
%	loading ; taking Bst=0.8 Wb/m2

%STATOR CORE
flux_in_stator_core=phim/2
Acs=flux_in_stator_core/0.8 %Area of stator core 
dcs=round(Acs/Li*1e3) %depth of stator core 
Bsc=(Acs/Li)/dcs*0.8*1e3 % flux density in core 
Do=round(D*1e3+2*dss+2*dcs) % outsise dia of stator lamination

% ROTOR DESIGN
lg=0.2+2*sqrt(L*D) % air gap length
%	as flux density is very less so we can use large air gap so that overload
%	capacity of motor increases
Dr=D*1e3-2*lg

% ROTOR SLOTS
%no. of rotor slots=one pole pair less than stator poles 
Sr=Ss-2
ysr=pi*Dr/Sr % rotor slot pitch at air gap

%ROTOR BARS
Ib=2*3*Kws*Ts*Is*0.83/Sr %bar current
%taking bar current density=6A/mm2 
ab=Ib/6 % area of rotor bar
%standard size= 7mm*4mm with area=27.1mm2 
Wsr=4+0.3 %width of rotor slot 
dsr=7+1+1+0.15+0.15 % depth of rotor slot 
swbr=pi*(Dr-2*dsr)/Sr 
Brt=phim*p/(Sr*Li*(swbr-Wsr)*1e-3)
%extending bars 15mm and extra 10mm due to skewing 
Lb=L*1e3+2*15+10 % length of bar in mm 
rb=0.021*Lb*1e-3/27.1 %resistance of each bar 
ohmlossB=Sr*(Ib^2)*rb %total copper losses in bars

%END RINGS
Ier=Sr*Ib/(pi*4) %end ring current
%taking current density in end ring=7 A/mm2 
ae=Ier/7 %area
%using a ring of 10*7mm
de=10 % depth of ring
te=7 %thickness of ring

%so area of each end ring=70 mm2 
Doer=Dr-2*dsr %outer dia of ring 
Dier=Doer-2*de %inner dia of ring 
Dme=(Doer+Dier)/2 %mean dia 
re=0.021*pi*Dme*1e-3/70 
ohmlossER=2*(Ier^2)*re %copper loss in end ring 
ohmlossTR=ohmlossB+ohmlossER %total copper loss in rotor

%	(rotor copper loss/rotor output)=s/(1-s) 
s=ohmlossTR/(rating*1e3+ohmlossTR) % full load slip

%	ROTOR CORE
%depth of rotor core is same as stator core
dcr=dcs
Di=Dr-2*dsr-2*dcr % inner dia of rotor lamination

%	NO LOAD CURRENT
%	MAGNETIZING CURRENT

% i) AIR GAP
Wos=2 % stator slot opening
ratios=Wos/0.4 %(slot opening/gap length)
%	for this value of ratio the carter's coffecient for semi closed slots
%	is
Kcss=0.6
Kgss=yss/(yss-Kcss*Wos) % gap contraction factor for stator slots

Wor=1.5 % rotor slot opening
ratior=Wor/0.4
%	for this value of ratio the carter's coffecient for semi closed slots
%	is
Kcsr=0.5
Kgsr=ysr/(ysr-Kcsr*Wor) % gap contraction factor for rotor slots

Kgs=Kgss*Kgsr % gap contraction factor for slots

%as there are no ducts so Kgd=1 for them
Kg=Kgs*1 % overall gap contraction factor

Ag=pi*D*L/4 %area of air gap per pole
Bg60=1.36*Bav
lge=lg*Kg % effective air gap length 
ATg=800000*Bg60*lg*1e-3*Kg %Ampere turns required for air gap

% ii) STATOR TEETH
Ast=(Ss/p)*6*1e-3*Li % tooth width=6mm; are aof stator teeth /pole 
Bst60=1.36*Bst % flux density of stator teeth
%	corresponding to this value of B ampere turns/meter required are 
ATst=300
mmfst=ATst*dss*1e-3

%	iii) STATOR CORE
Asc=Li*dcs*1e-3 % area of stator core
lcs=pi*(D+2*dss*1e-3+dcs*1e-3)/(3*pi) % length of magnetic path through stator
%	corresponding to this value of B ampere turns/meter required are 
ATsc=200
mmfsc=ATsc*lcs

%	iv) ROTOR TEETH
Wtr13=pi*(Dr-4*dsr/3)/Sr-Wsr % width of rotor teeth at 1/3 height from
% narrow end
Atr=Sr/4*Wtr13*1e-3*Li
Btr60=Brt*1.36*0.85 % width of rotor teeth at 1/3 height(0.85 for 1/3 hgt)
%	corresponding to this value of B ampere turns/meter required are 
ATrt=700
mmfrt=ATrt*dsr*1e-3

% v) ROTOR CORE
Acr=dcr*Li % rotor core area
%	corresponding to this value of B(Bsc) ampere turns/meter required are 
ATrc=200
lcr=pi*Di*1e-3/(3*p) 
mmfrc=ATrc*lcr

AT60=ATg+mmfst+mmfsc+mmfrt+mmfrc % total mmf required

Im=0.427*p*AT60/(Kws*Ts)

%	LOSSES
%	Iron Loss in stator teeth

Vst=Ast*p*dss*1e-3 % volume of stator teeth
WTst=Vst*7.6*1e3 % weight of teeth
Bmst=pi*Bst/2 % maximum B at stator teeth
%	for this value of Bmst, specific iron loss=5.5W/kg 
ironlossST=5.5*WTst % stator teeth

%	Iron Loss in stator core
Vsc=(D+2*dss*1e-3)*pi*L*dcs*1e-3
WTsc=Vsc*7.6*1e3

%	corresponding to the flux density=Bsc, specific iron loss=2.8 W/kg 
ironlossSC=2.8*WTsc % iron loss in stator core 
ironlossTS=2*(ironlossST+ironlossSC) % normally total loss is taken twice
                                     %	the calculated value
%	FRICTION & WINDAGE LOSS
%	with use of ball bearing FW losses are about 1.5% of output

FWloss=1.5*rating*1e3/100
NLloss=ironlossTS+FWloss % total no load losses
Inll=NLloss/(3*Es) % loss component of no load current per phase
Io=sqrt(Im^2+Inll^2)
no_load_current_as_percentage_of_full_load_currrent=Io/Is*100
no_load_pf=Inll/Io
phio=(acos(no_load_pf))*180/pi

%SHORT CIRCUIT CURRENT
%Leakage Reactance
%stator slot leakage

Pss=4*3.14*(10^-7)*((2*h/(3*(swts+swbs)))+(2*3/(swbs+2))+(1/2)) %specific slot permeance for a tapered slot 
                                        %(considering h2 also to be occupied by conductors)
%rotor slot leakage
Psr=(4*3.14*10^-7)*((7/(3*6.8))+(2*1/(6.8+1.5)+(1/1.5)))
%4*3.14*10^-7*((h1/3*Ws)+(h2/Ws)+(2*h3/(Ws+W0))+(h4/W0))
%specific slot permeance for a parallel sided slot
Psr1=Psr*Kws^2*Ss/(1^2*Sr) %referred to stator side 
Ps=Pss+Psr1 %total specific slot permeance 
xs=8*3.14*f*Ts^2*L*(Ps/(p*3)) %q=3

%overhang leakage
%coil span / pole pitch =1 so corresponding Ks=1 
LoPo=4*3.14*10^-7*1*t^2/(3.14*yss) 
x0=8*3.14*f*Ts^2*(LoPo/(p*3)) %overhang leakage reactance

%Zigzag leakage
Xm=Es/Im %magnetizing reactance
qr=Sr/(p*3)

xz=(5*Xm/(6*3^2))*(1/3^2+1/qr^2) %zigzag leakage reactance per phase 
%the differential leakage reactance can be ignored in case of squirrel

%cage induction motors
Xs=xs+x0+xz %total leakage reactance per phase referred to stator

%Resistance
rs=0.021*Ts*Lmts/ascn
ohmlossS=3*Is^2*rs
ohmlossTRpp=ohmlossTR/3
rr1=(ohmlossTRpp)/((Is*0.83)^2)
Rs=rs+rr1 %total resistance referred to stator

%Impedeance
Zs=(Xs^2+Rs^2)^0.5 %total impedance of the rotor at standstill 
Isc=Es/Zs %short circuited current per phase 
scpf=Rs/Zs %short circuit power factor
phase_angle_of_short_circuit_current=(acos(scpf))*180/pi

%Losses and efficiency
total_loss_at_full_load=ohmlossS+ohmlossTR+ironlossTS+FWloss
%total_loss_at_full_load=total_stator_copper_loss+total_rotor_copper_loss+
% total_iron_loss+friction_and_windage_loss 
input_at_full_load=rating*1e3+total_loss_at_full_load
%output_at_full_load+total_loss_at_full_load 
efficiency_at_full_load=rating*1e3*100/input_at_full_load

