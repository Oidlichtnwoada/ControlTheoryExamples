clear all
close all 
clc

syms kGSM iGSM dcGSM dvGSM uGSM RGSM iGSM dcP dvP Mext dqP wGSM wP phiGSMP cGSMP LGSM dGSMP JP JGSM iGSMR1 phiGSMPR1 wGSMR1 wPR1 iGSMR2 phiGSMPR2 wGSMR2 wPR2

assume(wP>0);
assume(wGSM>0);

%Zustandsdifferentialgleichungen
d1 = wGSM - wP;
d2 = (kGSM*(uGSM/RGSM-kGSM*wGSM/RGSM) - dcGSM - dvGSM*wGSM - (wGSM-wP)*dGSMP-phiGSMP*cGSMP)/JGSM;
d3 = ((wGSM-wP)*dGSMP+phiGSMP*cGSMP - dcP - dvP*wP - dqP*wP^2 - Mext)/JP;

eq1 = 0 == d1;
eq2 = 0 == d2;
eq3 = 0 == d3;

erg = solve(eq1,eq2,eq3, phiGSMP, wGSM, wP);
erg1 = [simplify(erg.phiGSMP(1)); simplify(erg.wGSM(1)); simplify(erg.wP(1))];
erg2 = [simplify(erg.phiGSMP(2)); simplify(erg.wGSM(2)); simplify(erg.wP(2))];

%Ruhelagen der ersten Lösung
phiGSMPR1 = erg1(1);
wGSMR1 = erg1(2);
wPR1 = erg1(3);

%Ruhelagen der zweiten Lösung
phiGSMPR2 = erg2(1);
wGSMR2 = erg2(2);
wPR2 = erg2(3);

%Linearisieren des Systems

%Ableiten
A1 = [diff(d1,phiGSMP),diff(d1,wGSM),diff(d1,wP)];
A2 = [diff(d2,phiGSMP),diff(d2,wGSM),diff(d2,wP)];
A3 = [diff(d3,phiGSMP),diff(d3,wGSM),diff(d3,wP)];
A = [A1;A2;A3];

%Ruhelagen einsetzen
ALsg1 = subs(subs(subs(A,phiGSMP,phiGSMPR1),wGSM,wGSMR1),wP,wPR1);
ALsg2 = subs(subs(subs(A,phiGSMP,phiGSMPR2),wGSM,wGSMR2),wP,wPR2);

%Ableiten
B1 = [diff(d1,uGSM), diff(d1, Mext)];
B2 = [diff(d2,uGSM), diff(d2, Mext)];
B3 = [diff(d3,uGSM), diff(d3, Mext)];
B = [B1;B2;B3];

%Ruhelagen einsetzen
BLsg1 = subs(subs(subs(B,phiGSMP,phiGSMPR1),wGSM,wGSMR1),wP,wPR1);
BLsg2 = subs(subs(subs(B,phiGSMP,phiGSMPR2),wGSM,wGSMR2),wP,wPR2);

%C und D wurden händisch abgeleitet
C = [0,0,1];
D = [0,0];

%Parameters
LGSM = 0.0014;
RGSM = 0.46;
kGSM = 0.1;
JGSM = 12.4E-3;
dcGSM = 0.152;
dvGSM = 1.8E-3;
JP = 32.5E-3;
dcP = 0.169;
dvP = 2.7E-3;
dqP = 1E-4;
cGSMP = 0.6822;
dGSMP = 1E-5;
uGSM = 5.6;
Mext = 0;

ALsg1 = double(subs(ALsg1));
ALsg2 = double(subs(ALsg2));
BLsg1 = double(subs(BLsg1));
BLsg2 = double(subs(BLsg2));

eig(ALsg1);
eig(ALsg2);

phiGSMP = double(subs(phiGSMPR2));
wGSM = double(subs(wGSMR2));
wP = double(subs(wPR2));

%zweite Lösung ist richtig, da nur bei dieser ein positiver Ausgang vorhanden ist
Ared = ALsg2;
Bred = BLsg2;
Cred = C;
Dred = D;

%System definieren
Ta = 10E-3;
linRedSysCont = ss(Ared,Bred, Cred, Dred);
linRedSysDisc = c2d(linRedSysCont,Ta,'zoh');

%Beobachtbarkeit 
rank(obsv(linRedSysDisc.A,linRedSysDisc.C));

Gq = tf(d2c(linRedSysDisc,'tustin'));
Gq = Gq(1);

%Matrizen für die weiteren Berechnungen
A = linRedSysDisc.A;
B = linRedSysDisc.B(:,1);
C = linRedSysDisc.C;
D = linRedSysDisc.D(:,1);

%Rang der Erreichbarkeitsmatrix == n => erreichbar
res = length(A) == rank(ctrb(A,B));

%komplexe Variable s erstellen
s = tf('s');

%Parameter
tr = 1;
wc = 1.2/tr;
ue = 0;
phi = 70-ue;

%Berechnungen
lambda0 = -5;
P = [lambda0,lambda0,lambda0];
k = -acker(A,B,exp(P*Ta));
g = 1/((C+D*k)*((eye(length(A))-A-B*k)\B)+D);

%geschlossener Regelkreis
AGR = A+B*k;
BGR = B*g;
CGR = C+D*k;
DGR = D*g;
GR = ss(AGR,BGR,CGR,DGR,Ta);

%Sprungantwort des geschlossenen Kreises
step(GR);

%Systemparameter um Integrator erweitern
AI = [A,zeros(3,1);-C,1];
BI = [B;0];
PI = [lambda0,lambda0,lambda0,lambda0];
kneu = -acker(AI,BI,exp(PI*Ta));

%Parameter für Regler
kI = kneu(4);
kP = 1/(C*inv((eye(length(A))-A))*B);
kX = kneu(1:3)+kP*C;

%Struct für Simulation
parZR.kI = kI;
parZR.kP = kP;
parZR.kX = kX;
parZR.C = C;

%Kompensationsregler
syms q
PoleGq = pole(Gq);
ZeroGq = zero(Gq);
res = simplify((q-PoleGq(1))*(q-PoleGq(2)));
res = double(coeffs(res));
F = res(1);
T = (res(3)/F)^0.5;
E = res(2)/F/2/T;

%gewählte Pole
qr1 = -10;
qr2 = -10;

%offener Kreis - VI und TI fehlen noch
L1 = Gq/s/(s-qr1)/(s-qr2)*(1+2*E*s*T+(s*T)^2);

%Auswerten der Phase und des Betrags an der Grenzfrequenz
[~,phase,~] = bode(L1, wc);

%Phase muss um a erhöht werden
a = mod(phi-180-phase,360);

%Berechnung von TI
TI = tand(a)/wc;

%offener Kreis ohne Verstärkung
L2 = L1*(1+s*TI);

%Auswerten der Phase und des Betrags an der Grenzfrequenz
[mag, phase, wout] = bode(L2, wc);

%Berechnung von VI
VI = mag^(-1);

%KRegler
KRegler = VI*(1+s*TI)*(1+2*E*s*T+(s*T)^2)/s/(s-qr1)/(s-qr1);
KompQ = (1+2*E*s*T+(s*T)^2)/(s-qr1)/(s-qr1);
InteQ = VI*(1+s*TI)/s;
KompZ = c2d(KompQ, Ta, 'tustin');
InteZ = c2d(InteQ, Ta, 'tustin');
[num,den] = tfdata(InteZ);
[r,p,o] = residue(cell2mat(num),cell2mat(den));
Vp = o;
Vi = r(1);
KReglerZ = c2d(KRegler, Ta, 'tustin');
KReglerZ = ss(KReglerZ);
[ARegler, BRegler, CRegler, DRegler] = ssdata(KReglerZ);

%vollständiger Beobachter
PBeobachter = [-20,-20,-20];
kBeobachter = -acker(linRedSysDisc.A',linRedSysDisc.C',exp(PBeobachter*Ta));
ABeo = linRedSysDisc.A + kBeobachter'*linRedSysDisc.C;
BBeo = [linRedSysDisc.B(:,1),-kBeobachter'];
CBeo = eye(3);
DBeo = zeros(3,2);

%vollständiger Beobachter Störung
PBeobachterS = [-20,-20,-20,-20];
AneuBeo = [linRedSysDisc.A,linRedSysDisc.B(:,2);zeros(1,4)];
kBeobachterS = -acker(AneuBeo',[linRedSysDisc.C,0]',exp(PBeobachterS*Ta));
ABeoS = AneuBeo + kBeobachterS'*[linRedSysDisc.C,0];
BBeoS = [[linRedSysDisc.B(:,1);0],-kBeobachterS'];
CBeoS = eye(4);
DBeoS = zeros(4,2);



