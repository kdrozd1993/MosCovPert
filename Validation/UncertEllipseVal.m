% Getting Uncertainty Ellipse Coverage on Observer Plane


clear all
close all
clc

startTime = '01 DEC 2018 05:00:00';
endTime   = '01 DEC 2018 12:18:09';


cspice_furnsh ('../Kernels/orx_v09.tf');
cspice_furnsh ('../Kernels/naif0012.tls');
cspice_furnsh ('../Kernels/orx_171006_190608_180710_od037-R-AM1-P-M3B_ort4_v1.bsp');

etStart = cspice_str2et (startTime);
etEnd = cspice_str2et (endTime);

SunStateStart = cspice_spkezr( 'Sun', etStart, 'J2000', 'NONE', '-64');                         
SunStateEnd = cspice_spkezr( 'Sun', etEnd, 'J2000', 'NONE', '-64');                                    

CovStart1SigStDev = [.4400; .1245; .1102];
CovEnd1SigStDev_1 = [0.4698; 0.1400; 0.1518];
CovEnd1SigStDev_2 = [0.4714; 0.1426; 0.1542];

min = 10/60;
sec = 50/60/60;

y1 = [CovEnd1SigStDev_1(1); CovEnd1SigStDev_2(1)];
y2 = [CovEnd1SigStDev_1(2); CovEnd1SigStDev_2(2)];
y3 = [CovEnd1SigStDev_1(3); CovEnd1SigStDev_2(3)];
x = [0; 1];
xv = min+sec;

interpResult1 = interp1(x,y1,xv);
interpResult2 = interp1(x,y2,xv);
interpResult3 = interp1(x,y3,xv);

CovEnd1SigStDev = [interpResult1; interpResult2; interpResult3];


% CovStart1SigStDev = CovEnd1SigStDev;

cspice_furnsh ('../Kernels/orx_181112_181212_pert_wrt_MPRevB_v2.bsp');

ScStateStart = cspice_spkezr( '-64', etStart, 'J2000', 'NONE', 'Bennu');          
rotMatJ2000_2_Rtn = cspice_pxform ('J2000','ORX_RTN_BENNU',etStart);


StartObs_Z = -ScStateStart(1:3,1)/norm(ScStateStart(1:3,1));
Cross_Y = cross(-ScStateStart(1:3,1),SunStateStart(1:3,1));
StartObs_Y = Cross_Y/norm(Cross_Y);
StartObs_X = cross(StartObs_Y,StartObs_Z);

rotMatStartJ2000_2_Obs = [StartObs_X'; StartObs_Y'; StartObs_Z'];

obs3Sig = rotMatStartJ2000_2_Obs*rotMatJ2000_2_Rtn'*3*CovStart1SigStDev;

r = .2925;
x = (abs(obs3Sig(1)) + r)*1
y = (abs(obs3Sig(2)) + r)*1

%%

fov = 13.8/1000; %rad

%start
raDec1 = [4.4296; 6.62949];
raDec2 = [4.44791; 6.25672];
raDec3 = [4.04802; 6.23304];

%end
% raDec1 = [4.0943; 6.63043];
% raDec2 = [4.12394; 6.17611];
% raDec3 = [3.64947; 6.14189];

angx = (sqrt((raDec3(2)-raDec2(2))^2+(raDec3(1)-raDec2(1))^2))*pi/180;
angy = (sqrt((raDec2(2)-raDec1(2))^2+(raDec2(1)-raDec1(1))^2))*pi/180;
degSumx = angx+fov;
degSumy = angy+fov;


cspice_furnsh ('../Kernels/orx_181112_181212_pert_wrt_MPRevB_v2.bsp');
ScStateStart = cspice_spkezr( '-64', etStart, 'J2000', 'NONE', 'Bennu');          

rangeP = norm(ScStateStart); % m
fovAreaPx = atan(degSumx/2)*rangeP;
fovAreaPy = atan(degSumy/2)*rangeP;

x1 = fovAreaPx;
y1 = fovAreaPy;
x0sig = 3;
y0sig = 3;
y0 = y;
x0 = x;

x1sig = (x1-r)/(x0-r)*x0sig
y1sig = (y1-r)/(y0-r)*y0sig

cspice_kclear;