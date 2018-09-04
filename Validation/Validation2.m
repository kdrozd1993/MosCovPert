% Validation 2

clear all
clc
close all

mcStart1SigStDevRtn = [3.18; 1.537; 1.66];
mcStartMeanPosJ2000 = [-63.32346427; -3.248127288; -7.450282064];
mcStartMeanVelJ2000 = [8.81E-05; 3.42E-05; 1.46E-05];
CovStart1SigStDev = [0.4635; 0.1297; 0.1423];

%%
mcEnd1SigStDevRtn_1 = [3.214; 1.558; 1.677];
mcEnd1SigStDevRtn_2 = [3.223; 1.563; 1.681];

min = 10/60;
sec = 50/60/60;

y1 = [mcEnd1SigStDevRtn_1(1); mcEnd1SigStDevRtn_2(1)];
y2 = [mcEnd1SigStDevRtn_1(2); mcEnd1SigStDevRtn_2(2)];
y3 = [mcEnd1SigStDevRtn_1(3); mcEnd1SigStDevRtn_2(3)];
x = [0; 1];
xv = min+sec;

interpResult1 = interp1(x,y1,xv);
interpResult2 = interp1(x,y2,xv);
interpResult3 = interp1(x,y3,xv);

mcEnd1SigStDevRtn = [interpResult1; interpResult2; interpResult3];

mcEndMeanPosJ2000_1 = [-62.05087042; -2.754437244; -7.23975495];
mcEndMeanPosJ2000_2 = [-61.73127359; -2.630770615; -7.186970015];

min = 10/60;
sec = 50/60/60;

y1 = [mcEndMeanPosJ2000_1(1); mcEndMeanPosJ2000_2(1)];
y2 = [mcEndMeanPosJ2000_1(2); mcEndMeanPosJ2000_2(2)];
y3 = [mcEndMeanPosJ2000_1(3); mcEndMeanPosJ2000_2(3)];
x = [0; 1];
xv = min+sec;

interpResult1 = interp1(x,y1,xv);
interpResult2 = interp1(x,y2,xv);
interpResult3 = interp1(x,y3,xv);

mcEndMeanPosJ2000 = [interpResult1; interpResult2; interpResult3];

mcEndMeanVelJ2000_1 = [8.87E-05; 3.43E-05; 1.47E-05];
mcEndMeanVelJ2000_2 = [8.89E-05; 3.44E-05; 1.47E-05];

min = 10/60;
sec = 50/60/60;

y1 = [mcEndMeanVelJ2000_1(1); mcEndMeanVelJ2000_2(1)];
y2 = [mcEndMeanVelJ2000_1(2); mcEndMeanVelJ2000_2(2)];
y3 = [mcEndMeanVelJ2000_1(3); mcEndMeanVelJ2000_2(3)];
x = [0; 1];
xv = min+sec;

interpResult1 = interp1(x,y1,xv);
interpResult2 = interp1(x,y2,xv);
interpResult3 = interp1(x,y3,xv);

mcEndMeanVelJ2000 = [interpResult1; interpResult2; interpResult3];

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


%%

% mcStart1SigStDevRtn = mcEnd1SigStDevRtn;
% mcStartMeanPosJ2000 = mcEndMeanPosJ2000;
% mcStartMeanVelJ2000 = mcEndMeanVelJ2000;
% CovStart1SigStDev = CovEnd1SigStDev;


startTime = '27 Nov 2018 05:00:00';


cspice_furnsh ('orx_v09.tf');
cspice_furnsh ('naif0012.tls');

etStart = cspice_str2et (startTime);

spkP = 'orx_181112_181212_pert_wrt_MPRevB_v1.bsp';
cspice_furnsh(spkP);

rotMatJ2000_2_Rtn = RotMatJ2000_2_Rtn (mcStartMeanPosJ2000,mcStartMeanVelJ2000);

stateP_diff = cspice_spkezr( '-64', etStart, 'J2000', ...
                                        'NONE', 'Bennu' );
                                  
                                    
RTNmeanPos = rotMatJ2000_2_Rtn*mcStartMeanPosJ2000;
RTNPertPos = rotMatJ2000_2_Rtn*stateP_diff(1:3);

actualSigKM = rotMatJ2000_2_Rtn*(mcStartMeanPosJ2000-stateP_diff(1:3))

actualSigSIG = (RTNmeanPos-RTNPertPos)./mcStart1SigStDevRtn

rotMatJ2000_2_Rtn'*actualSigKM

state = -rotMatJ2000_2_Rtn'*actualSigKM + mcStartMeanPosJ2000


%%

cspice_furnsh ('orx_171006_190608_180710_od037-R-AM1-P-M3B_ort4_v1.bsp');
SunStateStart = cspice_spkezr( 'Sun', etStart, 'J2000', 'NONE', '-64');
cspice_furnsh ('orx_170117_210312_refMPRevB_v1.bsp')
velstate = cspice_spkezr( '-64', etStart, 'J2000', ...
                                        'NONE', 'Bennu' );

StartObs_Z = -stateP_diff(1:3)/norm(stateP_diff(1:3));
Cross_Y = cross(-state,SunStateStart(1:3,1));
StartObs_Y = Cross_Y/norm(Cross_Y);
StartObs_X = cross(StartObs_Y,StartObs_Z);

rotMatStartJ2000_2_Obs = [StartObs_X'; StartObs_Y'; StartObs_Z'];

rotMatJ2000_2_Rtn = RotMatJ2000_2_Rtn (state,velstate(4:6));
%rotMatJ2000_2_Rtn = RotMatJ2000_2_Rtn (state,stateP_diff(4:6));

obs3Sig = rotMatStartJ2000_2_Obs*rotMatJ2000_2_Rtn'*3*CovStart1SigStDev;

r = .2925;
x = (abs(obs3Sig(1)) + r)*1
y = (abs(obs3Sig(2)) + r)*1

%%
fov = 13.8/1000; %rad

%end
raDec1 = [4.4296; 6.62949];
raDec2 = [4.44791; 6.25672];
raDec3 = [4.04802; 6.23304];

%start
% raDec1 = [4.0943; 6.63043];
% raDec2 = [4.12394; 6.17611];
% raDec3 = [3.64947; 6.14189];

angx = (sqrt((raDec3(2)-raDec2(2))^2+(raDec3(1)-raDec2(1))^2))*pi/180;
angy = (sqrt((raDec2(2)-raDec1(2))^2+(raDec2(1)-raDec1(1))^2))*pi/180;
degSumx = angx+fov;
degSumy = angy+fov;


rangeP = norm(state); % m
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