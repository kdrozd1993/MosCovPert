% Velocity Error 

clear all
clc
close all
addpath('..')

CovStart1SigStDev = [0.4635; 0.1297; 0.1423];
rI = .2925;
nTD = 1000;
fov = 13.8/1000;

raDec1 = [4.4296; 6.62949];
raDec2 = [4.44791; 6.25672];
raDec3 = [4.04802; 6.23304];
angx = (sqrt((raDec3(2)-raDec2(2))^2+(raDec3(1)-raDec2(1))^2))*pi/180;
angy = (sqrt((raDec2(2)-raDec1(2))^2+(raDec2(1)-raDec1(1))^2))*pi/180;
degSumStartx = angx+fov;
degSumStarty = angy+fov;

cspice_furnsh ('../Kernels/orx_v09.tf');
cspice_furnsh ('../Kernels/naif0012.tls');

Time = '01 Dec 2018 09:03:40';
et = cspice_str2et (Time);

cspice_furnsh ('../Kernels/orx_181112_181212_pert_wrt_MPRevB_v1.bsp');
ScState1 = cspice_spkezr( '-64', et, 'J2000', 'NONE', 'Bennu');

cspice_furnsh ('../Kernels/orx_181112_181212_pert_wrt_MPRevB_v2.bsp');
ScState2 = cspice_spkezr( '-64', et, 'J2000', 'NONE', 'Bennu');

cspice_furnsh ('../Kernels/orx_181112_181212_pert_wrt_MPRevB_v3.bsp');
ScState3 = cspice_spkezr( '-64', et, 'J2000', 'NONE', 'Bennu');

cspice_furnsh ('../Kernels/orx_181112_181212_pert_wrt_MPRevB_v4.bsp');
ScState4 = cspice_spkezr( '-64', et, 'J2000', 'NONE', 'Bennu');

velX = [ScState1(4); ScState2(4); ScState3(4); ScState4(4)];
velY = [ScState1(5); ScState2(5); ScState3(5); ScState4(5)];
velZ = [ScState1(6); ScState2(6); ScState3(6); ScState4(6)];

stdX = std(velX)
stdY = std(velY)
stdZ = std(velZ)

%%
cspice_furnsh ('../Kernels/orx_171006_190608_180710_od037-R-AM1-P-M3B_ort4_v1.bsp');
cspice_furnsh ('../Kernels/sb-101955-76.bsp');
ScState = cspice_spkezr( '-64', et, 'J2000', 'NONE', 'Bennu');

PertVelJ2000X = normrnd (ScState(4),stdX*1,1,nTD);
PertVelJ2000Y = normrnd (ScState(5),stdY*1,1,nTD);
PertVelJ2000Z = normrnd (ScState(6),stdZ*1,1,nTD);

cspice_furnsh ('../Kernels/orx_171006_190608_180710_od037-R-AM1-P-M3B_ort4_v1.bsp');
SunState = cspice_spkezr( 'Sun', et, 'J2000', 'NONE', '-64');                         


for k = 1:nTD
    
    PertStateVelX = ScState(4)+PertVelJ2000X(1,k);
    PertStateVelY = ScState(5)+PertVelJ2000Y(1,k);
    PertStateVelZ = ScState(6)+PertVelJ2000Z(1,k);
    PertState(:,k) = [ScState(1);ScState(2);ScState(3);PertStateVelX;PertStateVelY;PertStateVelZ];
    
    velocityVec = PertState(4:6,k);
    
    angle(k,1) = acosd(dot(velocityVec,ScState(4:6))/(norm(velocityVec)*norm(ScState(4:6))));
    
    rotMatStartPJ2000_2_Rtn(:,:,k) = RotMatJ2000_2_Rtn ...
    (PertState(1:3,k),PertState(4:6,k));
    
    rotMatStartJ2000_2_Obs(:,:,k) = RotMatJ2000_2_Obs ...
        (-PertState(1:3,k),SunState(1:3,1));
    
    obs3Sig = rotMatStartJ2000_2_Obs(:,:,k)*rotMatStartPJ2000_2_Rtn(:,:,k)'*3*CovStart1SigStDev;
    
    UEX(k,1) = (abs(obs3Sig(1)) + rI);
    UEY(k,1) = (abs(obs3Sig(2)) + rI);
    
    rangeStartP = norm(PertState(1:3,k));
    fovCovPx = atan(degSumStartx/2)*rangeStartP;
    fovCovPy = atan(degSumStarty/2)*rangeStartP;
    
    MosCovSigX(k,1) = (fovCovPx-rI)/(UEX(k,1)-rI)*3;
    MosCovSigY(k,1) = (fovCovPy-rI)/(UEY(k,1)-rI)*3;
    
    
    
end

aa = abs(3-MosCovSigX);

aa1 = sum(aa<.1);
aa2 = sum(aa<.2)-sum(aa<=.1);
aa3 = sum(aa>=.2);

prob1 = (aa1/nTD)*100
prob2 = (aa2/nTD)*100
prob3 = (aa3/nTD)*100
                       
