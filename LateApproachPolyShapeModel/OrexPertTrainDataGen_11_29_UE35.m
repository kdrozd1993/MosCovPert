%{
Main Script Title:
OREX Perturbed Trajectory Training Data Generator

Author:
Kristofer Drozd

Description:
This code generates training data for the blank script. The training data
is created such that supervised machine learning can be used to classify 
whether a mosaic covers 1.6-sigma to 3-sigma uncertainty ellipse in 
0.2-sigma increments in the observer-frame based on how much the spacecraft
trajectory is perturbed in the radial, transvers, and normal components of 
the RTN frame. Decision boundaries will be created such that any 
perturbation can be classified in terms of mosaic coverage. 


Features are how much the spacecraft is perturbed at the beginning and end 
of the observation in the radial, transverse, and normal direction in terms 
of sigma. Output is mosaic coverage in sigma based on the 3-sigma 
uncertainty ellipse. The training data will be represented in a text file. 
Each row represents a different a perturbation and the columns represent 
the features and output.

Inputs:
 - Monte Carlo standard deviations at start and end of plan
 - spacecraft reference ephemeris SPK
 - frames kernel
 - Bennu and sun ephemeris SPK
 - Area of Mosaic at start and end of plan

Version:
 - 6/11/18 - Beta - development

     ."".    ."",
     |  |   /  /
     |  |  /  /
     |  | /  /
     |  |/  ;-._
     }  ` _/  / ;
     |  /` ) /  /
     | /  /_/\_/\
     |/  /      |
     (  ' \ '-  |
      \    `.  /
       |      |
       |      |

%}

clc
clear all;
close all;

%% Inputs

%{
Start and end times of observation plan.
%}

startTime = '29 Nov 2018 05:00:00';
endTime   = '29 Nov 2018 12:18:09';


%{ 
1 sigma standard deviation uncertainties with respect to the mean of 1000 
monte carlo trajectories. Uncertainties are based off of the RTN frame.
(km)

Mean position of monte carlo trajectories. Given in RTN and J2000 frames.
(km)

Mean velocity of monte carlo trajectories. Given in J2000 frame. (km/s)

1 sigma standard deviation uncertainties with respect to the reference 
trajectory associated with covariance analysis. Uncertainties are based off 
of the RTN frame. (km)
%}

minStart = 30;
secStart = 0;
minEnd = 18;
secEnd = 09;

mcStart1SigStDevRtn_1 = [3.572; 1.855; 1.864];
mcStart1SigStDevRtn_2 = [3.572; 1.855; 1.864];
[mcStart1SigStDevRtn] = VecElementInterp (mcStart1SigStDevRtn_1,mcStart1SigStDevRtn_2, minStart, secStart);



mcEnd1SigStDevRtn_1 = [3.621; 1.921; 1.894];
mcEnd1SigStDevRtn_2 = [3.628; 1.931; 1.899];
[mcEnd1SigStDevRtn] = VecElementInterp (mcEnd1SigStDevRtn_1,mcEnd1SigStDevRtn_2, minEnd, secEnd);



mcStartMeanPosJ2000_1 = [-47.43053042; 2.785733787; -4.854858844];
mcStartMeanPosJ2000_2 = [-47.43053042; 2.785733787; -4.854858844];
[mcStartMeanPosJ2000] = VecElementInterp (mcStartMeanPosJ2000_1,mcStartMeanPosJ2000_2, minStart, secStart);



mcEndMeanPosJ2000_1 = [-44.99706549; 3.687021197; -4.462766699];
mcEndMeanPosJ2000_2 = [-44.64696122; 3.816250017; -4.406447259];
[mcEndMeanPosJ2000] = VecElementInterp (mcEndMeanPosJ2000_1,mcEndMeanPosJ2000_2, minEnd, secEnd);



mcStartMeanVelJ2000_1 = [9.60E-05; 3.57E-05; 1.55E-05];
mcStartMeanVelJ2000_2 = [9.60E-05; 3.57E-05; 1.55E-05];
[mcStartMeanVelJ2000] = VecElementInterp (mcStartMeanVelJ2000_1,mcStartMeanVelJ2000_2, minStart, secStart);



mcEndMeanVelJ2000_1 = [9.72E-05; 3.59E-05; 1.56E-05];
mcEndMeanVelJ2000_2 = [9.73E-05; 3.59E-05; 1.57E-05];
[mcEndMeanVelJ2000] = VecElementInterp (mcEndMeanVelJ2000_1,mcEndMeanVelJ2000_2, minEnd, secEnd);



mcStart1SigVelJ2000_1 = [3.0095e-06; 1.6563e-06; 1.8651e-06]*1;
mcStart1SigVelJ2000_2 = [3.0095e-06; 1.6563e-06; 1.8651e-06]*1;
[mcStart1SigVelJ2000] = VecElementInterp (mcStart1SigVelJ2000_1,mcStart1SigVelJ2000_2, minStart, secStart);



mcEnd1SigVelJ2000_1 = [3.0274e-06; 1.6503e-06; 1.8612e-06]*1;
mcEnd1SigVelJ2000_2 = [3.0274e-06; 1.6503e-06; 1.8612e-06]*1;
[mcEnd1SigVelJ2000] = VecElementInterp (mcEnd1SigVelJ2000_1,mcEnd1SigVelJ2000_2, minEnd, secEnd);



CovStart1SigStDev_1 = [0.5545; 0.2794; 0.2730];
CovStart1SigStDev_2 = [0.5545; 0.2794; 0.2730];
[CovStart1SigStDev] = VecElementInterp (CovStart1SigStDev_1,CovStart1SigStDev_2, minStart, secStart);



CovEnd1SigStDev_1 = [0.5632; 0.2944; 0.2851];
CovEnd1SigStDev_2 = [0.5654; 0.2982; 0.2882];
[CovEnd1SigStDev] = VecElementInterp (CovEnd1SigStDev_1,CovEnd1SigStDev_2, minEnd, secEnd);


%{
Region of interest. (km)
%}
rI = .2925;

%{
Camera fov. (rad)
%}
fov = 13.8/1000; 

%{
start
%}
raDec1 = [359.41959; 6.87197];
raDec2 = [359.71800; 4.85111];
raDec3 = [357.70602; 4.54105];

angx = (sqrt((raDec3(2)-raDec2(2))^2+(raDec3(1)-raDec2(1))^2))*pi/180;
angy = (sqrt((raDec2(2)-raDec1(2))^2+(raDec2(1)-raDec1(1))^2))*pi/180;
degSumStartx = angx+fov;
degSumStarty = angy+fov;

%{
end
%}
raDec1 = [358.45668; 6.90791];
raDec2 = [358.81575; 4.58705];
raDec3 = [356.49438; 4.21310];

angx = (sqrt((raDec3(2)-raDec2(2))^2+(raDec3(1)-raDec2(1))^2))*pi/180;
angy = (sqrt((raDec2(2)-raDec1(2))^2+(raDec2(1)-raDec1(1))^2))*pi/180;
degSumEndx = angx+fov;
degSumEndy = angy+fov;



%{
Loading spice kernels.
%}
cspice_furnsh ('Kernels/orx_v09.tf');
cspice_furnsh ('Kernels/naif0012.tls');
cspice_furnsh ('Kernels/orx_171006_190608_180710_od037-R-AM1-P-M3B_ort4_v1.bsp');

etStart = cspice_str2et (startTime);
etEnd = cspice_str2et (endTime);

ScStateStart = cspice_spkezr( '-64', etStart, 'J2000', 'NONE', 'Bennu');                         
ScStateEnd = cspice_spkezr( '-64', etEnd, 'J2000', 'NONE', 'Bennu');                                    

SunStateStart = cspice_spkezr( 'Sun', etStart, 'J2000', 'NONE', '-64');                         
SunStateEnd = cspice_spkezr( 'Sun', etEnd, 'J2000', 'NONE', '-64');                                    


%{
Number of positional perturbations.
%}
nTD = 100;

%{
Number of velocity perturbations.
%}
vTD = 100;

%{
UE planned to.

CHANGE AS NEEDED
%}
pUE = 3.5;

%{
Opening the perturbation mosaic coverage files and placing the headers.
%}
fileStartID = fopen('trainingDataStart.txt','w');
fprintf(fileStartID, '%14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s \n', ...
    'PosP R', 'PosP T', 'PosP N', 'VelP R', 'VelP T', 'VelP N', 'X S Cov', 'Y S Cov');

fileEndID = fopen ('trainingDataEnd.txt','w');
fprintf(fileEndID, '%14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s \n', ...
    'PosP R', 'PosP T', 'PosP N', 'VelP R', 'VelP T', 'VelP N', 'X E Cov', 'Y E Cov');

%% The cold dark bloody magic

%{
Loading spice kernels:

frames kernel
leap second kernel
reference trajectory spk
%}
cspice_furnsh ('Kernels/orx_v09.tf');
cspice_furnsh ('Kernels/naif0012.tls');
cspice_furnsh ('Kernels/orx_171006_190608_180710_od037-R-AM1-P-M3B_ort4_v1.bsp');

%{
Getting start and end times in seconds past J2000 (ephemeris time).
%}
etStart = cspice_str2et (startTime);
etEnd = cspice_str2et (endTime);

%{
Reference state of the SC in J2000 from Bennu.
%}
ScStateStart = cspice_spkezr ( '-64', etStart, 'J2000', 'NONE', 'Bennu');                         
ScStateEnd = cspice_spkezr ( '-64', etEnd, 'J2000', 'NONE', 'Bennu');                                    

%{
State of the Sun in J2000 from the spacecraft on the reference trajectory.
%}
SunStateStart = cspice_spkezr ( 'Sun', etStart, 'J2000', 'NONE', '-64');                         
SunStateEnd = cspice_spkezr ( 'Sun', etEnd, 'J2000', 'NONE', '-64');                                    

%{
Aquireing the positional perturbations from the mean in the RTN frame. 
Commented code is used for validation purposes.
 %}
 [trainPertStartRtn, tpSigPosStartRtn] = ...
     TrainingDataRtnPosPerturbation (mcStart1SigStDevRtn, nTD);

%  trainPertStartRtn = [mcStart1SigStDevRtn(1)*2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  tpSigPosStartRtn = [2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  trainPertStartRtn = [zeros(1,nTD); mcStart1SigStDevRtn(2)*2*ones(1,nTD); zeros(1,nTD)];
%  tpSigPosStartRtn = [zeros(1,nTD); 2*ones(1,nTD); zeros(1,nTD)];
%  trainPertStartRtn = [mcStart1SigStDevRtn(1)*(-0.3023)*ones(1,nTD); mcStart1SigStDevRtn(2)*0.6614*ones(1,nTD); mcStart1SigStDevRtn(3)*1.7471*ones(1,nTD)];
%  tpSigPosStartRtn = [-0.3023*ones(1,nTD); 0.6614*ones(1,nTD); 1.7471*ones(1,nTD)];

[trainPertEndRtn, tpSigPosEndRtn] = ...
   TrainingDataRtnPosPerturbation (mcEnd1SigStDevRtn, nTD);

%  trainPertEndRtn = [mcEnd1SigStDevRtn(1)*2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  tpSigPosEndRtn = [2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  trainPertEndRtn = [zeros(1,nTD); mcEnd1SigStDevRtn(2)*2*ones(1,nTD); zeros(1,nTD)];
%  tpSigPosEndRtn = [zeros(1,nTD); 2*ones(1,nTD); zeros(1,nTD)];
%  trainPertEndRtn = [zeros(1,nTD); zeros(1,nTD); mcEnd1SigStDevRtn(3)*2*ones(1,nTD)];
%  tpSigPosEndRtn = [zeros(1,nTD); zeros(1,nTD); 2*ones(1,nTD)];
 
%{
Aquireing the velocity perturbations from the mean in the J2000 frame. 
 %}
[tpPertVelStartJ2000, tpSigVelStartJ2000] = TrainingDataJ2000VelPerturbation ...
    (mcStart1SigVelJ2000, vTD);

[tpPertVelEndJ2000, tpSigVelEndJ2000] = TrainingDataJ2000VelPerturbation ...
    (mcEnd1SigVelJ2000, vTD);

%{
Obtaining rotation matrix that transforms a vecotr from J2000 to the mean's
RTN frame.
 %}
rotMatStartMcJ2000_2_Rtn = RotMatJ2000_2_Rtn (mcStartMeanPosJ2000, mcStartMeanVelJ2000);

rotMatEndMcJ2000_2_Rtn = RotMatJ2000_2_Rtn (mcEndMeanPosJ2000, mcEndMeanVelJ2000);

%{
Aquireing the velocity standard deviations in the RTN frame. 
 %}
mcStart1SigVelRtn = rotMatStartMcJ2000_2_Rtn * mcStart1SigVelJ2000;
mcEnd1SigVelRtn = rotMatEndMcJ2000_2_Rtn * mcEnd1SigVelJ2000;


tpPertVelStartRtn = zeros (3, vTD);
tpPertVelEndRtn = zeros (3, vTD);
tpSigVelStartRtn = zeros (3, vTD);
tpSigVelEndRtn = zeros (3, vTD);
tpStartVelJ2000 = zeros (3, vTD);
tpEndVelJ2000 = zeros (3, vTD);
tpStartStateJ2000 = zeros (6, nTD);
tpEndStateJ2000 = zeros (6, nTD);
rotMatStartPJ2000_2_Rtn = cell (nTD, vTD);
rotMatEndPJ2000_2_Rtn = cell (nTD, vTD);


for k = 1:nTD
    
    for l = 1:vTD
        
        %{
        This if statement implys each positional perturbation will have the
        same velocity perturbations.
        %}
        if k == 1
            
            %{
            Obtaining the velocity perturbations in the mean's RTN frame.
            (km)
            %}
            tpPertVelStartRtn(:, l) = rotMatStartMcJ2000_2_Rtn * tpPertVelStartJ2000(:, l);
            tpPertVelEndRtn(:, l) = rotMatEndMcJ2000_2_Rtn * tpPertVelEndJ2000(:, l);
            
            %{
            Obtaining the velocity perturbations in the mean's RTN frame.
            (sigma)
            %}
            tpSigVelStartRtn(:, l) = tpPertVelStartRtn(:, l) ./ mcStart1SigVelRtn;
            tpSigVelEndRtn(:, l) = tpPertVelEndRtn(:, l) ./ mcEnd1SigVelRtn;
            
            %{
            Obtaining the velocity state of the spacecraft after 
            perturbation. (km/s)
            %}
            tpStartVelJ2000(:, l) = ScStateStart(4:6, 1) + tpPertVelStartJ2000(:, l);
            tpEndVelJ2000(:, l) = ScStateEnd(4:6, 1) + tpPertVelEndJ2000(:, l);
            
        end
        
        %{
        Obtaining the positional perturbation in the mean's RTN to J2000. 
        (km)
        %}
        trainPertStartJ2000 = rotMatStartMcJ2000_2_Rtn' * trainPertStartRtn(:, k);
        trainPertEndJ2000 = rotMatEndMcJ2000_2_Rtn' * trainPertEndRtn(:, k);
        
        %{
        Obtaining the positional state of each perturbed point. (km)
        %}
        tpStartPosJ2000 = trainPertStartJ2000 + mcStartMeanPosJ2000;
        tpEndPosJ2000 = trainPertEndJ2000 + mcEndMeanPosJ2000;

        %{
        Obtaining the complete state of each perturbed point. [(km) (km/s)]
        %}
        tpStartStateJ2000(:,k) = [tpStartPosJ2000; tpStartVelJ2000(:, l)];
        tpEndStateJ2000(:,k) = [tpEndPosJ2000; tpEndVelJ2000(:, l)];  

        %{
        Obtaining the rotation matrix of J2000 to the each perturbed points
        RTN frame.
        %}
        rotMatStartPJ2000_2_Rtn{k,l} = RotMatJ2000_2_Rtn ...
            (tpStartStateJ2000(1:3,k), tpStartVelJ2000(:,l));

        rotMatEndPJ2000_2_Rtn{k,l} = RotMatJ2000_2_Rtn ...
            (tpEndStateJ2000(1:3,k), tpEndVelJ2000(:,l));
    
    end 
    
end

ct = nTD * vTD;

rotMatStartJ2000_2_Obs = zeros (3, 3, nTD);
rotMatEndJ2000_2_Obs = zeros (3, 3, nTD);
UEstartX = zeros (ct, 1);
UEstartY = zeros (ct, 1);
UEendX = zeros (ct, 1);
UEendY = zeros (ct, 1);
MosCovStartSigX = zeros (ct, 1);
MosCovStartSigY = zeros (ct, 1);
MosCovEndSigX = zeros (ct, 1);
MosCovEndSigY = zeros (ct, 1); 
rangePosSigStart = zeros (ct, 1); 
rangePosSigEnd = zeros (ct, 1); 
rangeVelSigStart = zeros (ct, 1); 
rangeVelSigEnd = zeros (ct, 1); 
rangePosKmStart = zeros (ct, 1);
rangePosKmEnd = zeros (ct, 1);
rangeVelKmStart = zeros (ct, 1);
rangeVelKmEnd = zeros (ct, 1);
distFromBennuKm = zeros (ct, 1);
StartPosSigRtn_R = zeros (ct, 1); 
StartPosSigRtn_T = zeros (ct, 1); 
StartPosSigRtn_N = zeros (ct, 1); 
StartVelSigRtn_R = zeros (ct, 1); 
StartVelSigRtn_T = zeros (ct, 1); 
StartVelSigRtn_N = zeros (ct, 1); 
EndPosSigRtn_R = zeros (ct, 1); 
EndPosSigRtn_T = zeros (ct, 1); 
EndPosSigRtn_N = zeros (ct, 1); 
EndVelSigRtn_R = zeros (ct, 1); 
EndVelSigRtn_T = zeros (ct, 1); 
EndVelSigRtn_N = zeros (ct, 1); 
StartPosKmRtn_R = zeros (ct, 1); 
StartPosKmRtn_T = zeros (ct, 1); 
StartPosKmRtn_N = zeros (ct, 1); 
StartVelKmRtn_R = zeros (ct, 1); 
StartVelKmRtn_T = zeros (ct, 1); 
StartVelKmRtn_N = zeros (ct, 1); 
EndPosKmRtn_R = zeros (ct, 1); 
EndPosKmRtn_T = zeros (ct, 1); 
EndPosKmRtn_N = zeros (ct, 1); 
EndVelKmRtn_R = zeros (ct, 1); 
EndVelKmRtn_T = zeros (ct, 1); 
EndVelKmRtn_N = zeros (ct, 1); 
 
ct = 0;

for k = 1:nTD
    
    for l = 1:vTD
        
        ct = ct+1;
        
        %{
        Rotation matrix for J2000 to every perturbed point's Observation 
        frame.
        %}
        rotMatStartJ2000_2_Obs(:, :, k) = RotMatJ2000_2_Obs (...
            -tpStartStateJ2000(1:3, k), SunStateStart(1:3, 1));

        rotMatEndJ2000_2_Obs(:,:,k) = RotMatJ2000_2_Obs (...
            -tpEndStateJ2000(1:3, k), SunStateEnd(1:3, 1));

        %{
        Getting the 3 sigma RTN spherical ellipsoid to observation frame.
        %}
        obs3SigStart = rotMatStartJ2000_2_Obs(:, :, k) * ...
            rotMatStartPJ2000_2_Rtn{k, l}' * 3 * CovStart1SigStDev;

        obs3SigEnd = rotMatEndJ2000_2_Obs(:, :, k) * ...
            rotMatEndPJ2000_2_Rtn{k, l}' * 3 * CovEnd1SigStDev;

        %{
        Getting the size of the uncertainty ellipse on the observer plane.
        %}
        UEstartX(ct, 1) = (abs (obs3SigStart(1)) + rI);
        UEstartY(ct, 1) = (abs (obs3SigStart(2)) + rI);

        UEendX(ct, 1) = (abs (obs3SigEnd(1)) + rI);
        UEendY(ct, 1) = (abs (obs3SigEnd(2)) + rI);

        rangeStartP = norm (tpStartStateJ2000(1:3, k)); 
        fovCovStartPx = atan (degSumStartx / 2) * rangeStartP;
        fovCovStartPy = atan (degSumStarty / 2) * rangeStartP;

        rangeEndP = norm (tpEndStateJ2000(1:3, k)); 
        fovCovEndPx = atan (degSumEndx / 2) * rangeEndP;
        fovCovEndPy = atan (degSumEndy / 2) * rangeEndP;   

        %{
        Obtaining the mosaic coverage after the perturbation
        %}
        MosCovStartSigX(ct, 1) = (fovCovStartPx - rI) / ...
            (UEstartX(ct, 1) - rI) * 3;
        MosCovStartSigY(ct, 1) = (fovCovStartPy - rI) / ...
            (UEstartY(ct, 1) - rI) * 3;

        MosCovEndSigX(ct, 1) = (fovCovEndPx - rI) / ...
            (UEendX(ct, 1) - rI) * 3;
        MosCovEndSigY(ct, 1) = (fovCovEndPy - rI) / ...
            (UEendY(ct, 1) - rI) * 3;


        fprintf (fileStartID, '%14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f \n', ...
            tpSigPosStartRtn(1:3, k)', tpSigVelStartRtn(1:3, l)', MosCovStartSigX(ct, 1), MosCovStartSigY(ct, 1));
        
        fprintf (fileEndID, '%14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f \n', ...
            tpSigPosEndRtn(1:3, k)', tpSigVelEndRtn(1:3, l)', MosCovEndSigX(ct, 1), MosCovEndSigY(ct, 1));

        rangePosSigStart(ct, 1) = norm (tpSigPosStartRtn(1:3, k));
        rangePosSigEnd(ct, 1) = norm (tpSigPosEndRtn(1:3, k));
        rangeVelSigStart(ct, 1) = norm (tpSigVelStartRtn(1:3, l));
        rangeVelSigEnd(ct, 1) = norm (tpSigVelEndRtn(1:3, l));

        rangePosKmStart(ct, 1) = norm (trainPertStartRtn(1:3, k));
        rangePosKmEnd(ct, 1) = norm (trainPertEndRtn(1:3, k));
        rangeVelKmStart(ct, 1) = norm (tpPertVelStartRtn(1:3, l));
        rangeVelKmEnd(ct, 1) = norm (tpPertVelEndRtn(1:3, l));      
        
        distFromBennuKm(ct, 1) = norm (tpStartStateJ2000(1:3,k));
        
        StartPosSigRtn_R(ct,1) = tpSigPosStartRtn(1, k);
        StartPosSigRtn_T(ct,1) = tpSigPosStartRtn(2, k);
        StartPosSigRtn_N(ct,1) = tpSigPosStartRtn(3, k);
        StartVelSigRtn_R(ct,1) = tpSigVelStartRtn(1, l);
        StartVelSigRtn_T(ct,1) = tpSigVelStartRtn(2, l);
        StartVelSigRtn_N(ct,1) = tpSigVelStartRtn(3, l);
        
        EndPosSigRtn_R(ct,1) = tpSigPosEndRtn(1, k);
        EndPosSigRtn_T(ct,1) = tpSigPosEndRtn(2, k);
        EndPosSigRtn_N(ct,1) = tpSigPosEndRtn(3, k);
        EndVelSigRtn_R(ct,1) = tpSigVelEndRtn(1, l);
        EndVelSigRtn_T(ct,1) = tpSigVelEndRtn(2, l);
        EndVelSigRtn_N(ct,1) = tpSigVelEndRtn(3, l);
        
        StartPosKmRtn_R(ct,1) = trainPertStartRtn(1, k);
        StartPosKmRtn_T(ct,1) = trainPertStartRtn(2, k);
        StartPosKmRtn_N(ct,1) = trainPertStartRtn(3, k);
        StartVelKmRtn_R(ct,1) = tpPertVelStartRtn(1, l);
        StartVelKmRtn_T(ct,1) = tpPertVelStartRtn(2, l);
        StartVelKmRtn_N(ct,1) = tpPertVelStartRtn(3, l);
        
        EndPosKmRtn_R(ct,1) = trainPertStartRtn(1, k);
        EndPosKmRtn_T(ct,1) = trainPertStartRtn(2, k);
        EndPosKmRtn_N(ct,1) = trainPertStartRtn(3, k);
        EndVelKmRtn_R(ct,1) = tpPertVelEndRtn(1, l);
        EndVelKmRtn_T(ct,1) = tpPertVelEndRtn(2, l);
        EndVelKmRtn_N(ct,1) = tpPertVelEndRtn(3, l);
        
    end
    
end



%% Setting up classifier table

ct = nTD * vTD;

StartMosCovLabelX = zeros (ct, 1);
StartMosCovLabelY = zeros (ct, 1);
EndMosCovLabelX = zeros (ct, 1);
EndMosCovLabelY = zeros (ct, 1);
plannedUE = zeros (ct, 1);

for k = 1:ct
    
    StartMosCovLabelX(k, 1) = MosaicCoverageLabeling (MosCovStartSigX(k, 1));
    StartMosCovLabelY(k, 1) = MosaicCoverageLabeling (MosCovStartSigY(k, 1));
    EndMosCovLabelX(k, 1) = MosaicCoverageLabeling (MosCovEndSigX(k, 1));
    EndMosCovLabelY(k, 1) = MosaicCoverageLabeling (MosCovEndSigX(k, 1));
    plannedUE(k, 1) = pUE;
    
end

CovStart1SigStDev_R = ones (ct,1) * CovStart1SigStDev(1);
CovStart1SigStDev_T = ones (ct,1) * CovStart1SigStDev(2);
CovStart1SigStDev_N = ones (ct,1) * CovStart1SigStDev(3);

classifierTableSig = table (StartPosSigRtn_R, StartPosSigRtn_T, StartPosSigRtn_N, ...
    StartVelSigRtn_R, StartVelSigRtn_T, StartVelSigRtn_N, ...
    CovStart1SigStDev_R, CovStart1SigStDev_T, CovStart1SigStDev_N, ...
    rangePosSigStart, rangeVelSigStart, distFromBennuKm, plannedUE,...
    StartMosCovLabelX, StartMosCovLabelY);

classifierTableKm = table (StartPosKmRtn_R, StartPosKmRtn_T, StartPosKmRtn_N, ...
    StartVelKmRtn_R, StartVelKmRtn_T, StartVelKmRtn_N, ...
    rangePosKmStart, rangeVelKmStart, distFromBennuKm, plannedUE, ...
    StartMosCovLabelX, StartMosCovLabelY);

Biggy = classifierTableKm{:, 1:9};

[ranks,weights] = relieff(Biggy,StartMosCovLabelX,10);


%% Recording Stats
%{
Here stats are recorded on the data
%}


%{
Obtaining the mean standard deviations for when velocity is held constant
and when position is held constant.
%}
ct = 0;

matStartX = zeros (nTD,vTD);
matStartY = zeros (nTD,vTD);
matEndX = zeros (nTD,vTD);
matEndY = zeros (nTD,vTD);

for k = 1:nTD
    
    for l = 1:vTD
        
        ct = ct+1;        
        matStartX(k,l) = MosCovStartSigX(ct, 1);
        matStartY(k,l) = MosCovStartSigY(ct, 1);
        matEndX(k,l) = MosCovEndSigX(ct, 1);
        matEndY(k,l) = MosCovEndSigY(ct ,1);
        
    end
    
end

meanStdVelConstStartX = mean (std (matStartX, 0, 1));
meanStdPosConstStartX = mean (std (matStartX, 0, 2));
meanStdVelConstStartY = mean (std (matStartY, 0, 1));
meanStdPosConstStartY = mean (std (matStartY, 0, 2));
meanStdVelConstEndX = mean (std (matEndX, 0, 1));
meanStdPosConstEndX = mean (std (matEndX, 0, 2));
meanStdVelConstEndY = mean (std (matEndY, 0, 1));
meanStdPosConstEndY = mean (std (matEndY, 0, 2));

%{
Obtaining the Probability the mosaic coverage is greater than each sigma
value.
%}

ct = nTD*vTD;

confI = 1.6:0.2:3.0;

probStartX = zeros (length(confI), 1);
probStartY = zeros (length(confI), 1);
probEndX = zeros (length(confI), 1);
probEndY = zeros (length(confI), 1);
standardNormal = zeros (length(confI), 1);
sNpSX = zeros (length(confI), 1);
sNpSY = zeros (length(confI), 1);
sNpEX = zeros (length(confI), 1);
sNpEY = zeros (length(confI), 1);

for k = 1:length(confI)
    
    conf = confI(k);
    aa1 = sum (MosCovStartSigX >= conf);
    aa2 = sum (MosCovStartSigY >= conf);
    bb1 = sum (MosCovEndSigX >= conf);
    bb2 = sum (MosCovEndSigY >= conf);

    probStartX(k) = (aa1 / ct);
    probStartY(k) = (aa2 / ct);
    probEndX(k) = (bb1 / ct);
    probEndY(k) = (bb2 / ct);
    
    standardNormal(k) = 1 - normcdf (0, conf) * 2;
    
    sNpSX(k) = probStartX(k) * standardNormal(k);
    sNpSY(k) = probStartY(k) * standardNormal(k);
    sNpEX(k) = probEndX(k) * standardNormal(k);
    sNpEY(k) = probEndY(k) * standardNormal(k);
    
end

%{
Obtaining the number of times mosaic coverage is in a bin.
%}

sigBin = 0:0.01:10;

sigBinPlot = sigBin(1:end - 1);

countBinStartX = zeros(length(sigBin) - 1, 1);
countBinStartY = zeros(length(sigBin) - 1, 1);
countBinEndX = zeros(length(sigBin) - 1, 1);
countBinEndY = zeros(length(sigBin) - 1, 1);

for k = 1:length(sigBin) - 1
    
    countBinStartX(k,1) = sum (MosCovStartSigX >= sigBin(k)) - ...
        sum (MosCovStartSigX > sigBin(k + 1));
    
    countBinStartY(k,1) = sum (MosCovStartSigY >= sigBin(k)) - ...
        sum (MosCovStartSigY > sigBin(k + 1));
    
    countBinEndX(k,1) = sum (MosCovEndSigX >= sigBin(k)) - ...
        sum (MosCovEndSigX > sigBin(k + 1));
    
    countBinEndY(k,1) = sum (MosCovEndSigY >= sigBin(k)) - ...
        sum (MosCovEndSigY > sigBin(k + 1));
    
end

%{
Obtaining probability the mosaic coverage is greater than 2,2.5, & 3.
%}

a = 2;

b = 2.5;

c = 3;

MosCovProbStartXa = sum (MosCovStartSigX >= a) / (nTD * vTD);
MosCovProbStartYa = sum (MosCovStartSigY >= a) / (nTD * vTD);
MosCovProbEndXa = sum (MosCovEndSigX >= a) / (nTD * vTD);
MosCovProbEndYa = sum (MosCovEndSigY >= a) / (nTD * vTD);
MosCovProbStartXb = sum (MosCovStartSigX >= b) / (nTD * vTD);
MosCovProbStartYb = sum (MosCovStartSigY >= b) / (nTD * vTD);
MosCovProbEndXb = sum (MosCovEndSigX >= b) / (nTD * vTD);
MosCovProbEndYb = sum (MosCovEndSigY >= b) / (nTD * vTD);
MosCovProbStartXc = sum (MosCovStartSigX >= c) / (nTD * vTD);
MosCovProbStartYc = sum (MosCovStartSigY >= c) / (nTD * vTD);
MosCovProbEndXc = sum (MosCovEndSigX >= c) / (nTD * vTD);
MosCovProbEndYc = sum (MosCovEndSigY >= c) / (nTD * vTD);




figure(1)
plot(confI,probStartX,'-b','linewidth',2); hold on
plot(confI,probStartY,'-r','linewidth',2); hold on
plot(confI,probEndX,'--b','linewidth',2); hold on
plot(confI,probEndY,'--r','linewidth',2); hold on
plot(confI,standardNormal,'g','linewidth',2);
legend('Start X Mos Cov','Start Y Mos Cov','End X Mos Cov', 'End Y Mos Cov','Stand Norm','location','southwest');
set(gca,'fontsize', 14)
set(gcf,'color','w');
grid
ylabel('Probability');
xlabel('Uncertainty Ellipse Sigma the Mosaic Coverage is \geq Than');
title('Probability a Mosaic is \geq to a Specified UE Size')
axis([1.6 3 0 1])

figure(2)
plot(confI,sNpSX,'-b','linewidth',2); hold on
plot(confI,sNpSY,'-r','linewidth',2); hold on
plot(confI,sNpEX,'--b','linewidth',2); hold on
plot(confI,sNpEY,'--r','linewidth',2); 
legend('Start X Mos Cov','Start Y Mos Cov','End X Mos Cov', 'End Y Mos Cov','location','southwest');
set(gca,'fontsize', 14)
set(gcf,'color','w');
grid
ylabel('Probability');
xlabel('Uncertainty Ellipse Sigma the Mosaic Coverage is \geq Than');
title('Probability a Mosaic is \geq to an UE that Contains Bennu')
axis([1.6 3 0 1])

toc
