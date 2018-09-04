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
tic

%% Inputs

%{
Start and end times of observation plan.
%}

startTime = '02 Dec 2018 05:00:00';
endTime   = '02 Dec 2018 09:10:50';


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
minEnd = 10;
secEnd = 50;

mcStart1SigStDevRtn_1 = [1.63; 1.507; 1.015];
mcStart1SigStDevRtn_2 = [1.63; 1.507; 1.015];
[mcStart1SigStDevRtn] = VecElementInterp (mcStart1SigStDevRtn_1,mcStart1SigStDevRtn_2, minStart, secStart);



mcEnd1SigStDevRtn_1 = [1.41; 1.394; .905];
mcEnd1SigStDevRtn_2 = [1.356; 1.364; .878];
[mcEnd1SigStDevRtn] = VecElementInterp (mcEnd1SigStDevRtn_1,mcEnd1SigStDevRtn_2, minEnd, secEnd);



mcStartMeanPosJ2000_1 = [-20.99748514; 12.34661639; -0.603521407];
mcStartMeanPosJ2000_2 = [-20.99748514; 12.34661639; -0.603521407];
[mcStartMeanPosJ2000] = VecElementInterp (mcStartMeanPosJ2000_1,mcStartMeanPosJ2000_2, minStart, secStart);



mcEndMeanPosJ2000_1 = [-19.43182944	12.89589636	-0.352196672];
mcEndMeanPosJ2000_2 = [-19.03861787	13.03345095	-0.289119365];
[mcEndMeanPosJ2000] = VecElementInterp (mcEndMeanPosJ2000_1,mcEndMeanPosJ2000_2, minEnd, secEnd);



mcStartMeanVelJ2000_1 = [0.000108328; 3.81E-05; 1.74E-05];
mcStartMeanVelJ2000_2 = [0.000108328; 3.81E-05; 1.74E-05];
[mcStartMeanVelJ2000] = VecElementInterp (mcStartMeanVelJ2000_1,mcStartMeanVelJ2000_2, minStart, secStart);



mcEndMeanVelJ2000_1 = [0.000109125; 3.82E-05; 1.75E-05];
mcEndMeanVelJ2000_2 = [0.000109326; 3.82E-05; 1.75E-05];
[mcEndMeanVelJ2000] = VecElementInterp (mcEndMeanVelJ2000_1,mcEndMeanVelJ2000_2, minEnd, secEnd);



mcStart1SigVelJ2000_1 = [1.8842e-05; 1.1925e-05; 1.0326e-05]*1;
mcStart1SigVelJ2000_2 = [1.8842e-05; 1.1925e-05; 1.0326e-05]*1;
[mcStart1SigVelJ2000] = VecElementInterp (mcStart1SigVelJ2000_1,mcStart1SigVelJ2000_2, minStart, secStart);



mcEnd1SigVelJ2000_1 = [1.8847e-05; 1.1959e-05; 1.0353e-05]*1;
mcEnd1SigVelJ2000_2 = [1.8847e-05; 1.1959e-05; 1.0353e-05]*1;
[mcEnd1SigVelJ2000] = VecElementInterp (mcEnd1SigVelJ2000_1,mcEnd1SigVelJ2000_2, minEnd, secEnd);



CovStart1SigStDev_1 = [0.3962; 0.0793; 0.04916];
CovStart1SigStDev_2 = [0.3962; 0.0793; 0.04916];
[CovStart1SigStDev] = VecElementInterp (CovStart1SigStDev_1,CovStart1SigStDev_2, minStart, secStart);



CovEnd1SigStDev_1 = [0.3936; 0.09431; 0.05509];
CovEnd1SigStDev_2 = [0.3928; 0.09835; 0.05659];
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
raDec1 = [333.43425; 2.68283];
raDec2 = [333.73689; 1.15049];
raDec3 = [331.64329; 0.72051];

angx = (sqrt((raDec3(2)-raDec2(2))^2+(raDec3(1)-raDec2(1))^2))*pi/180;
angy = (sqrt((raDec2(2)-raDec1(2))^2+(raDec2(1)-raDec1(1))^2))*pi/180;
degSumStartx = angx+fov;
degSumStarty = angy+fov;

%{
end
%}
raDec1 = [330.32758; 2.29920];
raDec2 = [330.68316; 0.50969];
raDec3 = [328.11091; -0.02160];

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
nTD = 1000;

%{
Number of velocity perturbations.
%}
vTD = 100;

%{
%}
fileStartID = fopen('trainingDataStart.txt','w');
fprintf(fileStartID, '%14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s \n', 'PosP R', 'PosP T', 'PosP N', 'VelP R', 'VelP T', 'VelP N', 'X S Cov', 'Y S Cov');

fileEndID = fopen('trainingDataEnd.txt','w');
fprintf(fileEndID, '%14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s %14.7s \n', 'PosP R', 'PosP T', 'PosP N', 'VelP R', 'VelP T', 'VelP N', 'X E Cov', 'Y E Cov');

 [trainPertStartRtn, tpSigStartRtn] = ...
     TrainingDataRtnPosPerturbation (mcStart1SigStDevRtn,nTD);

%  trainPertStartRtn = [mcStart1SigStDevRtn(1)*2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  tpSigStartRtn = [2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  trainPertStartRtn = [zeros(1,nTD); mcStart1SigStDevRtn(2)*2*ones(1,nTD); zeros(1,nTD)];
%  tpSigStartRtn = [zeros(1,nTD); 2*ones(1,nTD); zeros(1,nTD)];
%  trainPertStartRtn = [mcStart1SigStDevRtn(1)*(-0.3023)*ones(1,nTD); mcStart1SigStDevRtn(2)*0.6614*ones(1,nTD); mcStart1SigStDevRtn(3)*1.7471*ones(1,nTD)];
%  tpSigStartRtn = [-0.3023*ones(1,nTD); 0.6614*ones(1,nTD); 1.7471*ones(1,nTD)];

[trainPertEndRtn, tpSigEndRtn] = ...
   TrainingDataRtnPosPerturbation (mcEnd1SigStDevRtn,nTD);

%  trainPertEndRtn = [mcEnd1SigStDevRtn(1)*2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  tpSigEndRtn = [2*ones(1,nTD); zeros(1,nTD); zeros(1,nTD)];
%  trainPertEndRtn = [zeros(1,nTD); mcEnd1SigStDevRtn(2)*2*ones(1,nTD); zeros(1,nTD)];
%  tpSigEndRtn = [zeros(1,nTD); 2*ones(1,nTD); zeros(1,nTD)];
%  trainPertEndRtn = [zeros(1,nTD); zeros(1,nTD); mcEnd1SigStDevRtn(3)*2*ones(1,nTD)];
%  tpSigEndRtn = [zeros(1,nTD); zeros(1,nTD); 2*ones(1,nTD)];
 
[tpPertVelStartJ2000,tpSigVelStartJ2000] = TrainingDataJ2000VelPerturbation ...
    (mcStart1SigVelJ2000,vTD);

[tpPertVelEndJ2000,tpSigVelEndJ2000] = TrainingDataJ2000VelPerturbation ...
    (mcEnd1SigVelJ2000,vTD);


rotMatStartMcJ2000_2_Rtn = RotMatJ2000_2_Rtn (mcStartMeanPosJ2000,mcStartMeanVelJ2000);

rotMatEndMcJ2000_2_Rtn = RotMatJ2000_2_Rtn (mcEndMeanPosJ2000,mcEndMeanVelJ2000);

mcStart1SigVelRtn = rotMatStartMcJ2000_2_Rtn*mcStart1SigVelJ2000;
mcEnd1SigVelRtn = rotMatEndMcJ2000_2_Rtn*mcEnd1SigVelJ2000;

tpStartStateJ2000 = zeros(6,nTD);
tpEndStateJ2000 = zeros(6,nTD);
rotMatStartPJ2000_2_Rtn = cell(nTD,vTD);
rotMatEndPJ2000_2_Rtn = cell(nTD,vTD);


for k = 1:nTD
    for l = 1:vTD
        
        if k == 1
            tpPertVelStartRtn(:,l) = rotMatStartMcJ2000_2_Rtn*tpPertVelStartJ2000(:,l);
            tpPertVelEndRtn(:,l) = rotMatEndMcJ2000_2_Rtn*tpPertVelEndJ2000(:,l);
            tpSigVelStartRtn(:,l) = tpPertVelStartRtn(:,l)./mcStart1SigVelRtn;
            tpSigVelEndRtn(:,l) = tpPertVelEndRtn(:,l)./mcEnd1SigVelRtn;
            
            tpStartVelJ2000(:,l) = ScStateStart(4:6,1) + tpPertVelStartJ2000(:,l);
            tpEndVelJ2000(:,l) = ScStateEnd(4:6,1) + tpPertVelEndJ2000(:,l);
        end
        
        trainPertStartJ2000 = rotMatStartMcJ2000_2_Rtn'*trainPertStartRtn(:,k);
        trainPertEndJ2000 = rotMatEndMcJ2000_2_Rtn'*trainPertEndRtn(:,k);

        tpStartPosJ2000 = -trainPertStartJ2000 + mcStartMeanPosJ2000;
        tpEndPosJ2000 = -trainPertEndJ2000 + mcEndMeanPosJ2000;

        tpStartStateJ2000(:,k) = [tpStartPosJ2000; ScStateStart(4:6,1)];
        tpEndStateJ2000(:,k) = [tpEndPosJ2000; ScStateEnd(4:6,1)];  

        rotMatStartPJ2000_2_Rtn{k,l} = RotMatJ2000_2_Rtn ...
            (tpStartStateJ2000(1:3,k),tpStartVelJ2000(:,l));

        rotMatEndPJ2000_2_Rtn{k,l} = RotMatJ2000_2_Rtn ...
            (tpEndStateJ2000(1:3,k),tpEndVelJ2000(:,l));
    
    end    
end

ct = nTD*vTD;
rotMatStartJ2000_2_Obs = zeros(3,3,nTD);
rotMatEndJ2000_2_Obs = zeros(3,3,nTD);
UEstartX = zeros(ct,1);
UEstartY = zeros(ct,1);
UEendX = zeros(ct,1);
UEendY = zeros(ct,1);
MosCovStartSigX = zeros(ct,1);
MosCovStartSigY = zeros(ct,1);
MosCovEndSigX = zeros(ct,1);
MosCovEndSigY = zeros(ct,1);    
ct = 0;

for k = 1:nTD
    for l = 1:vTD
        
        ct = ct+1;

        rotMatStartJ2000_2_Obs(:,:,k) = RotMatJ2000_2_Obs ...
            (-tpStartStateJ2000(1:3,k),SunStateStart(1:3,1));

        rotMatEndJ2000_2_Obs(:,:,k) = RotMatJ2000_2_Obs ...
            (-tpEndStateJ2000(1:3,k),SunStateEnd(1:3,1));

        obs3SigStart = rotMatStartJ2000_2_Obs(:,:,k)*rotMatStartPJ2000_2_Rtn{k,l}'*3*CovStart1SigStDev;

        obs3SigEnd = rotMatEndJ2000_2_Obs(:,:,k)*rotMatEndPJ2000_2_Rtn{k,l}'*3*CovEnd1SigStDev;

        UEstartX(ct,1) = (abs(obs3SigStart(1)) + rI);
        UEstartY(ct,1) = (abs(obs3SigStart(2)) + rI);

        UEendX(ct,1) = (abs(obs3SigEnd(1)) + rI);
        UEendY(ct,1) = (abs(obs3SigEnd(2)) + rI);

        rangeStartP = norm(tpStartStateJ2000(1:3,k)); 
        fovCovStartPx = atan(degSumStartx/2)*rangeStartP;
        fovCovStartPy = atan(degSumStarty/2)*rangeStartP;

        rangeEndP = norm(tpEndStateJ2000(1:3,k)); 
        fovCovEndPx = atan(degSumEndx/2)*rangeEndP;
        fovCovEndPy = atan(degSumEndy/2)*rangeEndP;   

        MosCovStartSigX(ct,1) = (fovCovStartPx-rI)/(UEstartX(ct,1)-rI)*3;
        MosCovStartSigY(ct,1) = (fovCovStartPy-rI)/(UEstartY(ct,1)-rI)*3;

        MosCovEndSigX(ct,1) = (fovCovEndPx-rI)/(UEendX(ct,1)-rI)*3;
        MosCovEndSigY(ct,1) = (fovCovEndPy-rI)/(UEendY(ct,1)-rI)*3;


        fprintf(fileStartID, '%14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f \n', ...
            tpSigStartRtn(1:3,k)', tpSigVelStartRtn(1:3,l)', MosCovStartSigX(ct,1), MosCovStartSigY(ct,1));
        
        fprintf(fileEndID, '%14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f \n', ...
            tpSigEndRtn(1:3,k)', tpSigVelEndRtn(1:3,l)', MosCovEndSigX(ct,1), MosCovEndSigY(ct,1));

    end
end

%% Recording Stats



ct = 0;
for k = 1:nTD
    for l = 1:vTD
        ct = ct+1;        
        matStartX(k,l) = MosCovStartSigX(ct,1);
        matStartY(k,l) = MosCovStartSigY(ct,1);
        matEndX(k,l) = MosCovEndSigX(ct,1);
        matEndY(k,l) = MosCovEndSigY(ct,1);
        
    end
end

meanStdVelConstStartX = mean(std(matStartX,0,1));
meanStdPosConstStartX = mean(std(matStartX,0,2));
meanStdVelConstStartY = mean(std(matStartY,0,1));
meanStdPosConstStartY = mean(std(matStartY,0,2));
meanStdVelConstEndX = mean(std(matEndX,0,1));
meanStdPosConstEndX = mean(std(matEndX,0,2));
meanStdVelConstEndY = mean(std(matEndY,0,1));
meanStdPosConstEndY = mean(std(matEndY,0,2));

ct = nTD*vTD;
confI = 1.6:0.2:3.0;

for k = 1:length(confI)
    conf = confI(k);
    aa1 = sum(MosCovStartSigX>=conf);
    aa2 = sum(MosCovStartSigY>=conf);
    bb1 = sum(MosCovEndSigX>=conf);
    bb2 = sum(MosCovEndSigY>=conf);

    probStartX(k) = (aa1/ct);
    probStartY(k) = (aa2/ct);
    probEndX(k) = (bb1/ct);
    probEndY(k) = (bb2/ct);
    
    standardNormal(k) = 1-normcdf(0,conf)*2;
    
    sNpSX(k) = probStartX(k)*standardNormal(k);
    sNpSY(k) = probStartY(k)*standardNormal(k);
    sNpEX(k) = probEndX(k)*standardNormal(k);
    sNpEY(k) = probEndY(k)*standardNormal(k);
    
end

sigBin = 0:0.01:10;
sigBinPlot = sigBin(1:end-1);

for k = 1:length(sigBin)-1
    
    countBinStartX(k,1) = sum(MosCovStartSigX >= sigBin(k)) - sum(MosCovStartSigX > sigBin(k+1));
    countBinStartY(k,1) = sum(MosCovStartSigY >= sigBin(k)) - sum(MosCovStartSigY > sigBin(k+1));
    countBinEndX(k,1) = sum(MosCovEndSigX >= sigBin(k)) - sum(MosCovEndSigX > sigBin(k+1));
    countBinEndY(k,1) = sum(MosCovEndSigY >= sigBin(k)) - sum(MosCovEndSigY > sigBin(k+1));
    
end

a = 2;
b = 2.5;
c = 3;

MosCovProbStartX1 = sum(MosCovStartSigX >= a)/(nTD*vTD);
MosCovProbStartY1 = sum(MosCovStartSigY >= a)/(nTD*vTD);
MosCovProbEndX1 = sum(MosCovEndSigX >= a)/(nTD*vTD);
MosCovProbEndY1 = sum(MosCovEndSigY >= a)/(nTD*vTD);
MosCovProbStartX2 = sum(MosCovStartSigX >= b)/(nTD*vTD);
MosCovProbStartY2 = sum(MosCovStartSigY >= b)/(nTD*vTD);
MosCovProbEndX2 = sum(MosCovEndSigX >= b)/(nTD*vTD);
MosCovProbEndY2 = sum(MosCovEndSigY >= b)/(nTD*vTD);
MosCovProbStartX3 = sum(MosCovStartSigX >= c)/(nTD*vTD);
MosCovProbStartY3 = sum(MosCovStartSigY >= c)/(nTD*vTD);
MosCovProbEndX3 = sum(MosCovEndSigX >= c)/(nTD*vTD);
MosCovProbEndY3 = sum(MosCovEndSigY >= c)/(nTD*vTD);




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
title('Probability Mosaic \geq a Specified UE Size')
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
title('Probability Mosaic \geq an UE Size that Contains Bennu')
axis([1.6 3 0 1])

figure(3)
plot(sigBinPlot,countBinStartX,'-b','linewidth',2); hold on
plot(sigBinPlot,countBinStartY,'-r','linewidth',2);
legend('X Mos Cov','Y Mos Cov','location','northEast');
set(gca,'fontsize', 14)
set(gcf,'color','w');
grid
xlabel('Sigma')
ylabel('# in Bins')
title('Start Histogram');
axis([0 10 0 1000])
xticks(0:1:10);
yo3 = ['P(X Mos Cov \geq 2.0) = ',num2str(MosCovProbStartX1,'%2.3f');...
    'P(Y Mos Cov \geq 2.0) = ',num2str(MosCovProbStartY1,'%2.3f');...
    'P(X Mos Cov \geq 2.5) = ',num2str(MosCovProbStartX2,'%2.3f');...
    'P(Y Mos Cov \geq 2.5) = ',num2str(MosCovProbStartY2,'%2.3f');...
    'P(X Mos Cov \geq 3.0) = ',num2str(MosCovProbStartX3,'%2.3f');...
    'P(Y Mos Cov \geq 3.0) = ',num2str(MosCovProbStartY3,'%2.3f')];
annotation(figure(3),'textbox',...
    [0.15 .8 0.2 0.1],...
    'String',yo3,'FitBoxToText','on');

figure(4)
plot(sigBinPlot,countBinEndX,'-b','linewidth',2); hold on
plot(sigBinPlot,countBinEndY,'-r','linewidth',2); 
legend('X Mos Cov','Y Mos Cov','location','northEast');
set(gca,'fontsize', 14)
set(gcf,'color','w');
grid
xlabel('Sigma')
ylabel('# in Bins')
title('End Histogram');
axis([0 10 0 1000])
xticks(0:1:10);
yo4 = ['P(X Mos Cov \geq 2.0) = ',num2str(MosCovProbEndX1,'%2.3f');...
    'P(Y Mos Cov \geq 2.0) = ',num2str(MosCovProbEndY1,'%2.3f');...
    'P(X Mos Cov \geq 2.5) = ',num2str(MosCovProbEndX2,'%2.3f');...
    'P(Y Mos Cov \geq 2.5) = ',num2str(MosCovProbEndY2,'%2.3f');...
    'P(X Mos Cov \geq 3.0) = ',num2str(MosCovProbEndX3,'%2.3f');...
    'P(Y Mos Cov \geq 3.0) = ',num2str(MosCovProbEndY3,'%2.3f')];
annotation(figure(4),'textbox',...
    [0.15 .8 0.2 0.1],...
    'String',yo4,'FitBoxToText','on');

toc
