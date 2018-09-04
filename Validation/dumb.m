clear all 
clc
close all

%%

r1 = 1008;
r2 = 1511;
rP = 1356;

rr = r2-r1;

% Loading frame and leap second kernel
cspice_furnsh ('orx_v09.tf');
cspice_furnsh ('naif0012.tls');

% Specifying UTC time and putting it in ephemeris time
utc1 = '12 Nov 2018 17:00:00';
utc2 = '03 Dec 2018 16:00:00';
et1 = cspice_str2et (utc1);
et2 = cspice_str2et (utc2);
timeRange = linspace(et1,et2,rr+1);
j2000Sec = datenum('01/01/2000 12:00:00')*24*60*60;

%% Perturbation_v2

spkP = 'orx_181112_181212_pert_wrt_MPRevB_v2.bsp';

cspice_furnsh(spkP);
id = cspice_spkobj(spkP,10);
timeWin = cspice_spkcov(spkP,id(1),100);
cspice_timout(timeWin(1),'YYYY-MM-DD THR:MN:SC.### ::RND')
cspice_timout(timeWin(2),'YYYY-MM-DD THR:MN:SC.### ::RND')

% Obtaining position in RTN frame of reference trajectory at time specified
stateP = cspice_spkpos( '-64', timeRange, 'ORX_RTN_BENNU', ...
                                        'NONE', 'Bennu' );
                                    
stateP_diff = cspice_spkpos( '-64', timeRange, 'J2000', ...
                                        'NONE', 'Bennu' );

data = xlsread('Approach_Deliv_Summary_RTN.xls');
dataJ2000 = xlsread('J2000_Summary.xls');

meanPosJ2000 = [ dataJ2000(r1:r2,6), dataJ2000(r1:r2,7), dataJ2000(r1:r2,8)];
meanVelJ2000 = [ dataJ2000(r1:r2,3), dataJ2000(r1:r2,4), dataJ2000(r1:r2,5)];

twoSigP_R = [data(r1:r2,4)*2, data(r1:r2,6)*0, data(r1:r2,8)*0];

funTime = zeros(length(timeRange),3);
PertPosJ2000 = zeros(length(timeRange),3);

for n = 1:length(timeRange)
    
    rotMatStartJ2000_2_Rtn = RotMatJ2000_2_Rtn (meanPosJ2000(n,:)',meanVelJ2000(n,:)');
    
    twoSigP_R_J2000 = rotMatStartJ2000_2_Rtn'*twoSigP_R(n,:)';
    
    PertPosJ2000(n,:) = meanPosJ2000(n,:)' + twoSigP_R_J2000;
    
    rotMatStartJ2000_2_Rtn_PP = RotMatJ2000_2_Rtn (PertPosJ2000(n,:)',meanVelJ2000(n,:)');
    
    funTime(n,:) = (stateP(:,n) - rotMatStartJ2000_2_Rtn_PP*meanPosJ2000(n,:)')';
        
end

figure(1)
plot((timeRange+j2000Sec)/(24*60*60),funTime(:,1),'b','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),twoSigP_R(:,1),'r','LineWidth',2)
grid
xlabel('Time');
ylabel('Perturbation (km)');
title('R Perturbation vs Time (2018)');
set(gca,'FontSize',12)
datetick('x', 'mmm dd')
legend('MC Traj Pert', 'Mean Pert')

figure(2)
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,1)'-stateP_diff(1,:),'b','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,2)'-stateP_diff(2,:),'r','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,3)'-stateP_diff(3,:),'g','LineWidth',4); hold on
grid
xlabel('Time');
ylabel('Perturbation (km)');
title('R Traj Diff vs Time (2018)');
set(gca,'FontSize',12)
datetick('x', 'mmm dd')
legend('x','y','z')


%% Perturbation_v3

spkP = 'orx_181112_181212_pert_wrt_MPRevB_v3.bsp';
cspice_furnsh(spkP);

% Obtaining position in RTN frame of reference trajectory at time specified
stateP = cspice_spkpos( '-64', timeRange, 'ORX_RTN_BENNU', ...
                                        'NONE', 'Bennu' );
stateP_diff = cspice_spkpos( '-64', timeRange, 'J2000', ...
                                        'NONE', 'Bennu' );

                                    
data = xlsread('Approach_Deliv_Summary_RTN.xls');
dataJ2000 = xlsread('J2000_Summary.xls');

meanPosJ2000 = [ dataJ2000(r1:r2,6), dataJ2000(r1:r2,7), dataJ2000(r1:r2,8)];
meanVelJ2000 = [ dataJ2000(r1:r2,3), dataJ2000(r1:r2,4), dataJ2000(r1:r2,5)];

twoSigP_T = [data(r1:r2,4)*0, data(r1:r2,6)*2, data(r1:r2,8)*0];

funTime = zeros(length(timeRange),3);
PertPosJ2000 = zeros(length(timeRange),3);

for n = 1:length(timeRange)
    
    rotMatStartJ2000_2_Rtn = RotMatJ2000_2_Rtn (meanPosJ2000(n,:)',meanVelJ2000(n,:)');
    
    twoSigP_T_J2000 = rotMatStartJ2000_2_Rtn'*twoSigP_T(n,:)';
    
    PertPosJ2000(n,:) = meanPosJ2000(n,:)' - twoSigP_T_J2000;
    
    rotMatStartJ2000_2_Rtn_PP = RotMatJ2000_2_Rtn (PertPosJ2000(n,:)',meanVelJ2000(n,:)');
    
    funTime(n,:) = (stateP(:,n) + rotMatStartJ2000_2_Rtn_PP*meanPosJ2000(n,:)')';
       
end




figure(3)
plot((timeRange(2:end)+j2000Sec)/(24*60*60),funTime(2:end,2),'b','LineWidth',4); hold on
plot((timeRange(2:end)+j2000Sec)/(24*60*60),twoSigP_T(2:end,2),'r','LineWidth',2)
datetick('x', 'mmm dd')
grid
xlabel('Time');
ylabel('Perturbation (km)');
title('T Perturbation vs Time (2018)');
set(gca,'FontSize',12)
datetick('x', 'mmm dd')
legend('MC Traj Pert', 'Mean Pert')

figure(4)
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,1)'-stateP_diff(1,:),'b','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,2)'-stateP_diff(2,:),'r','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,3)'-stateP_diff(3,:),'g','LineWidth',4); hold on
grid
xlabel('Time');
ylabel('Perturbation (km)');
title('T Traj Diff vs Time (2018)');
set(gca,'FontSize',12)
datetick('x', 'mmm dd')
legend('x','y','z')

%% Perturbation_v4

spkP = 'orx_181112_181212_pert_wrt_MPRevB_v4.bsp';
cspice_furnsh(spkP);

% Obtaining position in RTN frame of reference trajectory at time specified
stateP = cspice_spkpos( '-64', timeRange, 'ORX_RTN_BENNU', ...
                                        'NONE', 'Bennu' );
                                    
stateP_diff = cspice_spkpos( '-64', timeRange, 'J2000', ...
                                        'NONE', 'Bennu' );

data = xlsread('Approach_Deliv_Summary_RTN.xls');
dataJ2000 = xlsread('J2000_Summary.xls');

meanPosJ2000 = [ dataJ2000(r1:r2,6), dataJ2000(r1:r2,7), dataJ2000(r1:r2,8)];
meanVelJ2000 = [ dataJ2000(r1:r2,3), dataJ2000(r1:r2,4), dataJ2000(r1:r2,5)];

twoSigP_N = [data(r1:r2,4)*0, data(r1:r2,6)*0, data(r1:r2,8)*2];

funTime = zeros(length(timeRange),3);
PertPosJ2000 = zeros(length(timeRange),3);

for n = 1:length(timeRange)
    
    rotMatStartJ2000_2_Rtn = RotMatJ2000_2_Rtn (meanPosJ2000(n,:)',meanVelJ2000(n,:)');
    
    twoSigP_N_J2000 = rotMatStartJ2000_2_Rtn'*twoSigP_N(n,:)';
    
    PertPosJ2000(n,:) = meanPosJ2000(n,:)' - twoSigP_N_J2000;
    
    rotMatStartJ2000_2_Rtn_PP = RotMatJ2000_2_Rtn (PertPosJ2000(n,:)',meanVelJ2000(n,:)');
    
    funTime(n,:) = (stateP(:,n) + rotMatStartJ2000_2_Rtn_PP*meanPosJ2000(n,:)')';
    
end

figure(5)
plot((timeRange+j2000Sec)/(24*60*60),funTime(:,3),'b','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),twoSigP_N(:,3),'r','LineWidth',2)
grid
xlabel('Time');
ylabel('Perturbation (km)');
title('N Perturbation vs Time (2018)');
set(gca,'FontSize',12)
datetick('x', 'mmm dd')
legend('MC Traj Pert', 'Mean Pert')

figure(6)
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,1)'-stateP_diff(1,:),'b','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,2)'-stateP_diff(2,:),'r','LineWidth',4); hold on
plot((timeRange+j2000Sec)/(24*60*60),PertPosJ2000(:,3)'-stateP_diff(3,:),'g','LineWidth',4); hold on
grid
xlabel('Time');
ylabel('Perturbation (km)');
title('N Traj Diff vs Time (2018)');
set(gca,'FontSize',12)
datetick('x', 'mmm dd')
legend('x','y','z')

cspice_kclear;

