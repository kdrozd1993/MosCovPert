fun = vertcat(classifierTableSig_11_27_UE30{:,:}, ...
    classifierTableSig_11_27_UE35{:,:}, ...
    classifierTableSig_11_27_UE40{:,:}, ...
    classifierTableSig_11_29_UE30{:,:}, ...
    classifierTableSig_11_29_UE35{:,:});

StartPosSigRtn_R = fun(:,1);
StartPosSigRtn_T = fun(:,2);
StartPosSigRtn_N = fun(:,3);
StartVelSigRtn_R = fun(:,4);
StartVelSigRtn_T = fun(:,5);
StartVelSigRtn_N = fun(:,6);
CovStart1SigStDev_R = fun(:,7);
CovStart1SigStDev_T = fun(:,8);
CovStart1SigStDev_N = fun(:,9);
rangePosSigStart = fun(:,10);
rangeVelSigStart = fun(:,11);
distFromBennuKm = fun(:,12);
plannedUE = fun(:,13);
StartMosCovLabelX = fun(:,14);
StartMosCovLabelY = fun(:,15);


classifierTable = table (StartPosSigRtn_R, StartPosSigRtn_T, StartPosSigRtn_N, ...
    StartVelSigRtn_R, StartVelSigRtn_T, StartVelSigRtn_N, ...
    CovStart1SigStDev_R, CovStart1SigStDev_T, CovStart1SigStDev_N, ...
    rangePosSigStart, rangeVelSigStart, distFromBennuKm, plannedUE, ...
    StartMosCovLabelX, StartMosCovLabelY);