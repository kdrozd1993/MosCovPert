%-------------------------------------------------------------------------%
function [unitPert, sigmaPert] = ...
    TrainingDataRtnPosPerturbation (SigStDevRtn,nTD)
%-------------------------------------------------------------------------%

%{
Author:
Kristofer Drozd

Description:
This function generates Rtn perturbations from an initial vector that gives 
the 1 sigma standard devations in each direction. The perturbations are
obtained from a normal distrubtion with a mean of 0 and standard deviation
equal to the initial vector for each Rtn direction.

Inputs:
mcStart1SigStDevRtn - [3xnTD] matrix where each colunm is a normally
distributed random Rtn perturnation.

nTD - Amount of perturbations desired.

Version:

 - 6/11/18 - Beta - development

%}

trainPertRtnR = normrnd (0,SigStDevRtn(1),1,nTD);
trainPertRtnT = normrnd (0,SigStDevRtn(2),1,nTD);
trainPertRtnN = normrnd (0,SigStDevRtn(3),1,nTD);

unitPert = [trainPertRtnR; trainPertRtnT; trainPertRtnN];

sigmaPert = unitPert./SigStDevRtn;

end