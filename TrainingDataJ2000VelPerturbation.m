%-------------------------------------------------------------------------%
function [unitPert, sigmaPert] = ...
    TrainingDataJ2000VelPerturbation (SigStDevJ2000,nTD)
%-------------------------------------------------------------------------%

%{
Author:
Kristofer Drozd

Description:
This function generates J2000 Velocity perturbations from an initial vector
that gives the 1 sigma standard devations in each direction. The 
perturbations are obtained from a normal distrubtion with a mean of 0 and 
standard deviation equal to the initial vector for the J2000 Rtn direction.

Inputs:
MeanVel - [3x1] mean velocity in J2000 cartesian

mcStart1SigStDevRtn - [3xnTD] matrix where each colunm is a normally
distributed random Rtn perturnation.

nTD - Amount of perturbations desired.

Version:

 - 6/11/18 - Beta - development

%}

trainPertJ2000X = normrnd (0,SigStDevJ2000(1),1,nTD);
trainPertJ2000Y = normrnd (0,SigStDevJ2000(2),1,nTD);
trainPertJ2000Z = normrnd (0,SigStDevJ2000(3),1,nTD);

unitPert = [trainPertJ2000X; trainPertJ2000Y; trainPertJ2000Z];

sigmaPert = unitPert./SigStDevJ2000;

end