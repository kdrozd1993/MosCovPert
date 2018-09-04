%-------------------------------------------------------------------------%
function [rotMat] = RotMatJ2000_2_Obs (ScPos,sunPos)
%-------------------------------------------------------------------------%

%{
Author:
Kristofer Drozd

Description:
Obtains the rotation matrix from J2000 to Obs frame.

Inputs:
ScPos - [3x1] vector that is the position in J2000 of the body of interest 
        (Bennu) from the SC.

sunPos - [3x1] vector that is the position in J2000 of the Sun from the SC.

Version:

 - 6/11/18 - Beta - development

%}

StartObs_Z = ScPos/norm (ScPos);
Cross_Y = cross (ScPos,sunPos);
StartObs_Y = Cross_Y/norm (Cross_Y);
StartObs_X = cross (StartObs_Y,StartObs_Z);

rotMat = [StartObs_X'; StartObs_Y'; StartObs_Z'];


end