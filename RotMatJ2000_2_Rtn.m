%-------------------------------------------------------------------------%
function [rotMat] = RotMatJ2000_2_Rtn (pos,vel)
%-------------------------------------------------------------------------%

%{
Author:
Kristofer Drozd

Description:
Obtains the rotation matrix from J2000 to Rtn frame.

Inputs:
pos - [3x1] vector that is the position in J2000.

vel - [3x1] vector that is the velocity in J2000.

Version:

 - 6/11/18 - Beta - development

%}

uR = pos/norm (pos);
uCrossN = cross (pos,vel);
uN = uCrossN/norm (uCrossN);
uT = cross (uN,uR);

rotMat = [uR'; uT'; uN'];

end