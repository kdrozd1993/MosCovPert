%-------------------------------------------------------------------------%
function [label] = UeSizeLabeling (PertUeSize, RefUeSize)
%-------------------------------------------------------------------------%

%{
Author:
Kristofer Drozd

Description:
Obtains the bins for UE size change due to state perturbation

Inputs:
PertUeSize - Scalar representing Uncertainty Ellipse component size on the
             perturbed trajectory

RefUeSize - Scalar representing the Uncertainty Ellipse component size on
            the reference trajectory

Outputs:
label - Scalar placing the UeSize variable in a bin. Bins represent that
        the UE size change is above, below, or above and below a certain 
        number

Version:

 - 08/24/18 - Code Completed

%}

if PertUeSize <= (RefUeSize - RefUeSize * 0.1)
    
    label = 0;
    
elseif PertUeSize > (RefUeSize - RefUeSize * 0.1) ...
        && PertUeSize <= RefUeSize
    
    label = 1;
    
elseif PertUeSize > RefUeSize ...
        && PertUeSize <= (RefUeSize + RefUeSize * 0.1)
    
    label = 2;
       
else
    
    label = 3;
    
end