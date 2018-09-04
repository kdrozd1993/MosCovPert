%-------------------------------------------------------------------------%
function [label] = MosaicCoverageLabeling (MosCov)
%-------------------------------------------------------------------------%

%{
Author:
Kristofer Drozd

Description:
Obtains the bins for mosaic coverage. Bins are below 1, above 1 and below
2, above 2 and below 3, and above 3.

Inputs:
MosCov - Scalar representing mosaic coverage in sigma

Outputs:
label - Scalar placing the MosCov variable in a bin. Bins represent that
        the coverage is above, below, or above and below a certain number

Version:

 - 08/02/18 - Code Completed

%}

% if MosCov < 1
%     
%     label = 0;
%     
% elseif MosCov >= 1 && MosCov < 1.2
%     
%     label = 1;
%     
% elseif MosCov >= 1.2 && MosCov < 1.4
%     
%     label = 1.2;
%     
% elseif MosCov >= 1.4 && MosCov < 1.6
%     
%     label = 1.4;
% 
% elseif MosCov >= 1.6 && MosCov < 1.8
%     
%     label = 1.6;
% 
% elseif MosCov >= 1.8 && MosCov < 2.0
%     
%     label = 1.8;
%     
% elseif MosCov >= 2.0 && MosCov < 2.2
%     
%     label = 2.0;
%     
% elseif MosCov >= 2.2 && MosCov < 2.4
%     
%     label = 2.2;
% 
% elseif MosCov >= 2.4 && MosCov < 2.6
%     
%     label = 2.4;
% 
% elseif MosCov >= 2.6 && MosCov < 2.8
%     
%     label = 2.6;
%     
% elseif MosCov >= 2.8 && MosCov < 3.0
%     
%     label = 2.8;    
%     
% else
%     
%     label = 3;
%     
% end

if MosCov < 2
    
    label = 0;
    
elseif MosCov >= 2 && MosCov < 2.5
    
    label = 2;
    
elseif MosCov >= 2.5 && MosCov < 3.0
    
    label = 2.5;
        
else
    
    label = 3.0;
    
end
