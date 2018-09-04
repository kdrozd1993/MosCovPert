% Approach Perturbation Validation Script

clear all
clc

%% Obtaining perturbation in km from kernels

% Loading frame and leap second kernel
cspice_furnsh ('orx_v09.tf');
cspice_furnsh ('naif0012.tls');
cspice_furnsh('bennu_refdrmc_v1.bsp');

% Specifying UTC time and putting it in ephemeris time
utc = '12 Nov 2018 18:00:00';
et = cspice_str2et (utc);

% Loading in 2-sigma perturbed spk in the radial direction 
cspice_furnsh ('orx_181112_181212_pert_wrt_MPRevB_v2.bsp');

% Obtaining position in RTN frame of perturbed trajectory at time specified
stateP = cspice_spkpos( '-64', et, 'J2000', ...
                                        'NONE', 'Bennu' )

% Unloading perturbed spk                                    
cspice_unload ('orx_181112_181212_pert_wrt_MPRevB_v2.bsp');

% Loading reference spk
cspice_furnsh ('orx_170117_210312_refMPRevB_v1.bsp');

% Obtaining position in RTN frame of reference trajectory at time specified
stateRef = cspice_spkpos( '-64', et, 'J2000', ...
                                        'NONE', 'Bennu' )
                      
rotMatJ20002RtnStart = cspice_pxform ('J2000','ORX_RTN_BENNU',et);                                    

statePRtn = rotMatJ20002RtnStart*stateP
stateRefRtn = rotMatJ20002RtnStart*stateRef

kernelPerturbation = statePRtn-stateRefRtn                                   

%% Attempting to get same perturbation, but from Monte Carlo Analysis file
% Error and ref taken from then monte carlo file are in spherical SAM
% coordinates taken on the utc date.

% Error from Monte Carlo file (1-sig) [range (km), lat (deg), long(deg)]
error = [3.18; 1.510; 1.338];

% Ref from Monte Carlo file (1-sig) [range (km), lat (deg), long(deg)]
ref = mcSph2SamCoords ([63.883; 4.565; 8.378]); 

% Used 3 equations with 3 unkowns to obtain the 1-sigma standard deviation 
% errors in cartesian SAM coordinates from spherical SAM coordinates. See
% picture for how these equations were derived. Did it this way because I
% wasn't convinced the typical conversion between spherical and cartesian
% would work with standard deviation. I could be wrong.
sigmax = sqrt(error(1)^2/(1+tand(error(3))^2+tand(error(2))^2+tand(error(3))^2+tand(error(2))^2));
sigmay = tand(error(3))*sigmax;
sigmaz = tand(error(2))*sqrt(sigmax^2+sigmay^2);

% Getting the rotation matrix to go from SAM to RTN frame
rotMatSam2RtnStart = cspice_pxform ('ORX_BENNU_SAM','ORX_RTN_BENNU',et);

% Reference RTN vector at the utc from monte carlo file
refRTN = rotMatSam2RtnStart*ref;

% 1-sig error RTN vector at the utc from monte carlo file
errorRTN = rotMatSam2RtnStart*[sigmax;sigmay;sigmaz];

errorSam = [3.2094;1.4485;1.6555];
ey = rotMatSam2RtnStart*errorSam;

%%

% Unloading spk                                    
%cspice_unload ('orx_171006_190608_180125_od030-R-DSM2-P-M28D_v1.bsp');

% Loading reference spk
cspice_furnsh ('orx_170117_210312_refMPRevB_v1.bsp');                               
                                    
                                    
errorRTNreal = [3.18*2;1.537*0;1.66*0];
refFromMC = [63.883;0;0];

ha = rotMatJ20002RtnStart'*refFromMC;
ho = rotMatJ20002RtnStart'*errorRTNreal;

yo = ha+ho;

rotMatJ20002RtnStart*yo;