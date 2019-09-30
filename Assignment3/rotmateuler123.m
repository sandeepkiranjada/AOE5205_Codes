function R = rotmateuler123(psi,theta,phi)
%
%  Copyright (c) 2019 Mark L. Psiaki.  All rights reserved.  
%
%  This function computes the 1-2-3 Euler-angles 3-by-3 direction
%  cosines rotation matrix for the transformation from
%  local-level north-east-down coordinates to aircraft
%  body-frame coordinates:
%
%    R = R3(psi)*R2(theta)*R1(phi)
%
%  where Rm is the simple rotation matrix about the mth
%  coordinate axis, psi is the yaw angle, theta is the
%  pitch angle, and phi is the roll angle.
%  
%
%  Inputs:
%
%    psi                   The scalar yaw angle, in radians.
%
%    theta                 The scalar pitch angle, in radians.
%
%    phi                   The scalar roll angle, in radians.
%
%  Outputs:
%
%    R                     The 3-by-3 direction cosines matrix
%                          for the rotation transformation from
%                          north-east-down coordinates to 
%                          body cordinates.  It is non-dimensional.
%

%
%  Compute the yaw rotation matrix about the khat axis.
%
   R3atpsi = [ cos(psi) sin(psi) 0 ;...
              -sin(psi) cos(psi) 0 ;
                    0        0   1];
%
%  Compute the pitch rotation matrix about the jhat axis.
%
   R2attheta = [ cos(theta) 0 -sin(theta) ;...
                       0    1        0    ;
                 sin(theta) 0  cos(theta)];
%
%  Compute the roll rotation matrix about the ihat axis.
%
   R1atphi = [1       0        0   ;...
              0  cos(phi) sin(phi) ;...
              0 -sin(phi) cos(phi)];              
%
%  Form the total rotation matrix as a sequence of
%  the three rotations:
%
   R = R3atpsi*R2attheta*R1atphi;