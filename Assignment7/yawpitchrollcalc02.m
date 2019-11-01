function [phi,theta,psi] = yawpitchrollcalc02(q)
%
%  Copyright (c) 2017 Mark L. Psiaki.  All rights reserved.  
%
%
%  This function computes the roll, pitch, and yaw angles that
%  yield the same rotation matrix as a quaternion that
%  transforms from local-level orbit-following coordinates,
%  in which ihata points along the velocity vector, khata
%  points down, and jhata = cross(khata,ihata) points 
%  "out the right wing".  The roll, pitch, and yaw
%  angles are defined using a 3-1-2 Euler rotation so that
%
%    R(q) = R2(theta)*R1(phi)*R3(psi)
%
%  with phi being the roll angle, theta being the pitch angle,
%  and psi being the yaw angle.
%
%  Inputs:
%
%    q                       The 4-by-1 unit-normalized quaternion
%                            that parameterizes the transformation
%                            from local-level orbit-following
%                            coordinates to spacecraft
%                            body-fixed coordinates.
%
%  Outputs:
%
%    phi                     The roll angle of the spacecraft,
%                            in radians.
%
%    theta                   The pitch angle of the spacecraft,
%                            in radians.
%
%    psi                     The yaw angle of the spacecraft,
%                            in radians.

%
%  Compute the roll angle.
%
   tanthetanum = q(1,1)*q(3,1) - q(2,1)*q(4,1);
   tanthetaden = 0.5*(q(1,1)^2 + q(2,1)^2 - ...
                        (q(3,1)^2) - (q(4,1)^2));
   theta = atan2(tanthetanum,tanthetaden);
%
%  Determine the qphipsi = qrtmul(qX(phi),qZ(psi))
%
   thetao2 = 0.5*theta;
   costhetao2 = cos(thetao2);
   sinthetao2 = sin(thetao2);
   Qythetatr = ...
       [  costhetao2,            0, sinthetao2,             0;...
                     0, costhetao2,            0, -sinthetao2;...
         -sinthetao2,            0, costhetao2,             0;...
                     0, sinthetao2,            0,  costhetao2];
   qphipsi = Qythetatr*q;
   qphipsi = qphipsi*(1/sqrt(sum(qphipsi.^2)));
%
%  Use the elements of qphipsi to compute phi and psi.
%
   psi_raw = 2*atan2(qphipsi(3,1),qphipsi(4,1));
   twopi = 2*pi;
   ootwopi = 1/twopi;
   psi = psi_raw - twopi*round(psi_raw*ootwopi);
%
   psi_rawo2 = 0.5*psi_raw;
   cospsi_rawo2 = cos(psi_rawo2);
   sinpsi_rawo2 = sin(psi_rawo2);
   sinphio2 = qphipsi(1,1)*cospsi_rawo2 + ...
                     qphipsi(2,1)*sinpsi_rawo2;
   sinphio2 = max([-1;min([1;sinphio2])]);
   phi = 2*asin(sinphio2);
%
%  Double check that the angles produce the correct quaternion
%
   psio2 = 0.5*psi;
   qzpsiA = [0;0;sin(psio2);cos(psio2)];
   phio2 = 0.5*phi;
   qxphiA = [sin(phio2);0;0;cos(phio2)];
   qythetaA = [0;sinthetao2;0;costhetao2];
   qtot = qrtmul(qythetaA,qrtmul(qxphiA,qzpsiA));
   if sum(qtot.*q) < 0
      qtot = - qtot;
   end
   if norm(qtot - q) > 1.e-09
      disp('Warning in yawpitchrollcalc02.m: First roll/pitch/yaw')
      disp(' angular solution does not reproduce q to 9')
      disp(' significant digits.')
   end
%
%  Try alternate roll, pitch, and yaw angles.
%
   theta_alt = theta + pi;
   theta_alt = theta_alt - twopi*round(theta_alt*ootwopi);
%
   theta_alto2 = 0.5*theta_alt;
   costheta_alto2 = cos(theta_alto2);
   sintheta_alto2 = sin(theta_alto2);
   Qythetatr_alt = ...
       [  costheta_alto2,                0, sintheta_alto2, ...
                                                  0;...
                         0, costheta_alto2,                0, ...
                                  -sintheta_alto2;...
         -sintheta_alto2,                0, costheta_alto2,      ...
                                                  0;...
                         0, sintheta_alto2,                0, ...
                                   costheta_alto2];
   qphipsi_alt = Qythetatr_alt*q;
   qphipsi_alt = qphipsi_alt*(1/sqrt(sum(qphipsi_alt.^2)));
%
%  Use the elements of qphipsi to compute phi and psi.
%
   psi_raw_alt = 2*atan2(qphipsi_alt(3,1),qphipsi_alt(4,1));
   psi_alt = psi_raw_alt - twopi*round(psi_raw_alt*ootwopi);
%
   psi_rawo2_alt = 0.5*psi_raw_alt;
   cospsi_rawo2_alt = cos(psi_rawo2_alt);
   sinpsi_rawo2_alt = sin(psi_rawo2_alt);
   sinphio2_alt = qphipsi_alt(1,1)*cospsi_rawo2_alt + ...
                         qphipsi_alt(2,1)*sinpsi_rawo2_alt;
   sinphio2_alt = max([-1;min([1;sinphio2_alt])]);
   phi_alt = 2*asin(sinphio2_alt);
%
%  Double check that the angles produce the correct quaternion
%
   psio2_alt = 0.5*psi_alt;
   qzpsiA = [0;0;sin(psio2_alt);cos(psio2_alt)];
   phio2_alt = 0.5*phi_alt;
   qxphiA = [sin(phio2_alt);0;0;cos(phio2_alt)];
   qythetaA = [0;sintheta_alto2;0;costheta_alto2];
   qtot_alt = qrtmul(qythetaA,qrtmul(qxphiA,qzpsiA));
   if sum(qtot_alt.*q) < 0
      qtot_alt = - qtot_alt;
   end
   if norm(qtot_alt - q) > 1.e-09
      disp('Warning in yawpitchrollcalc02.m: Second roll/pitch/yaw')
      disp(' angular solution does not reproduce q to 9')
      disp(' significant digits.')
   end
%
%  Choose as the output angles those that yield the smallest
%  absolute value of the pitch angle.
%
   if abs(phi_alt) < abs(phi)
      theta = theta_alt;
      phi = phi_alt;
      psi = psi_alt;
   end