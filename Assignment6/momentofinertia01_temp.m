function [Mtot,rcmtot,IMoItot] = momentofinertia01(mvec,rcmmat,IMoIarray)
%
%  Copyright (c) 2019 Mark L. Psiaki.  All rights reserved.  
%
%  This function computes the total mass, the center of mass,
%  and the total moment of inertia of a collection of
%  rigid bodies that, taken together, form a larger rigid body.
%
%  All position vectors, those of the individual rigid bodies'
%  centers of mass, rcmi = rcmmat(:,i) for i = 1:N, and that of 
%  the system center of mass, rcmtot, are given in a common 
%  coordinate system as are the individual moment-of-inertia 
%  matrices, IMoIi = IMoIarray(:,:,i), and the final total 
%  system moment-of-inertia matrix, IMoItot. 
%  
%
%  Inputs:
%
%    mvec                  The 1-by-N vector that contains the
%                          masses of the individual rigid-body
%                          components, in kg units.  mi = mvec(1,i)
%                          is the mass of the ith rigid body.
%
%    rcmmat                The 3-by-N matrix that contains the
%                          positions of the centers of mass
%                          of the individual rigid bodies,
%                          given in meters units and along the
%                          common axes that are used 
%                          throughout these calculations.
%                          rcmi = rcmmat(:,i) is the center-
%                          of-mass position of the ith 
%                          rigid body.
%
%    IMoIarray             The 3-by-3-by-N array that contains
%                          the moment-of-inertia matrices of the
%                          individual rigid bodies about their
%                          respective centers of mass, in
%                          kg-m^2 units and along the common 
%                          axes that are used throughout these 
%                          calculations. IMoIi =  IMoIarray(:,:,i)
%                          is the moment-of-inertia matrix of
%                          the ith rigid body about its own
%                          center of mass.
%
%  Outputs:
%
%    Mtot                  The total mass of the composite
%                          rigid body, in kg.
%
%    rcmtot                The 3-by-1 vector that gives the
%                          center of mass of the composite rigid
%                          body, in meters and along the common 
%                          axes that are used throughout these 
%                          calculations. 
%
%    IMoItot               The 3-by-3 moment-of-inertia matrix
%                          of the composite rigid body about its
%                          center of mass, in kg-m^2 and along 
%                          the common axes that are used 
%                          throughout these calculations.
%

%
%  Compute the total mass.
%
   Mtot = ????;
%
%  Compute the composite rigid body's center of mass.
%
   N = size(mvec,2);
   Mtot_rcmtot = zeros(3,1);
   for i = 1:N
      mi = mvec(1,i);
      rcmi = rcmmat(:,i);
      Mtot_rcmtot = Mtot_rcmtot + ????;
   end
   rcmtot = ????;
%
%  Compute the composite rigid body's moment-of-inertia
%  matrix about its center of mass.
%
   IMoItot = zeros(3,3);
   for i = 1:N
      mi = mvec(1,i);
      rcmi = rcmmat(:,i);
      deltarcmi = rcmi - rcmtot;
      IMoIi =  IMoIarray(:,:,i);
      deltaIMoIi = ????;
      IMoItot = ????;
   end