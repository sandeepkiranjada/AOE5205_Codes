%% Problem 3
function [A,B] = linearizedmodelgravgradsc01(norbit,Ib11,Ib22,Ib33)
%
%  Copyright (c) 2019 Mark L. Psiaki.  All rights reserved.  
%
%  This function computes the linearized version of the 
%  nonlinear dynamics model of the rigid-body attitude 
%  dynamics of a spacecraft that is subject to the gravity-
%  gradient torques caused by a spherical central attracting 
%  body.  The spacecraft orbits this body in a circular 
%  orbit with mean motion norbit radians/sec and orbital period
%  Torbit = 2*pi/norbit.  The principal moments of inertia
%  are Ib11, Ib22, and Ib33, with Ib11 being the moment-of-
%  inertia about the principal axis that is nominally
%  aligned with the velocity vector (i.e., nominally the
%  roll axis), Ib22 being the moment-of-inertia
%  about the principal axis that nominally points out
%  the "right wing" (i.e., nominally the pitch axis),
%  and Ib33 being the moment-of-inertia about the
%  principal axis that nominally points towards nadir/
%  the center of the Earth (i.e., nominally the yaw axis).  
%  This is a linearization of the nonlinear model that is 
%  contained in ffunctgravgradsc03.m.
%
%  The state vector of the linearized dynamic model has only 6
%  elements despite the corresponding nonlinear model having a
%  7-element state vector.  This model's 6-element state vector is:
%
%     Deltaxtil = [Deltaq1;Deltaq2;Deltaq3;Deltaomegab1;...
%
%                  Deltaomegab2;Deltaomegab3]
%
%  where Deltaq1, Deltaq2, and Deltaq3 are the perturbations
%  of the first three elements of the actual quaternion from
%  the nonlinear system's equilibrium quaternion value
%  qeq = [0;0;0;1] and where Deltaomegab1, Deltaomegab2, and
%  Deltaomegab3 are the perturbations of the components of the
%  actual inertial angular rate along body axes (which are
%  principal axes) from the equilibrium value omegabeq = ...
%  [0;-norbit;0].  Thus, xeq = [0;0;0;1;0;-norbit;0] is the
%  equilibrium state from which perturbations are measured.
%
%  Recall that q = x(1:4,1) in the original
%  nonlinear system state vector is the unit-normalized
%  attitude quaternion for the rotation from local-level
%  orbit-following coordinates to spacecraft body-axes
%  coordinates and that omegab = x(5:7,1) in the original
%  nonlinear system state vector is the spin-rate vector
%  of the body-axis coordinate system relative to inertial
%  coordinates and resolved into components that are defined
%  along the body-fixed axes.
%
%  Note that the control input is the net external torque
%  in addition to the gravity-gradient torque.  It is
%  defined along spacecraft body-fixed axes.  Thus, u = Tb.  Note
%  that the equilibrium value is ueq = Tbeq = [0;0;0].
%  
%
%  Inputs:
%
%    norbit                The mean orbital motion in radians/sec.
%                          Note that the orbital period is Torbit = ...
%                          2*pi/norbit.
%
%    Ib11                  The moment of inertia about the principal
%                          axis that is nominally aligned with
%                          the roll axis (the velocity axis),
%                          in kg-m^2.
%
%    Ib22                  The moment of inertia about the principal
%                          axis that is nominally aligned with
%                          the pitch axis (out the "right wing"),
%                          in kg-m^2.
%
%    Ib33                  The moment of inertia about the principal
%                          axis that is nominally aligned with
%                          the yaw axis (the nadir-pointing axis),
%                          in kg-m^2.
%
%  Outputs:
%
%    A                     The 6-by-6 state coefficient matrix
%                          of the linearized model about xeq and
%                          ueq.
%
%    B                     The 6-by-3 control coefficient matrix
%                          of the linearized model about xeq and
%                          ueq.
%
%                          The linearized dynamics model takes
%                          the form
%
%                            Deltaxtildot(t) = A*Deltaxtil(t) + B*Deltau(t)
%
%                          where Deltaxtil = x([1:3,5:7],1) - ...
%                          xeq([1:3,5:7],1) and Deltau = u - ueq,
%                          with xeq and ueq defined above.
%                          Thus, Deltaxtil has had the 4th element
%                          of Deltax = x - xeq deleted from it
%                          because this fourth element, Deltaq4
%                          is known to equal 0 to first-order
%                          in the linearized perturbations due to
%                          the quaternion unit normalization 
%                          constraint.
%

%
%  Initialize the output arrays.
%
   A = zeros(6,6);
   B = zeros(6,3);
%
%  Assign the individual non-zero elements of these two arrays.
%
   A(1,3) = norbit;
   A(1,4) = 0.5;
   A(2,5) = 0.5;
   A(3,1) = -norbit;
   A(3,6) = 0.5;
   norbitsq = norbit^2;
   sixnorbitsq = 6*norbitsq;
   Iratio_row4 = (Ib33 - Ib22)/Ib11;
   A(4,1) = sixnorbitsq*Iratio_row4;
   A(4,6) = norbit*Iratio_row4;
   Iratio_row5 = (Ib33 - Ib11)/Ib22;
   A(5,2) = sixnorbitsq*Iratio_row5;
   Iratio_row6 = (Ib22 - Ib11)/Ib33;
   A(6,4) = norbit*Iratio_row6;
   B(4,1) = 1/Ib11;
   B(5,2) = 1/Ib22;
   B(6,3) = 1/Ib33; 