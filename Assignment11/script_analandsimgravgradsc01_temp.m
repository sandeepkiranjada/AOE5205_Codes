%script_analandsimgravgradsc01.m
%
%  Copyright (c) 2019 Mark L. Psiaki.  All rights reserved.  
%
%  This Matlab script performs various tasks
%  related to analysis of spacecraft attitude
%  dynamics acting under the influence of gravity
%  gradient torque and additional external
%  torque in the neighborhood of an equilibrium,
%  to design of an observer for this motion
%  using a linearized system model, and to do linear
%  and nonlinear simulation of the motion of the
%  resulting system and the observer.
%
%  This script works with the attitude dynamics model and
%  the steady-motion conditions that are given in the
%  file linearizedmodelgravgradsc13_data.mat. 
%
%  This script also makes plots of attitude time histories.
%
%  Clear the Matlab workspace and set up for long-format
%  Matlab display.
%
   clear
   format long
%
%  Load the spacecraft parameters and the mean orbital
%  motion.  Note that an
%  alternate data file can be loaded in order to yield
%  test data for which displayed answers are given.
%
   load linearizedmodelgravgradsc13_data
%  load linearizedmodelgravgradsc16_data   
%
%  Compute the linearized model's A and B matrices and
%  display them.
%
   [A,B] = linearizedmodelgravgradsc01(norbit,Ib11,Ib22,Ib33)
%
%  Compute the eigenvalues of the A matrix in order
%  to test stability.  Let the results be displayed
%  after sorting them in ascending order of their absolute
%  values.
%
   lambdavec = eig(A);
   [~,idumvec] = sort(abs(lambdavec));
   lambdavec = lambdavec(idumvec,1)
   clear dum idumvec
%
%  Determine whether the linearized open-loop system stable, 
%  unstable, or neutrally stable.
%
   maxreallambda = max(real(lambdavec))
   lambdatestval = max(abs(lambdavec))*1.e-10;
   iflagneutralstabtest = 0;
   if maxreallambda < -lambdatestval
      disp('This system appears to be stable because all of its')
      disp(' eigenvalues have negative real parts.')
   elseif maxreallambda > lambdatestval
      disp('This system appears to be unstable because one or more')
      disp(' of its eigenvalues has a positive real part.')
   else
      disp('This system may be neutrally stable because the maximum')
      disp(' eigenvalue real part maximized over all of its eigenvalues')
      disp(' appears to be zero to within machine precision.')
      iflagneutralstabtest = 1;
   end
%
%  Find all repeated eigenvalues that have
%  zero real parts to within approximately machine
%  precision.
%
   if iflagneutralstabtest == 1
      irealzerosvec = find(abs(real(lambdavec)) <= lambdatestval);
      nrealzerosvec = size(irealzerosvec,1);
      lambdavec_realzerovec = lambdavec(irealzerosvec,1);
      lambdavec_realzerovecunique = lambdavec_realzerovec(1,1);
      nrepeat_realzerovecunique = 1;
      lambdavec_realzerovecrem = lambdavec_realzerovec;
      lambdavec_realzerovecrem(1,:) = [];
      nleft = nrealzerosvec - 1;
      nunique = 1;
      while nleft > 0
         testvec = abs(lambdavec_realzerovecrem - ...
                       lambdavec_realzerovecunique(nunique,1));
         iequalvec = find(testvec <= lambdatestval);
         nequalvec = size(iequalvec,1);
         if nequalvec > 0
            nrepeat_realzerovecunique(nunique,1) = ...
                nrepeat_realzerovecunique(nunique,1) + nequalvec;
            lambdavec_realzerovecrem(iequalvec,:) = [];
            nleft = nleft - nequalvec;
         end
         if nleft > 0
            nunique = nunique + 1;
            lambdavec_realzerovecunique = ...
                  [lambdavec_realzerovecunique;...
                   lambdavec_realzerovecrem(1,1)];
            nrepeat_realzerovecunique = ...
                  [nrepeat_realzerovecunique;1];
            lambdavec_realzerovecrem(1,:) = [];
            nleft = nleft - 1;
         end
      end
   else
      lambdavec_realzerovecunique = [];
      nrepeat_realzerovecunique = [];
      nunique = 0;
   end
   clear irealzerosvec nrealzerosvec lambdavec_realzerovec ...
         lambdavec_realzerovecrem nleft testvec iequalvec ...
         nequalvec iflagneutralstabtest lambdatestval
%
%  Test all repeated eigenvalues with real parts to check
%  whether they do not admit neutral stability.
%
   n = size(A,1);
   if nunique > 0
      for jj = 1:nunique
         nrepeatjj = nrepeat_realzerovecunique(jj,1);
         if nrepeatjj > 1
            lambdajj = lambdavec_realzerovecunique(jj,1);
            rankjj = rank(eye(n)*lambdajj - A);
            rankjj_neutralstab = n - nrepeat_realzerovecunique(jj,1);
            if rankjj > rankjj_neutralstab
               disp(' ')
               disp(['Warning: For eigenvalue lambda = ',...
                      num2str(lambdajj),...
                      ', the rank of (lambda*eye(n) - A) is'])
               disp([' ',int2str(rankjj),' but it should be smaller,',...
                     ' it should be ',int2str(rankjj_neutralstab),' in',...
                     ' order for neutral'])
               disp([' stability to hold true, because this eigenvalue',...
                     ' is repeated ',int2str(nrepeatjj),' times.'])
               disp(' Therefore, this system is unstable.')
            end
         end
      end
   end
   clear jj lambdajj rankjj nrepeatjj rankjj_neutralstab
%
%  Define a linearized system output matrix where
%  the first output is the measured roll angle (actually,
%  twice the roll component of the attitude quaternion,
%  which is equivalent for a linearized model),
%  the second output is the measured pitch angle (actually,
%  twice the pitch component of the attitude quaternion,
%  which is equivalent for a linearized model),
%  and the third outout is the measured yaw rate
%  relative to inertial coordinates.  The first
%  two outputs could be measured by a horizon
%  scannere.  The third output could be measured
%  by a yaw axis rate-gyro.
%
   C = [ 2, 0, 0, 0, 0, 0;...
         0, 2, 0, 0, 0, 0;...
         0, 0, 0, 0, 0, 1];
%
%  Check system observability.  Compute the observability
%  matrix and compute its rank.  Also, compute and display 
%  the singular values.  The number of non-zero singular
%  values equals the rank of the matrix.  Therefore,
%  they should all be non-zero in order for the system
%  to be controllable.
%  
   Observabilitymat = ????;
   nrank_Observabilitymat = rank(Observabilitymat);
   svsObservabilitymat = svd(Observabilitymat)
   if nrank_Observabilitymat == n
      disp(' ')
      disp('The system is observable.')
   else
      disp(' ')
      disp('The system is not observable.  The rank of the')
      disp([' observability matrix is only ',...
            int2str(nrank_Observabilitymat),'.'])
   end
   clear nrank_Observabilitymat
%
%  Design a full-state observer for this sytem
%  using the place.m function.  Design the
%  observer so that the linearized error dynamics
%  exhibit 3 distinct pairs of complex-conjugate
%  eigenvalues, each with a damping ratio of 1/sqrt(2).
%  The first pair wil have an undamped natural 
%  frequency of 0.0035 radians/sec, the second pair
%  0.0060 radians/sec, and the third pair 0.0085 radians/sec.
%  Display the observer gain.
% 
   omegana = 0.0035;
   omeganb = 0.0060;
   omeganc = 0.0085;
   oosqrt2 = 1/sqrt(2);
   observereigenvalues = oosqrt2*...
         [(-omegana + sqrt(-1)*omegana);...
          (-omegana - sqrt(-1)*omegana);...
          (-omeganb + sqrt(-1)*omeganb);...
          (-omeganb - sqrt(-1)*omeganb);...
          (-omeganc + sqrt(-1)*omeganc);...
          (-omeganc - sqrt(-1)*omeganc)];
   L = ????
%
%  Prepare to do a nonlinear simulation of the 
%  response of this system and its observer 
%  to a perturbed initial condition.
%
%  This is the equilibrium state and control.
%
   xeq = [0;0;0;1;0;(-norbit);0];
   ueq = zeros(3,1);
%
%  This is the perturbed initial condition.
%  The alternate perturbed initial condition reduces
%  the perturbation by a factor 10 as a test case 
%  to find out whether agreement between the linear
%  and nonlinear models improves for a smaller
%  perturbation from the steady-motion trajectory.
%
   Deltax0 = [-(2*pi/180);(3*pi/180);(-0.5*pi/180);0;...
               (-0.005*norbit);(0.005*norbit);(-0.001*norbit)];
%  Deltax0 = Deltax0*(1/10);        
   x0 = xeq + Deltax0;
%
%  Normalize the attitude quaterion part of the
%  perturbed initial state.
%
   x0(1:4,1) = x0(1:4,1)*(1/sqrt(sum(x0(1:4,1).^2)));
%
%  Extend the measurement model matrix and the 
%  observer gain to include the 4th quaternion state.
%
   Cextended = [C(:,1:3),zeros(3,1),C(:,4:6)];
   Lextended = [L(1:3,:);zeros(1,3);L(4:6,:)];
%
%  Define the body-axis moment-of-inertia matrix.
%
   IMoIbody = diag([Ib11;Ib22;Ib33]);
%
%  Define the nonlinear dynamics model for the observer.
%  This model effectively takes the form:
%
%    xobsdot = f(xbos,u) + L*(C*(x - xobs))
%
%  There are added complications not found in standard
%  observers that are caused by the quaternion's 4th
%  state.  Some of the complication comes from the fact
%  that the 4th quaternion element does not enter the 
%  linearized model.  That is why Lextended and
%  Cextended are used in the model below rather
%  than L and C.  Additionally, there is
%  a need to enforce the quaternion's unit
%  normalization in the observer by ensuring
%  that the time rate of change of the observer
%  quaternion's magnitude is zero.  Note the 
%  matrix involving (eye(4) - xhatarga(1:4,1)*...
%  (xhatarga(1:4,1)')) that pre-multiplies the observer
%  feedback term that is part of the quaternion
%  rate model, i.e., part of the first 4 elements of
%  the function ffunctobserver.  This term
%  enforces the constraint that the time rate
%  of change of the quaternion magnitude must
%  be zero.  It has no impact to first order
%  in small perturbations from the equilibrium,
%  but for larger perturbations it will modify
%  the nonlinear dynamics of the observer.
%
   ffunctobserver = @(xhatarga,xarga) ...
          ffunctgravgradsc03(xhatarga,ueq,IMoIbody,norbit) + ...
              [(eye(4) - xhatarga(1:4,1)*(xhatarga(1:4,1)')),...
                 zeros(4,3);...
               zeros(3,4),eye(3)]*...
                  (Lextended*(Cextended*(xarga - xhatarga)));
%
%  Define the nonlinear model for the full system that
%  consists of the true system state as its first
%  7 elements and the observer state is its last 7 elements.
%  This function is needed in order to simulate the
%  observer because its dynamics model must be integrated
%  simultaneously with that of the original system.
%  In effect, this model stacks the 7 states of the
%  true spacecraft together with the 7 states of the
%  nonlinear observer into a 14-dimensional state
%  vector, and it computes the dynamics model
%  for it.
%
   ffunctaugmented = @(targb,xaugmentedargb) ...
     [ffunctgravgradsc03(xaugmentedargb(1:7,1),ueq,IMoIbody,norbit);...
      ffunctobserver(xaugmentedargb(8:14,1),xaugmentedargb(1:7,1))];
%
%  Set up the inital state for the full system, with the
%  observer starting at the equilibrium state.
%
   xaugmented0 = [x0;xeq];
%
%  Set up the remaining inputs to ode45.m.
%
   Torbit = 2*pi/norbit;
   tspan = ((0:1000)')*(5*Torbit/1000);      
   optionsode45 = odeset('RelTol',1.e-10);
%
%  Run ode45.m.  Afterwards, compute the perturbation from the
%  nominal steady-motion trajectory and the control time history,
%  both absolute and perturbed.
%
   [thist,xaugmentedhist] = ode45(ffunctaugmented,tspan,...
                                  xaugmented0,optionsode45);
   Deltaxerrhist = xaugmentedhist(:,[8 9 10 12 13 14]) - ...
                   xaugmentedhist(:,[1 2 3 5 6 7]);
%
%  Redo the simulation using the linearized observer error
%  dynamics model form
%  loop linearized model takes the form:
%
%    Deltaxerrdot = Aclobs*Deltaxerr
%
   Aclobs = ????;
   ffunctobserrlin = @(targdumc,xargdumc) Aclobs*xargdumc;
   Deltaxerr0 = zeros(6,1) - Deltax0([1 2 3 5 6 7],1);
   [thist_alt,Deltaxerrhist_alt] = ...
          ode45(ffunctobserrlin,tspan,Deltaxerr0,optionsode45);
%
%  Plot the true attitude time history.
%
   figure(1)
   hold off
   plot(thist/Torbit,xaugmentedhist(:,1:3)*2*(180/pi),'LineWidth',3)
   grid
   xlabel('Time (orbits)')
   ylabel('Attitude Angle (deg)')
   title('True attitude response from analandsimgravgradsc01.mat')
   legend('Roll','Pitch','Yaw')
%
%  Plot nonlinear and linearized observer attitude error
%  time histories.
%
   figure(2)
   hold off
   plot(thist/Torbit,Deltaxerrhist(:,1:3)*2*(180/pi),'LineWidth',3)
   hold on
   plot(thist_alt/Torbit,Deltaxerrhist_alt(:,1:3)*2*(180/pi),':',...
        'LineWidth',2)
   hold off
   grid
   xlabel('Time (orbits)')
   ylabel('Attitude Angle Estimation Error (deg)')
   title(['Attitude estimation error response from',...
          ' analandsimgravgradsc01.mat'])
   legend('Roll, nonlinear','Pitch, nonlinear',...
          'Yaw, nonlinear','Roll, linearized',...
          'Pitch, linearized','Yaw, linearized')
   xlim([0 1])      
%
%  Save the results.
%
   textcommands = ['These data have been generated by the',...
                   ' commands in script_analandsimgravgradsc01.m'];
   save analandsimgravgradsc01 
%
%
%  Displayed outputs when the alternate input file
%  linearizedmodelgravgradsc16_data.mat is used to
%  compute the model, analyze the model, design
%  an observer, simulate the system and
%  observer response, and plot results for this
%  file's different steady-motion condition:
%
%  A =
%  
%    Columns 1 through 3
%  
%                     0                   0   0.001000000000000
%                     0                   0                   0
%    -0.001000000000000                   0                   0
%    -0.000005684210526                   0                   0
%                     0  -0.000004090909091                   0
%                     0                   0                   0
%  
%    Columns 4 through 6
%  
%     0.500000000000000                   0                   0
%                     0   0.500000000000000                   0
%                     0                   0   0.500000000000000
%                     0                   0  -0.000947368421053
%                     0                   0                   0
%     0.000750000000000                   0                   0
%  
%  B =
%  
%                     0                   0                   0
%                     0                   0                   0
%                     0                   0                   0
%     0.010526315789474                   0                   0
%                     0   0.009090909090909                   0
%                     0                   0   0.050000000000000
%  
%  lambdavec =
%  
%   -0.000000000000000 + 0.000864158930789i
%   -0.000000000000000 - 0.000864158930789i
%    0.000000000000000 + 0.001430193883868i
%    0.000000000000000 - 0.001430193883868i
%   -0.000000000000000 + 0.001950861584348i
%   -0.000000000000000 - 0.001950861584348i
%  
%  maxreallambda =
%  
%       0
%  
%  This system may be neutrally stable because the maximum
%   eigenvalue real part maximized over all of its eigenvalues
%   appears to be zero to within machine precision.
%  
%  svsObservabilitymat =
%  
%     2.000000000014762
%     2.000000000004184
%     1.000002281253503
%     1.000000001385294
%     1.000000000002092
%     0.000001499996580
%   
%  The system is observable.
%  
%  L =
%  
%     0.005763471178738   0.000233029506572  -0.002251544799168
%    -0.000832712644307   0.004477259670438  -0.021283755332446
%     0.006156867617960   0.000251034124737  -7.733846556400598
%     0.000050191522923  -0.000005909101133   0.015434377892544
%    -0.000005271476577   0.000037306408243   0.000021484987600
%     0.000004650584426   0.000002736458145   0.004974382424362