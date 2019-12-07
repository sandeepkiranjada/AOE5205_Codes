%script_analandsimaircraft01.m
%
%  Copyright (c) 2019 Mark L. Psiaki.  All rights reserved.  
%
%  This Matlab script performs various tasks
%  related to analysis of aircraft steady motion
%  design of a feedback controller for this motion
%  using a linearized system model, and linear
%  and nonlinear simulation of the motion of the
%  resulting system.
%
%  This script works with the airplane model and
%  the steady-motion conditions that are given in the
%  file linearizedmodelaircraft01_data03.mat. 
%
%  This script also makes plots of flight time histories.
%
%  Clear the Matlab workspace and set up for long-format
%  Matlab display.
%
   clear
   format long
%
%  Load the aircraft parameters and the equilibrium state
%  and control vectors that apply to this analysis.  Note that
%  alternate data file can be loaded in order to yield
%  test data for which displayed answers are given.
%
   load linearizedmodelaircraft01_data03
%  load linearizedmodelaircraft01_data04   
%
%  Compute the linearized model's A and B matrices and
%  display them.
%
   [A,B] = linearizedmodelaircraft01(xeq,ueq,m,S,CLalpha,...
                                     CD0,oneoverpiARe)
%
%  Compute the eigenvalues of the A matrix in order
%  to test stability.  Let the results be displayed
%  after sorting them in ascending order of their real
%  parts.
%
   lambdavec = eig(A);
   [~,idumvec] = sort(real(lambdavec));
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
%  Check system controllability.  Compute the controllability
%  matrix and compute its rank.  Also, compute and display 
%  the singular values.  The number of non-zero singular
%  values equals the rank of the matrix.  Therefore,
%  they should all be non-zero in order for the system
%  to be controllable.
%  
   Controllabilitymat = ctrb(A,B);
   nrank_Controllabilitymat = rank(Controllabilitymat);
   svsControllabilitymat = svd(Controllabilitymat)
   if nrank_Controllabilitymat == n
      disp(' ')
      disp('The system is controllable.')
   else
      disp(' ')
      disp('The system is not controllable.  The rank of the')
      disp([' controllability matrix is only ',...
            int2str(nrank_Controllabilitymat),'.'])
   end
   clear nrank_Controllabilitymat
%
%  Design a full-state feedback controller for this sytem
%  using the place.m function.  Have one set of 3 eigenvalues
%  with an undampled natural frequency of 0.02 radians/second
%  with one real eigenvalue at negative this value and two 
%  complex-conjugate eigenvalues using this undamped
%  natural frequency and a damping ratio of 0.5.  Similarly
%  implement another set of 3 eigenvalues with
%  a similar pattern, except using an undampled natural
%  frequency of 0.4 radians/second.  Display the gain.
% 
   omegana = 0.02;
   omeganb = 0.4;
   closedloopeigenvalues = ...
         [(-omegana);...
          (-0.5*omegana + sqrt(-1)*omegana*sqrt(1 - 0.5^2));...
          (-0.5*omegana - sqrt(-1)*omegana*sqrt(1 - 0.5^2));...
          (-omeganb);...
          (-0.5*omeganb + sqrt(-1)*omeganb*sqrt(1 - 0.5^2));...
          (-0.5*omeganb - sqrt(-1)*omeganb*sqrt(1 - 0.5^2))];
   K = ????
%
%  Prepare to do a nonlinear simulation of the closed-loop
%  response of this system to a perturbed initial condition.
%
%  These are the steady-motion rates of the northward
%  travel and of the eastward travel.
%
   XdotSM = xeq(4,1)*cos(xeq(5,1))*cos(xeq(6,1));
   YdotSM = xeq(4,1)*cos(xeq(5,1))*sin(xeq(6,1));
%
%  This is the perturbed initial condition.
%  The alternate perturbed initial condition reduces
%  the perturbation by a factor 20 as a test case 
%  to find out whether agreement between the linear
%  and nonlinear models improves for a smaller
%  perturbation from the steady-motion trajectory.
%
    x0 = xeq + [500;-700;-20;4;-(1*pi/180);(2*pi/180)];
%   x0 = xeq + [500;-700;-20;4;-(1*pi/180);(2*pi/180)]*(1/20);
%
%  This is the linear form of the feedback control law.
%  It computes the perturbation of the state in xargdum from the
%  steady-motion state.  It multiplies this state perturbation 
%  by the feedback control gain.  It then adds the result to the
%  equilibrium control input.
%
   ulinearfeedbackfunct = @(targdum,xargdum) ...
         ???? - ????*(xargdum - xeq - [XdotSM;YdotSM;zeros(4,1)]*targdum);
%
%  This is the full nonlinear form of the feedback control law.
%  It implements a minimum thrust limit
%  of 200 N, a maximum thrust limit of 22000N, a minimum angle
%  of attack limit of -3 deg, a maximum angle of attack limit of
%  12 deg, and minimum and maximum roll angle limits of -/+ 45 deg.
%  to avoid ridiculous negative thrust values.
%
   uminvec = [200;(-3*pi/180);(-45*pi/180)];
   umaxvec = [22000;(12*pi/180);(45*pi/180)];
   unonlinearfeedbackfunct = @(targduma,xargduma) ...
       min([umaxvec,max([uminvec,...
              ulinearfeedbackfunct(targduma,xargduma)],[],2)],[],2);
%
%  This is the nonlinear dynamics model function that will be input
%  to ode45.m.  It includes the effects of the above feedback
%  control law.
% 
   ffunctclosedloopaircraft = @(targdumb,xargdumb) ...
        ffunctaircraft04(xargdumb,...
                         unonlinearfeedbackfunct(targdumb,xargdumb),...
                         m,S,CLalpha,CD0,oneoverpiARe);
%
%  Set up the remaining inputs to ode45.m.
%
   tspan = (0:600)';      
   optionsode45 = odeset('RelTol',1.e-10);
%
%  Run ode45.m.  Afterwards, compute the perturbation from the
%  nominal steady-motion trajectory and the control time history,
%  both absolute and perturbed.
%
   [tclhist,xclhist] = ode45(ffunctclosedloopaircraft,tspan,...
                             x0,optionsode45);
   N = size(tclhist,1);
   xSMhist = ones(N,1)*(xeq') + tclhist*[XdotSM,YdotSM,zeros(1,4)];
   Deltaxclhist = xclhist - xSMhist;
   uclhist = zeros(N,3);
   for k = 1:N
      uclhist(k,:) = unonlinearfeedbackfunct(tclhist(k,1),xclhist(k,:)')';
   end
   Deltauclhist = uclhist - ones(N,1)*(ueq');
   clear k
%
%  Redo the simulation using the linearized closed-loop model.
%  Acl is the closed-loop dynamic matrix such that the closed-
%  loop linearized model takes the form:
%
%    Deltaxdot = Acl*Deltax
%
   Acl = ????;
   ffunctclosedloopaircraftlin = @(targdumb,xargdumb) Acl*xargdumb;
   Deltax0 = x0 - xeq;
   [tclhist_alt,Deltaxclhist_alt] = ...
          ode45(ffunctclosedloopaircraftlin,tspan,Deltax0,optionsode45);
   N = size(tclhist_alt,1);
   xSMhist = ones(N,1)*(xeq') + tclhist_alt*[XdotSM,YdotSM,zeros(1,4)];
   xclhist_alt = xSMhist + Deltaxclhist_alt;
%
%  Redo the nonlinear and linear simulations for open-loop motion,
%  i.e., for motion without feedback.
%
   ffunctopenloopaircraft = @(targdumc,xargdumc) ...
        ffunctaircraft04(xargdumc,ueq,m,S,CLalpha,CD0,oneoverpiARe);
   [tolhist,xolhist] = ode45(ffunctopenloopaircraft,tspan,x0,optionsode45);
   N = size(tolhist,1);
   xSMhist = ones(N,1)*(xeq') + tolhist*[XdotSM,YdotSM,zeros(1,4)];
   Deltaxolhist = xolhist - xSMhist;
   clear N
   ffunctopenloopaircraftlin = @(targdumd,xargdumd) A*xargdumd;
   [tolhist_alt,Deltaxolhist_alt] = ...
          ode45(ffunctopenloopaircraftlin,tspan,Deltax0,optionsode45);
   N = size(tolhist_alt,1);
   xSMhist = ones(N,1)*(xeq') + tolhist_alt*[XdotSM,YdotSM,zeros(1,4)];
   xolhist_alt = xSMhist + Deltaxolhist_alt;
   clear Deltax0 N
%
%  Plot the ground tracks for the open-loop and closed-loop systems.
%
   dumx = [xclhist(:,2);xclhist_alt(:,2);xolhist(:,2);...
           xolhist_alt(:,2);xSMhist(:,2)];
   xmin = min(dumx)*0.001;
   xmin = 10*floor(xmin/10);
   xmax = max(dumx)*0.001;
   xmax = 10*ceil(xmax/10);
   dumy = [xclhist(:,1);xclhist_alt(:,1);xolhist(:,1);...
           xolhist_alt(:,1);xSMhist(:,1)];
   ymin = min(dumy)*0.001;
   ymin = floor(ymin);
   ymax = max(dumy)*0.001;
   ymax = ceil(ymax);
   axisvec = [xmin xmax ymin ymax];
   clear dumx dumy xmin xmax ymin ymax
%   
   figure(1)
   subplot(121)
   hold off
   plot(xclhist(:,2)*0.001,xclhist(:,1)*0.001,'b-','LineWidth',3)
   hold on
   plot(xclhist_alt(:,2)*0.001,xclhist_alt(:,1)*0.001,'r-.','LineWidth',2)
   plot(xSMhist(:,2)*0.001,xSMhist(:,1)*0.001,'g:','LineWidth',1.5)
   hold off
   axis(axisvec)
   grid
   xlabel('Eastward Displacement (km)')
   ylabel('Northward Displacment (km)')
   title('Closed-Loop Ground Track for analandsimaircraft0.mat')
   legend('Nonlinear','Linearized','Target')
   subplot(122)
   hold off
   plot(xolhist(:,2)*0.001,xolhist(:,1)*0.001,'b-','LineWidth',3)
   hold on
   plot(xolhist_alt(:,2)*0.001,xolhist_alt(:,1)*0.001,'r-.','LineWidth',2)
   plot(xSMhist(:,2)*0.001,xSMhist(:,1)*0.001,'g:','LineWidth',1.5)
   hold off
   axis(axisvec)
   grid
   xlabel('Eastward Displacement (km)')
   ylabel('Northward Displacment (km)')
   title('Open-Loop Ground Track for analandsimaircraft0.mat')
   legend('Nonlinear','Linearized','Target')
   clear axisvec xSMhist
%
%  Plot the altitude, airspeed, flight-path angle,
%  and heading angle perturbation time histories from
%  their steady-motion values.  Do this only for the
%  closed-loop system.  Plot both the nonlinear and the
%  linearized responses.
%
   figure(2)
   subplot(411)
   hold off
   plot(tclhist,-Deltaxclhist(:,3),'b-','LineWidth',3)
   hold on
   plot(tclhist_alt,-Deltaxclhist_alt(:,3),'r-.','LineWidth',2)
   hold off
   grid
   ylabel('Altitude Perturbation (m)')
   title(['Closed-loop state perturbation time histories for',...
          ' analandsimaircraft01.mat'])
   legend('Nonlinear','Linearized')
   subplot(412)
   hold off
   plot(tclhist,Deltaxclhist(:,4),'b-','LineWidth',3)
   hold on
   plot(tclhist_alt,Deltaxclhist_alt(:,4),'r-.','LineWidth',2)
   hold off
   grid
   ylabel('Airspeed Perturbation (m/sec)')
   subplot(413)
   hold off
   plot(tclhist,Deltaxclhist(:,5)*(180/pi),'b-','LineWidth',3)
   hold on
   plot(tclhist_alt,Deltaxclhist_alt(:,5)*(180/pi),'r-.','LineWidth',2)
   hold off
   grid
   ylabel('Flight Path Angle Perturbation (deg)')
   subplot(414)
   hold off
   plot(tclhist,Deltaxclhist(:,6)*(180/pi),'b-','LineWidth',3)
   hold on
   plot(tclhist_alt,Deltaxclhist_alt(:,6)*(180/pi),'r-.','LineWidth',2)
   hold off
   grid
   ylabel('Heading Angle Perturbation (deg)')
   xlabel('Time (seconds)')
%
%  Plot the thrust, angle-of-attack, and roll/bank-angle
%  time histories.
%
   figure(3)
   subplot(311)
   hold off
   plot(tclhist,uclhist(:,1),'LineWidth',1.5)
   grid
   ylabel('Thrust (N)')
   title('Control input time histories for analandsimaircraft01.mat')
   subplot(312)
   hold off
   plot(tclhist,uclhist(:,2)*(180/pi),'LineWidth',1.5)
   grid
   ylabel('Angle-of-Attack (deg)')
   subplot(313)
   hold off
   plot(tclhist,uclhist(:,3)*(180/pi),'LineWidth',1.5)
   grid
   ylabel('Roll/Bank-Angle (deg)')
   xlabel('Time (seconds)')
%
%  Save the results.
%
   textcommands = ['These data have been generated by the',...
                   ' commands in script_analandsimaircraft01.m'];
   save analandsimaircraft01 
%
%
%  Displayed outputs when the alternate input file
%  linearizedmodelaircraft01_data04.mat is used to
%  compute the model, analyze the model, design
%  a feedback controller, simulate the closed-loop and
%  open-loop response, and plot results for this
%  file's different steady-motion condition:
%
%  
%  A =
%  
%     1.0e+02 *
%  
%    Columns 1 through 3
%  
%                     0                   0                   0
%                     0                   0                   0
%                     0                   0                   0
%                     0                   0  -0.000000862951937
%                     0                   0   0.000000090753043
%                     0                   0                   0
%  
%    Columns 4 through 6
%  
%     0.000000000000000                   0   1.050000000000000
%    -0.010000000000000                   0   0.000000000000000
%                     0  -1.050000000000000                   0
%    -0.000168152348923  -0.097972100000000                   0
%     0.000017683878724                   0                   0
%                     0                   0                   0
%  
%  B =
%  
%                     0                   0                   0
%                     0                   0                   0
%                     0                   0                   0
%     0.000190111366426  -4.523413162537767                   0
%     0.000000100439162   1.683726275945098                   0
%                     0                   0   0.093306761904762
%  
%  lambdavec =
%  
%   -0.008407617446159 + 0.134935118747691i
%   -0.008407617446159 - 0.134935118747691i
%   -0.000000000000000 + 0.000000000000000i
%    0.000000000000000 + 0.000000000000000i
%    0.000000000000000 + 0.000000000000000i
%    0.000000000000000 + 0.000000000000000i
%  
%  maxreallambda =
%  
%       0
%  
%  This system may be neutrally stable because the maximum
%   eigenvalue real part maximized over all of its eigenvalues
%   appears to be zero to within machine precision.
%   
%  Warning: For eigenvalue lambda = -7.1657e-19, the rank of (lambda*eye(n) - A) is
%   4 but it should be smaller, it should be 2 in order for neutral
%   stability to hold true, because this eigenvalue is repeated 4 times.
%   Therefore, this system is unstable.
%  
%  svsControllabilitymat =
%  
%     1.0e+02 *
%  
%     1.776403628034741
%     0.164446174143162
%     0.097972100000000
%     0.048074657590039
%     0.000933067619048
%     0.000000822563116
%   
%  The system is controllable.
%  
%  K =
%  
%     1.0e+05 *
%  
%    Columns 1 through 3
%  
%    -0.000509667071436   0.000381386641789   0.000038260219039
%     0.000000000079468  -0.000000000100308  -0.000000000417663
%    -0.000000001934614  -0.000000004902788   0.000000002402519
%  
%    Columns 4 through 6
%  
%     0.002614022501539  -0.276005884747732  -1.791666041962433
%     0.000000008278612   0.000002501356656   0.000000272839670
%     0.000000480936388  -0.000000689040680   0.000038458753045