clc;
clear;

%% Some Functions

Weighting_func = @(m_vec,x_vec) (x_vec*m_vec')./sum(m_vec);

AngMom_func = @(m_vec,r_vec,v_vec,rcom) ...
    [ ((r_vec(2,:)-rcom(2)).*v_vec(3,:) - (r_vec(3,:)-rcom(3)).*v_vec(2,:)) *m_vec' ; ...
      ((r_vec(3,:)-rcom(3)).*v_vec(1,:) - (r_vec(1,:)-rcom(1)).*v_vec(3,:)) *m_vec' ; ...
      ((r_vec(1,:)-rcom(1)).*v_vec(2,:) - (r_vec(2,:)-rcom(2)).*v_vec(1,:)) *m_vec' ];

%% Checking with data set at time t0

load mrvdata0_2019

Mtot = sum(mvec);
rcm0 = Weighting_func(mvec,rmat0);
vcm0 = Weighting_func(mvec,vmat0);
ptot0 = vcm0.*Mtot;
h0  = AngMom_func(mvec,rmat0,vmat0,rcm0);

rcm0_c =[ 4.135068659460397; 4.232438883625573; -2.649670988093459];
ptot0_c =[ -3.304955397587092; 2.615681258483552; -12.378512402556995];
vcm0_c =[ -0.448630861116144; 0.355065407616567; -1.680320007510681];
h0_c =[ 0.158538578051034; -0.117644152647727; -0.423922316388524];


format long
disp('Errors are:')
disp('Center of mass')
disp(norm(rcm0-rcm0_c))
disp('Total momentum')
disp(norm(ptot0-ptot0_c))
disp('Velocity of center of mass')
disp(norm(vcm0-vcm0_c))
disp('Angular velocity about center of mass')
disp(norm(h0-h0_c))


%% Solution for data set at time t1

load mrvdata1_2019

rcm1 = Weighting_func(mvec,rmat1);
vcm1 = Weighting_func(mvec,vmat1);
ptot1 = vcm1.*Mtot;
h1  = AngMom_func(mvec,rmat1,vmat1,rcm1);


disp('Solution for Data set at t1:')
disp('Center of mass')
disp(rcm1')
disp('Total momentum')
disp(vcm1')
disp('Velocity of center of mass')
disp(ptot1')
disp('Angular velocity about center of mass')
disp(h1')

%% Forces and Moments

dt = t1 - t0;

Fexttotavg = (ptot1 - ptot0)./dt;

Texttotavg = (h1 - h0)./dt;

disp('Solution for Forces and Moments:')
disp('Total average exteral forces')
disp(Fexttotavg')
disp('Total average exteral moments')
disp(Texttotavg')
