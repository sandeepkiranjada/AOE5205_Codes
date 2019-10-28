%% Problem 3
clc; clear;
format long

load moicalcs01_data

[my_Mtot,my_rcmtot,my_IMoItot] = momentofinertia01(mvec,rcmmat,IMoIarray);

Mtot = 2.058306536059932e+02;

rcmtot = [0.409101895949622;0.526651850364819;0.058154823743388];

IMoItot = 1.0e+04 *[...
1.047260719697208 0.028550325166975 0.040407207761532;...
0.028550325166975 1.063729452526229 0.010880405279456;...
0.040407207761532 0.010880405279456 1.146269471862332];

disp('Checking with data set moicalcs01_data')
disp('Error in Mtot')
disp(Mtot-my_Mtot);
disp('Error in rcmtot')
disp(rcmtot-my_rcmtot);
disp('Error in IMoItot')
disp(IMoItot-my_IMoItot);


clear;

load moicalcs02_data
disp('Results with data set moicalcs02_data')
[Mtot,rcmtot,IMoItot] = momentofinertia01(mvec,rcmmat,IMoIarray)

%% Problem 4
clear;
a = 0.4;
b = 0.2;
c = 0.8;

l = 2.1;
w = 0.6;

M = 20;
m = 0.6;
theta = 0.523599;

rcm_panel_left = [0;-b/2-l/2;c/2];
rcm_box = [0;0;0];
rcm_panel_right = [0;b/2+l/2;c/2];

rcmmat = [rcm_panel_left,rcm_box,rcm_panel_right];
mvec = [m,M,m];

I_polar_panel = (m/12)*(l^2+w^2);
I_alongl_panel = (m/12)*(w^2);
I_alongw_panel = (m/12)*(l^2);

R2theta_pr = [cos(-theta) 0 -sin(-theta);...
              0           1            0;...
              sin(-theta) 0  cos(-theta)];

Ipanel_left_principle = diag([I_alongw_panel I_alongl_panel I_polar_panel]);
Ipanel_right_principle = diag([I_alongw_panel I_alongl_panel I_polar_panel]);

Ipanel_left_b = R2theta_pr*Ipanel_left_principle*R2theta_pr';
Ipanel_right_b = R2theta_pr*Ipanel_right_principle*R2theta_pr';

Ibox_i_b = (M/12)*(b^2+c^2);
Ibox_j_b = (M/12)*(c^2+a^2);
Ibox_k_b = (M/12)*(a^2+b^2);

Ibox_b = diag([Ibox_i_b Ibox_j_b Ibox_k_b]);

IMoIarray(:,:,1) = Ipanel_left_b;
IMoIarray(:,:,2) = Ibox_b;
IMoIarray(:,:,3) = Ipanel_right_b;

[my_Mtot,my_rcmtot,my_IMoItot] = momentofinertia01(mvec,rcmmat,IMoIarray);



Mtot = 21.200000000000003;
rcmtot = [0 ; 0;0.022641509433962];

IMoItot =[...
3.351465415801186 0 0.015588461307349;...
0 1.550465408805032 0;...
0.015588461307349 0 2.388333326337180];

disp('Checking with test data')
disp('Error in Mtot')
disp(Mtot-my_Mtot);
disp('Error in rcmtot')
disp(rcmtot-my_rcmtot);
disp('Error in IMoItot')
disp(IMoItot-my_IMoItot);

clear;
a = 0.3;
b = 0.4;
c = 0.6;

l = 1.1;
w = 0.5;

M = 15;
m = 0.8;
theta = 0.34906585;

rcm_panel_left = [0;-b/2-l/2;c/2];
rcm_box = [0;0;0];
rcm_panel_right = [0;b/2+l/2;c/2];

rcmmat = [rcm_panel_left,rcm_box,rcm_panel_right];
mvec = [m,M,m];

I_polar_panel = (m/12)*(l^2+w^2);
I_alongl_panel = (m/12)*(w^2);
I_alongw_panel = (m/12)*(l^2);

R2theta_pr = [cos(-theta) 0 -sin(-theta);...
              0           1            0;...
              sin(-theta) 0  cos(-theta)];

Ipanel_left_principle = diag([I_alongw_panel I_alongl_panel I_polar_panel]);
Ipanel_right_principle = diag([I_alongw_panel I_alongl_panel I_polar_panel]);

Ipanel_left_b = R2theta_pr*Ipanel_left_principle*R2theta_pr';
Ipanel_right_b = R2theta_pr*Ipanel_right_principle*R2theta_pr';

Ibox_i_b = (M/12)*(b^2+c^2);
Ibox_j_b = (M/12)*(c^2+a^2);
Ibox_k_b = (M/12)*(a^2+b^2);

Ibox_b = diag([Ibox_i_b Ibox_j_b Ibox_k_b]);

IMoIarray(:,:,1) = Ipanel_left_b;
IMoIarray(:,:,2) = Ibox_b;
IMoIarray(:,:,3) = Ipanel_right_b;

[my_Mtot,my_rcmtot,my_IMoItot] = momentofinertia01(mvec,rcmmat,IMoIarray);



Mtot = 21.200000000000003;
rcmtot = [0 ; 0;0.022641509433962];

IMoItot =[...
3.351465415801186 0 0.015588461307349;...
0 1.550465408805032 0;...
0.015588461307349 0 2.388333326337180];

disp('Results:')
disp('Mtot')
disp(my_Mtot);
disp('rcmtot')
disp(my_rcmtot);
disp('IMoItot')
disp(my_IMoItot);