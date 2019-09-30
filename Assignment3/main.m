%% Question 2
clear;clc;
close all;
format long

load rotmateuler123_data01 

R01 = rotmateuler123(psi,theta,phi);

disp('Checking Result for Euler angles data set 01')
disp('Error is:')
disp(norm(R01-R_true));

clear
load rotmateuler123_data02

R_true = rotmateuler123(psi,theta,phi);
disp('R_true from Euler angles for data set 02 is:')
disp(R_true);

%% Question 3

clear;

load rotmatquaternion_data05

R01 = rotmatquaternion(q);

disp('Checking Result for quaternion data set 01')
disp('Error is:')
disp(norm(R01-R_true));

clear
load rotmatquaternion_data06

R_true = rotmatquaternion(q);
disp('R_true from quaternion for data set 02 is:')
disp(R_true);
format short
%% Question 4
%%
% 
% $$R_a = R_1(\frac{\pi}{2})$$
%
% $$R_b = R_2(\frac{\pi}{2})$$
%
% $$R_aR_b = R_1(\frac{\pi}{2})R_2(\frac{\pi}{2})$$
%
% $$R_bR_a = R_2(\frac{\pi}{2})R_1(\frac{\pi}{2})$$
%

Ra = rotmateuler123(pi/2,0,0);
disp('Ra is:')
disp(Ra)

Rb = rotmateuler123(0,pi/2,0);
disp('Rb is:')
disp(Rb)

RaRb = Ra*Rb;
RbRa = Rb*Ra;

disp('Ra*Rb is:')
disp(RaRb)
disp('Rb*Ra is:')
disp(RbRa)

disp('These examples of Ra and Rb show that Ra*Rb is not always equal to Rb*Ra')