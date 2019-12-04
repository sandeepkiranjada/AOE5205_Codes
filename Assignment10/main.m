clc; clear; close all;

format long

% load steadystateaircraft01_data01
% 
% [gammaeq,Teq,alphaeq,phieq,XdotSM,YdotSM,iflagterm,niter] = ...
%              solvesteadystateaircraft01(Zeq,Veq,psieq,m,S,CLalpha,...
%                                         CD0,oneoverpiARe,...
%                                         nitermax);
%                                     
%                                     
% gammaeq_chk = 0;
% Teq_chk = 3.240358191011175e+03;
% alphaeq_chk = 0.135886048735789;
% phieq_chk = 0;
% XdotSM_chk = -22.231309316168467;
% YdotSM_chk = 61.080020351084045;
% 
% disp(gammaeq_chk-gammaeq);
% disp(Teq_chk-Teq)
% disp(alphaeq_chk-alphaeq)
% disp(phieq_chk-phieq)
% disp(XdotSM_chk-XdotSM)
% disp(YdotSM_chk-YdotSM)

clear

load steadystateaircraft01_data02

[gammaeq,Teq,alphaeq,phieq,XdotSM,YdotSM,iflagterm,niter] = ...
             solvesteadystateaircraft01(Zeq,Veq,psieq,m,S,CLalpha,...
                                        CD0,oneoverpiARe,...
                                        nitermax);
                                    
                                    
disp('gammaeq');
disp(gammaeq);
disp('Teq')
disp(Teq)
disp('alphaeq')
disp(alphaeq)
disp('phieq')
disp(phieq)
disp('XdotSM')
disp(XdotSM)
disp('YdotSM')
disp(YdotSM)