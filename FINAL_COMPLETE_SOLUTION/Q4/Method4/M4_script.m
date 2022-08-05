clear all;
close all;
%% Matrix partioning
%q1-q4
load('C:\Users\Gianciccetto\Desktop\INGEGNERD\Magistrale\1 anno\1 semestre\dynamics of mechanical system\esercizi\matlab_labs\Yearwork-20211118\FINAL_COMPLETE_SOLUTION\Q1\initial_structure_mkr.mat');
%q5
% load('C:\Users\Gianciccetto\Desktop\INGEGNERD\Magistrale\1 anno\1 semestre\dynamics of mechanical system\esercizi\matlab_labs\Yearwork-20211118\FINAL_COMPLETE_SOLUTION\Q5\q5_55n_82s_mkr.mat');
% q1-q4
ntot=3*53; %159
ndof=153;
% %q5
% ntot=3*55;
% ndof=159;
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);
CFF=R(1:ndof,1:ndof);
% MFC=M(1:ndof,ndof+1:ntot);
% CFC=R(1:ndof,ndof+1:ntot);
% KFC=K(1:ndof,ndof+1:ntot);
MCF=M(ndof+1:ntot,1:ndof);
CCF=R(ndof+1:ntot,1:ndof);
KCF=K(ndof+1:ntot,1:ndof);
% MCC=M(ndof+1:ntot,ndof+1:ntot);
% CCC=R(ndof+1:ntot,ndof+1:ntot);
% KCC=K(ndof+1:ntot,ndof+1:ntot);

dof_ya=idb(49,2); %ya (A=node 49)
dof_xa=idb(49,1); %xa (A=node 49)

extract_ya=zeros(ndof,1); extract_ya(dof_ya)=1;%vector that defines the position at which the approx. impulse force is applied

sim('q4_ode');
% Output definition
ya=ans.x(:,dof_ya);
xa=ans.x(:,dof_xa);
Rx=ans.Rp(:,1);
Ry=ans.Rp(:,2);
Mp=ans.Rp(:,3);
% PLOTTING %
t=ans.tout;
%ya
figure; plot(t,ya); grid; title('Vertical displacement of point A')
%xa
figure; plot(t,xa); grid; title('Horizontal displacement of point A')
%Rx
figure; plot(t,Rx); grid; title('Horizontal constraint force at point P')
%Ry
figure; plot(t,Ry); grid; title('Vertical constraint force at point P')
%Mp
figure; plot(t,Mp); grid; title('Constraint moment at point P')