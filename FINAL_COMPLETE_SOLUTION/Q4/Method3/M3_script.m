clear all;
close all;

%% Matrix partioning
%q1-q4
load('C:\Users\Gianciccetto\Desktop\INGEGNERD\Magistrale\1 anno\1 semestre\dynamics of mechanical system\esercizi\matlab_labs\Yearwork-20211118\FINAL_COMPLETE_SOLUTION\Q1\initial_structure_mkr.mat');
% q5
% load('C:\Users\Gianciccetto\Desktop\INGEGNERD\Magistrale\1 anno\1 semestre\dynamics of mechanical system\esercizi\matlab_labs\Yearwork-20211118\FINAL_COMPLETE_SOLUTION\Q5\q5_55n_82s_mkr.mat');
%q1-q4
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

%% Impulse (half-sine) definition

fs=100; % sample frequency
dt=1/fs; % time step
T=15; %in the range [6,24] s
t=0:dt:T; 
Fmax=4000;
Ti=2*0.3; % period of the half-sine function
fi=1/Ti; % freq. of the half-sine function
for j=1:length(t)
    if t(j)<=0.3
        imp(j)=Fmax*sin(2*pi*fi*t(j));
    else
        imp(j)=0;
    end
   
end

% % PLOTTING %
% plot(t,imp);grid
% title('Input force: half-sine approx. of impulse')
%% 4.1. DFT %%

%% DFT of input force=approx.impulse
N = length(imp); % number of samples in the signal being analyzed by the DFT

% define the frequency components of the series 
bin_num = 0:N-1;  %freq.expressed in samples
freq = bin_num*fs/N; %freq. expressed in Hz

Fk = fft(imp, N); 
fmax=5;

% PLOTTING %
% Showing just till fmax
figure; subplot(2,1,1); plot(freq,abs(Fk));grid; xlim([0 fmax]); title('DFT Magnitude up to 5Hz')
subplot(2,1,2); plot(freq,angle(Fk)*180/pi);grid; xlim([0 fmax]); title('DFT Phase up to 5Hz')
% % Showing the extended spectrum
figure; subplot(2,1,1); plot(freq,abs(Fk));grid; title('DFT Magnitude')
subplot(2,1,2); plot(freq,angle(Fk)*180/pi);grid; title('DFT Phase')
% % Showing just half of the symmetric spectrum (!! Attention will cut also Fk vector in 2 !!)
% cutOff = ceil(N/2); Fk = Fk(1:cutOff); freq = freq(1:cutOff);
% figure; subplot(2,1,1); plot(freq,abs(Fk));grid; title('DFT Magnitude symmetric half')
% subplot(2,1,2); plot(freq,angle(Fk)*180/pi);grid; title('DFT Phase symmetric half')
% pause

%% 4.2 TIME HISTORY %%
%% Removing the elements of vector Fk related to a frequency span higher than 2*pi*fmax and normalizing the amplitude

for k=1:floor(N/2)+1
        Omk=(k-1)*2*pi/T;
    
    if Omk>2*pi*fmax
        Fk(k)=0;
    else
        Fk(k)=Fk(k)/(N/2+1);
    end
end

%% Creating the Fourier series expressing a fictious periodic input force 
dof_ya=idb(49,2); %ya (A=node 49)
dof_xa=idb(49,1); %xa (A=node 49)

% Initializing of time domain vectors for outputs
xa_t=zeros(1,N);
ya_t=zeros(1,N);
Rx_t=zeros(1,N);
Ry_t=zeros(1,N);
Mp_t=zeros(1,N);
F=zeros(ndof,1);
% Fourier series
for k=1:N/2+1
    ome=(k-1)*2*pi/T;
    Om0=2*pi/T;
    F(dof_ya)=abs(Fk(k))*exp(i*angle(Fk(k)));
    A=-ome^2*MFF+i*ome*CFF+KFF;
    x=A\F;
    A_cf=-ome^2*MCF+i*ome*CCF+KCF;
    Rp=A_cf*x;
    % Output in freq. domain
    xa=x(dof_xa);
    ya=x(dof_ya);
    Rx=Rp(1);
    Ry=Rp(2);
    Mp=Rp(3);
    % Output in time domain
    xa_t=xa_t+abs(xa)*cos(ome*t+angle(xa));
    ya_t=ya_t+abs(ya)*cos(ome*t+angle(ya));
    Rx_t=Rx_t+abs(Rx)*cos(ome*t+angle(Rx));
    Ry_t=Ry_t+abs(Ry)*cos(ome*t+angle(Ry));
    Mp_t=Mp_t+abs(Mp)*cos(ome*t+angle(Mp));

end       


% PLOTTING %
%ya
figure; plot(t,ya_t); grid; title('Vertical displacement of point A')
%xa
figure; plot(t,xa_t); grid; title('Horizontal displacement of point A')
%Rx
figure; plot(t,Rx_t); grid; title('Horizontal constraint force at point P')
%Ry
figure; plot(t,Ry_t); grid; title('Vertical constraint force at point P')
%Mp
figure; plot(t,Mp_t); grid; title('Constraint moment at point P')