clear all;
% close all;

%% Matrix partioning
%q1-q4
load('C:\Users\Gianciccetto\Desktop\INGEGNERD\Magistrale\1 anno\1 semestre\dynamics of mechanical system\esercizi\matlab_labs\Yearwork-20211118\FINAL_COMPLETE_SOLUTION\Q1\initial_structure_mkr.mat');
%q5
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
% figure; subplot(2,1,1); plot(freq,abs(Fk));grid; title('DFT Magnitude')
% subplot(2,1,2); plot(freq,angle(Fk)*180/pi);grid; title('DFT Phase')
% % Showing just half of the symmetric spectrum (!! Attention will cut also Fk vector in 2 !!)
% cutOff = ceil(N/2); Fk = Fk(1:cutOff); freq = freq(1:cutOff);
% figure; subplot(2,1,1); plot(freq,abs(Fk));grid; title('DFT Magnitude symmetric half')
% subplot(2,1,2); plot(freq,angle(Fk)*180/pi);grid; title('DFT Phase symmetric half')
% pause

%% 4.2 TIME HISTORY %%
%% Removing the elements of vector Fk related to a frequency span higher than 2*pi*fmax

for k=2:N
    if k<N/2+1
        Omk=(k-1)*2*pi/T;
    else 
        Omk=(N-k+1)*2*pi/T;
    end
    if Omk>2*pi*fmax
        Fk(k)=0;
    end
end



%% Creating the Fourier coefficient Fkt an phases to use in the third form of the Fourier series
Fkt=zeros(size(Fk));
phik=zeros(size(Fk));
for k=1:N/2+1
    Fkt(k)=2*abs(Fk(k))/N;
    phik(k)=angle(Fk(k)/N);
    if k==1
        Fkt(k)=Fk(k)/N;
        phik(k)=0;
    end
end
%% Creating the Fourier series expressing a fictious periodic input force 
dof_ya=idb(49,2); %ya (A=node 49)
dof_xa=idb(49,1); %xa (A=node 49)

% Initializing of time domain vectors for outputs
xa=zeros(1,N);
ya=zeros(1,N);
Rx=zeros(1,N);
Ry=zeros(1,N);
Mp=zeros(1,N);

% Fourier series
for k=1:N/2+1
    ome=(k-1)*2*pi/T;
    Om0=2*pi/T;
%     ck=Fk(k)/N;
    A=-ome^2*MFF+i*ome*CFF+KFF;
    A_cf=-ome^2*MCF+i*ome*CCF+KCF;
    G_a=A^(-1); % FRF matrix for displacements
    G_r=A_cf*A^(-1); % FRF matrix for displacements
    
    % computing the different components of the sum based on the value of intex k
    if k==1 %x0 computation
        %ome=0 %A=KFF %A_cf=KCF 
        G_a=KFF^(-1);
        G_r=KCF*KFF^(-1);
    end 
    % Fourier coeff. for each output
    xk_xa=abs(G_a(dof_xa,dof_ya))*Fkt(k);
    xk_ya=abs(G_a(dof_ya,dof_ya))*Fkt(k);
    xk_Rx=abs(G_r(1,dof_ya))*Fkt(k); %looking at idb the #dof of clamp P are in the end but still before than the ones of clamp Q
    xk_Ry=abs(G_r(2,dof_ya))*Fkt(k); 
    xk_Mp=abs(G_r(3,dof_ya))*Fkt(k); 
    % Definig the phase contributions
    psik_xa=angle(G_a(dof_xa,dof_ya))+phik(k);
    psik_ya=angle(G_a(dof_ya,dof_ya))+phik(k);
    psik_Rx=angle(G_r(1,dof_ya))+phik(k);
    psik_Ry=angle(G_r(2,dof_ya))+phik(k);
    psik_Mp=angle(G_r(3,dof_ya))+phik(k);
    
    if 1<k<=N/2+1
%        % in time domain outputs
        xa=xa+xk_xa*cos((k-1)*Om0*t+psik_xa); 
        ya=ya+xk_ya*cos((k-1)*Om0*t+psik_ya);
        Rx=Rx+xk_Rx*cos((k-1)*Om0*t+psik_Rx);
        Ry=Ry+xk_Ry*cos((k-1)*Om0*t+psik_Ry);
        Mp=Mp+xk_Mp*cos((k-1)*Om0*t+psik_Mp);
    end

end       


% PLOTTING %
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