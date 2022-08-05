clear all;
%close all;

%% Computation of alpha and beta for structural damping
%q4
ome1=1.698848*2*pi;
ome2=3.101310*2*pi;
ome3=3.317977*2*pi;
% %q5
% ome1=1.633141*2*pi;
% ome2=3.312939*2*pi;
% ome3=4.244738*2*pi;
a1=1/(2*ome1);
b1=ome1/2;
a2=1/(2*ome2);
b2=ome2/2;
a3=1/(2*ome3);
b3=ome3/2;
h1=7e-3;
h2=5e-3;
h3=6e-3;

A=[a1 b1;
    a2 b2;
    a3 b3];
b=[h1; h2; h3];
x=(A'*A)^(-1)*A'*b;
alpha=x(1)
beta=x(2)
    
