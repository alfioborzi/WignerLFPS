
% Stability and Accuracy of a Pseudospectral Scheme for the Wigner Function Equation
% Andrea Thomann, Alfio Borzì
% Numer. Methods Partial Differential Equations, 33 (2017), pp. 62–87.
%
%% Description main file
% calculate initial structure, potential, hamiltonian function
% apply LFPSqm.m to solve the Wigner function equation

%% Constants
m = 1; 
omega = 1;
hquer = 0.5;  

%% settings LFPS
pointsq = 64; 
pointsp = 64; 
L = 5; 
K = 5;
a = 2*L; 
b = 2*K; 
deltaq = a/pointsq; %stepsize q-dim
deltap = b/pointsp; %stepsize p-dim
delta_t = 10^-5; %stepsize in time
q = -L + deltaq*(0:pointsq-1); %q coordinate
p = -K + deltap*(0:pointsp-1); %p coordinate

%% initial structure 
W0=wigner_harm_osc(p,q,-delta_t,hquer); 
W1=wigner_harm_osc(p,q,0,hquer); 

%% symmetric double well potential
[V,DV,~,D3V] = sym_double_well(q,0.1,1,2.5); 
% plot(q,V); 
% ylabel('$V(q)$','Interpreter','LaTex','FontSize',15); 
% xlabel('$q$', 'Interpreter','LaTex','FontSize',15);
%% Hamiltonian function
H = hamilton(p,m,V);
%% solve Wigner function equation
val = LFPSqm(p,q,W0,W1,delta_t,m,hquer,DV,D3V,a,b,H,2);  
%