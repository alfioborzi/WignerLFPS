%% Description main file
% calculate initial structure, potential, hamiltonian function
% apply LFPSerr.m to solve Wigner function equation

%% Constants
m = 1; 
omega = 1;   

%% settings LFPS
pointsq =20; 
pointsp =20; 
L = 4; 
K = 4;
a = 2*L; 
b = 2*K; 
deltaq = a/pointsq; %stepsize q-dim
deltap = b/pointsp; %stepsize p-dim
delta_t = 0.00001; %stepsize in time
q = -L + deltaq*(0:pointsq-1); 
p = -K + deltap*(0:pointsp-1); 
 

%% initial structure
W0=wigner_harm_osc(p,q,-delta_t,1); 
W1=wigner_harm_osc(p,q,0,1);   

%% symmetric double well potential
[V,DV,~,~] = sym_double_well(q,0,-0.5,0); 

%% Hamilton function
H = hamilton(p,m,V);

%% solve Wigner function equation
% fftw('planner','measure'); 
[error,norm,E,val]=LFPSerr(p,q,W0,W1,delta_t,m,DV,deltaq,deltap,a,b,H,0.1);

% dlmwrite('error_harm_01_m6_128_cl.txt',error,'precision',11); 
% dlmwrite('norm_harm_01_m6_128_cl.txt',norm,'precision',11); 
% dlmwrite('en_harm_1_m4_128_cl.txt',E,'precision',11); 
fprintf('%.10d,%.10d, %.10d\n',error(end), norm(end),E(end)); 
