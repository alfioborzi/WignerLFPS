%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotqm.m plots data generated in LFPSqm.m in a surface plot
% on the domain (-L,L)x(-K,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pointsq = 64; 
pointsp = 64; 
L = 5; %length of interval (0,L) in q direction
K = 5; %length of interval (O,K) in p direction
a = 2*L; 
b = 2*K; 
deltaq = a/pointsq; %step-size q-dim
deltap = b/pointsp; %step-size p-dim
delta_t = 10^-6; %step-size in time
q = -L + deltaq*(0:pointsq-1); % q coordinate
p = -K + deltap*(0:pointsp-1); % p coordinate

W = dlmread('structure.txt'); %read in data
[r,c]=size(W);

f=0; 
for k = 1:(r/c)
    figure(k)
    surf(p,q,W(f+1:k*c,:)); shading interp; view([0 90]); 
    caxis([-0.3,0.6]); 
    ylabel('$q$','Interpreter','LaTex','FontSize',15); 
    xlabel('$p$', 'Interpreter','LaTex','FontSize',15); 
    f=f+c; 
end
