function [val] = LFPSqm(p,q,W0,W1,deltat,m,hquer,gradV,grad3V,a,b,H,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function LFPSqm.m solves the Wigner function equation with a quartic
% symmetric double-well potential with a leap-frog pseudospectral scheme,
% writes data in a txt file, new data appended, at time steps dt*Nstep
% and displays normalization and energy expectation value at console
%
% p,q coordinates (vector)
% W0, W1 initial structures at t=-deltat, t=0 (array)
% deltat time step-size (constant)
% m mass (constant)
% hquer reduced plank constant (constant)
% gradV first derivative of potential (vector)
% grad3V third derivative of potential (vector)
% a length of interval in q direction (constant)
% b length of interval in p direction (constant)
% H hamiltonian function (array)
% T final point of time (constant)
% val structure at time T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
PM = diag(-p/m); 
GRADV= diag(gradV); 
GRAD3V = diag(-grad3V/24); 

tol = 10^-8; 

val_m1 = W0; 
val_0 = W1; 

disp('  time        normalization       energy')
Nend = T/deltat; %number of iterations in time 
Nstep=(Nend-1)/10;
for s = 1:Nstep:Nend 
    t = s*deltat; %time
    %% LFPS scheme
    val = val_m1 + ...
        2*deltat*DNq(val_0,a)*PM + ... % -(p/m) d_q W
        2*deltat*GRADV*DNp(val_0,b) + ...  % +d_q H * d_p W 
        2*deltat*hquer^2*GRAD3V*DN3p(val_0,b);  % - d_q^3 H * d_p^3 W / 24
    val_m1 = val_0; 
    val_0 = val; 
    %% calculation of normalization, energy, saving data, new data appended
        dlmwrite('structure2.txt',val, '-append');
        fprintf(' %d  %e  %e  \n',t,trapz(p,trapz(q,val)), trapz(p,trapz(q,val.*H))); 

end
end

