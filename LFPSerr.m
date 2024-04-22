function [error,normalize,E,val] = LFPSerr(p,q,W0,W1,deltat,m,gradV,deltaq,deltap,a,b,H,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function LFPSqm.m solves the Wigner function equation with a harmonic 
% oscillator with a leap-frog pseudospectral scheme, calculates error
% between exact solution and numerical solution in L2 norm, normalization 
% and energy expectation value at every time-step. 
%
% p,q coordinates (vector)
% W0, W1 initial structures at t=-deltat, t=0 (array)
% deltat time step-size (constant)
% m mass (constant)
% gradV first derivative of potential (vector)
% deltaq step-size q direction
% deltap step-size p direction
% a length of interval in q direction (constant)
% b length of interval in p direction (constant)
% H hamiltonian function (array)
% T final point of time (constant)
% error error values (vector)
% normalize normalization values (vector)
% E energy expectation values (vector)
% val structure at time T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
PM = diag(-p/m); 
GRADV= diag(gradV);  

val_m1 = W0; 
val_0 = W1; 
Tend = T/deltat;
error = zeros(1,Tend); 
normalize = error;  
E = error; 
for s = 1:Tend 
    t = s*deltat; %time
    %% LFPS 
    val = val_m1 + ...
        2*deltat*DNq(val_0,a)*PM + ...
        2*deltat*GRADV*DNp(val_0,b);
    %% exact solution
    valex = wigner_harm_osc(p,q,t,1);
    %% error in L2 norm
    error(s)=sqrt(deltap*deltaq)*norm((valex-val),'fro'); %L2norm(valex,val,t,deltaq,deltap,q,p);
    %% normalization, energy expectation value
    normalize(s)= trapz(p,trapz(q,val)); 
    E(s)= trapz(p,trapz(q,val.*H));
    val_m1 = val_0; 
    val_0 = val; 
end 

% E = trapz(p,trapz(q,val.*H)); 

end

