function W = wigner_harm_osc(p,q,t,hquer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function wigner_harm_osz.m calculates the normalized Wigner function 
% based on a superpostition of the ground and first state of a harmonic 
% oszillator
% 
% p,q coordinates (vector)
% t point of time (constant)
% hquer reduced plank constant (constant)
% W Wigner function (array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(q); 
M = length(p); 

Q = zeros(N,N);
Q(:,1) = q; 
Q(:,end) = ones(N,1); 

P = zeros(N,M); 
P(1,:) = ones(1,M); 
P(end,:) = -p;  

Qterm = zeros(N,N); 
Qterm(:,1)= 2/sqrt(hquer)*cos(t)*q + 1/hquer*q.^2; 
Qterm(:,end)= ones(N,1); 
Pterm = zeros(N,M); 
Pterm(1,:) = ones(1,M); 
Pterm(end,:) = -2/(hquer)*sin(t)*p + 1/hquer*p.^2; 

EXP = exp(-1/hquer*(Q.^2*P.^2)); 

W = 1/(pi*hquer)*(Qterm*Pterm).*EXP; 
 %trapz(p,trapz(q,W)) 
end 