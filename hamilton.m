function H = hamilton(p,m,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the function hamilton.m calculates the Hamiltonian function H = p^2/2m+V(q)
% m mass (constant)
% v potential (vector)
% H Hamiltonian function (array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(v); 
M=length(p); 

V = zeros(N,N); 
V(:,1)=v; 
V(:,end)=ones(N,1); 
 
T =zeros(N,M);  
T(1,:)=ones(1,M); 
T(end,:)=p.^2/(2*m); 

% the H(p,q) = p^2/2m+V(q) function on the phase space
H=V*T; 

end