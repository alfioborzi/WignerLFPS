function [V,DV,D2V,D3V] = sym_double_well(q,a,b,c) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function sym_double_well.m calculates a quartic symmetric double-well
% potential aq^4 - bq^2 + c and its derivatives up to the third one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V =a*q.^4 - b*q.^2+c; 
DV = 4*a*q.^3 - 2*b*q; 
D2V = 12*a*q.^2 - 2*b; 
D3V = 24*a*q; 
end