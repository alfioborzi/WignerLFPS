function df=DNq(f,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the function DNq.m calculates the first partial pseudospectral derivative
% in q
% a length intervall q direction (constant)
% f 2D function (array)
% df first pseuodospectral partial derivative of f in q (array) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M]=size(f); 
fftq = fft2(f); %2D Fourier transform 
k = 2*pi/a * [0:(N/2-1) 0 (-N/2+1):-1].'; 
K = repmat(k,1,M);
dffftq = 1i*K.*fftq; %2D inverse Fourier transform
df = ifft2(dffftq); 
end