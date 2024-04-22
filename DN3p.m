function df=DN3p(f,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the function DN3p.m calculates the third partial pseudospectral derivative
% in p
% b length intervall p direction (constant)
% f 2D function (array)
% df third pseuodospectral partial derivative of f in p (array) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M]=size(f);  
fftq = fft2(f); %2D Fourier transform
k = (2*pi/b)^3 * ([0:(M/2-1) 0 (-M/2+1):-1]).^3; 
K = repmat(k,N,1);
dffftq = -1i*K.*fftq; 
df = ifft2(dffftq); %2D inverse Fourier transform
end

