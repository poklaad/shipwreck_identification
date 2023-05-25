function [y] = spectrplotn(x)
h = 4; 
tau = 4.8*sqrt(h);
A = 0.28*(2*pi).^4*h.^2*tau.^(-4);
B = 0.44*(2*pi).^4*tau.^(-4);
h1 = 3;
tau1 = 4.8*sqrt(h1);
A1 = 0.28*(2*pi).^4*h1.^2*tau1.^(-4);
B1 = 0.44*(2*pi).^4*tau1.^(-4);
k = 5;
n = 4;
p = 9;
m = 8;
y = A.*x.^(-k).*exp(-B.*x.^(-n)) + A1.*x.^(-p).*exp(-B1.*x.^(-m));
end
