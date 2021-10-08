
 L=2^16;

 N=512;


D = load('hw2_3.mat')
test = D.d;
x = test(:,2);
D = load('hw2_3.mat')
test = D.d;
x = test(:,2);

delta = test(2,1)-test(1,1)
 t=[0:L-1]'*delta;
fs=1/delta;
fmax=fs/2;
%Autocorrelation estimate
Rx=zeros(size(t));
for n=1:L-1
Rx(n)=x(1:L-n+1)'*x(n:L)/(L-n);
end
%Correlation coefficient estimate
rhox=Rx/(std(x)^2);
%Spectral estimates
Sx=abs(fft([0;Rx(N:-1:2);Rx(1:N-1);0]));
Sx=Sx*Rx(1)/sum(Sx);
Sxb=abs(fft(blackman(2*N).*...
[0;Rx(N:-1:2);Rx(1:N-1);0]));
Sxb=Sxb*Rx(1)/sum(Sxb);
Sxb= fftshift(Sxb);
freq=linspace(-fmax,fmax,length(Sxb));
plot(freq,10.*log(Sxb))
 title('Blackman Tukey With Blackman window N = 512')
 xlabel('freq in hertz')
 ylabel('magnitude in db')