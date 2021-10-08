L=2^16;

M= 40

D = load('hw2_3.mat')
test = D.d;
x = test(:,2);
delta = test(2,1)-test(1,1)
deltaf=1/(2*L*delta);
fs=1/delta;
fmax=fs/2;
t=[0:L-1]'*delta;
Rx=zeros(size(t));
for i=1:L-1, Rx(i)=x(1:L-i+1)'*x(i:L)/(L-i); end

M=50; w=blackman(2*L);
Rxx=toeplitz(w(L+1:2*L-1).*...
Rx(1:L-1),w(L+1:L+M).*Rx(1:M));
Rxs=w(L+2:end).*Rx(2:L);
a_fit=-inv(Rxx'*Rxx)*Rxx'*Rxs;
a_fit=[Rx(1)+Rx(2:M+1)'*a_fit;a_fit];
z=exp(-1i*2*pi*linspace(0,0.5,L+1)');
zz=z;
for j=1:M-1
zz=[z,zz.*repmat(z,1,j)];
end
SxYW=a_fit(1)./(abs((1+zz*a_fit(2:end))).^2);
SxYW=[SxYW(end:-1:2);SxYW(1:end-1)];
SxYW=Rx(1)*SxYW/sum(SxYW);

freq=linspace(-fmax,fmax,length(SxYW));
plot(freq,20.*log(SxYW))
 title('All pole model')
 xlabel('freq in hertz')
 ylabel('magnitude in db')