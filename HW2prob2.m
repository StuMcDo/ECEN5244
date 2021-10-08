Hdatafile = load('hw2_2.mat');
Hdata = Hdatafile.d;
omega = Hdata(:,1);
Hrd = Hdata(:,2);
Hid = Hdata(:,3);
sigma = 0.000000001

H = Hrd+i.*Hid;
H = H +sigma*(randn(size(H))+j*randn(size(H)))
s=i*omega;
%Form design matrix...
L=1; K=2;
A=[]; Acol=ones(size(H));
for nnn=0:L
A=[A,Acol];
Acol=Acol.*s;
end
Acol=H.*s;
for nnn=1:K
A=[A,Acol];
Acol=Acol.*s;
end
%SVD LS solution (xi=0)...
[u,w,v]=svd([real(A);imag(A)],0);
xr=v*inv(w)*u.'*[real(H);imag(H)];
num=xr(L+1:-1:1);
den=[-xr(L+K+1:-1:L+2);1];
zest=roots(num)
pest=roots(den)


 ht = @(s) (num(1,1).*s+num(2,1))./(den(1,1).*s.^2+den(2,1).*s+den(3,1))
% 
 Ht = ht(omega.*i);

Mag = 10*log10(Hrd.^2+Hid.^2);

Magt = 10.*log10(real(Ht).^2+imag(Ht).^2);

Phase = atan(Hid./Hrd);
subplot(2,1,1);

semilogx(omega,Mag) %,omega,Magt
title('Sampled data')
xlabel('frequency rads/s')
ylabel('magnitude')
 subplot(2,1,2);
%   semilogx(omega,Phase)
semilogx(omega,Magt)
title('estimated transfer function')
xlabel('frequency rads/s')
ylabel('magnitude')
   %a = tf(transpose(num),transpose(den));
  % bode(a)
  
kai2 = sum(norm(H-Ht)./sigma)