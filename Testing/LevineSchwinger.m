clear all; close all;

%LS Model
% kmax=2*pi*fmax*a/c;
[rho,c]=thermconstants(26.85);
kmax=3.82;Nka=100;a=0.01;
% if kmax>3.82
%     kmax=3.82;
% end
% ka=[0:Nka-1]*kmax/Nka;
ka=linspace(0.001,kmax,Nka);
magR=zeros(1,Nka);
l=zeros(1,Nka);
x0=1e-15;
fun1=@(x,k) atan2(besselj(1,x),(-bessely(1,x)))./(x.*sqrt(k.^2-x.^2));
fun2=@(x,k) log(pi*besselj(1,x).*sqrt((besselj(1,x).^2+bessely(1,x).^2)))./(x.*sqrt(k.^2-x.^2));
fun3=@(x,k) log(1./(2*besseli(1,x,1).*besselk(1,x,1)))./(x.*sqrt(x.^2+k.^2));

for nn=1:Nka
    f1int=integral(@(x)fun1(x,ka(nn)),x0,ka(nn));
    f2int=integral(@(x)fun2(x,ka(nn)),x0,ka(nn));
    f3int=integral(@(x)fun3(x,ka(nn)),x0,inf);
    
    magR(nn)=exp(-2*ka(nn)*f1int/pi);
    l(nn)=(1/pi)*(f2int+f3int);
end
RLS=-magR.*exp(-2i*ka.*l);
ZLS=rho*c*(1+RLS)./(1-RLS);
error('a');

Ns=50;Nk=100;
fmax=kmax*c/a;
sigmamax=0.5*fmax;
sigma=linspace(-sigmamax,sigmamax,Ns);

om=linspace(-2*fmax,fmax,Nk);

magR2=zeros(Ns,Nk);
l2=zeros(Ns,Nk);
s=zeros(Ns,Nk);

fun1=@(x,k) atan(besselj(1,x)./(-bessely(1,x)))./(x.*sqrt(k.^2-x.^2));

for mm=1:Ns
    for nn=1:Nk
        s(mm,nn)=sigma(mm)+1j*om(nn);
        ka=(a/c)*-1j*(sigma(mm)+1j*om(nn));
        wp=[0+1j*imag(ka),real(ka)+1j*imag(ka),real(ka)];
        f1int=integral(@(x)fun1(x,ka),x0,ka);
        f2int=integral(@(x)fun2(x,ka),x0,ka,'RelTol',1e-12,'AbsTol',1e-3);
        f3int=integral(@(x)fun3(x,ka),x0,inf);

        magR2(mm,nn)=exp(-2*ka*f1int/pi);
        l2(mm,nn)=(1/pi)*(f2int+f3int);
     
    end
end
RLS2=-magR2.*exp(-2i*ka.*l2);
ZLS2=rho*c*(1+RLS2)./(1-RLS2);

surf(om,sigma,abs(s.*ZLS2),'Edgecolor','none')
shading interp