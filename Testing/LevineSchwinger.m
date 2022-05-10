function [ka,ZLS,Zrlc,RLS,Rrlc,l]=LSRadImp(rho,c,a,Nka,fmax)
%gets radiation impedance from Levine and Schwinger model (ZLS) and lumped
%element model (ZRLC) as a function of wavenumber ka=2pi*f*a/c. Inputs are
%air density rho, speed of sound c and tube radius a and number of
%frequency points Nka. fmax is sets a frequency (Hz) limit on to reduce
%compute time for high resolution. May be adjusted if makes ka>3.82

%LS Model
kmax=2*pi*fmax*a/c;
% if kmax>3.82
%     kmax=3.82;
% end
% ka=[0:Nka-1]*kmax/Nka;
ka=linspace(0,kmax,Nka);
magR=zeros(1,Nka);
l=zeros(1,Nka);

deltax=0.25*(ka(2)-ka(1));

xinf=[0.00000001:0.001:1];
for nn=0:20%powers of 10
    xinf=[xinf,[1:0.001:9]*10^nn];
end
numinf=log(0.5./(besseli(1,xinf,1).*besselk(1,xinf,1)));%besseli and besselk normalised but normalisation constant cancels in this case
numinf(1)=0;
LSprog=Nka/10;LSprog0=LSprog;
for nn=1:Nka
    x=[deltax:deltax:ka(nn)-deltax];
%     num=unwrap(2*atan(-besselj(1,x)./bessely(1,x)))/2;
    num=atan2(besselj(1,x),-bessely(1,x));
    den=(x.*sqrt(ka(nn)^2-x.^2));
    f=num./den;
    if(nn != 1)
         magR(nn)=exp(-(2*ka(nn)/pi)*trapz(x,f));
     else
         magR(1) = 1;
     endif
     
    
    num=log(pi*besselj(1,x).*sqrt(besselj(1,x).^2+bessely(1,x).^2));
    den=x.*sqrt(ka(nn)^2-x.^2);
    f=num./den;
    
    deninf=xinf.*sqrt(xinf.^2+ka(nn)^2);
    finf=numinf./deninf;
    if(nn != 1)
        l(nn)=(1/pi)*(trapz(x,f)+trapz(xinf,finf));
    else
        l(1) = 0.6133;
    endif
    
    if nn>LSprog
        {LSprog/Nka}
        LSprog=LSprog+LSprog0;
    end
end
RLS=-magR.*exp(-2i*ka.*l);
ZLS=rho*c*(1+RLS)./(1-RLS);
%RLC Model
[R1,R2,C,L]=RLCconstants(rho,c,a);
omega=1i*ka*c/a;
Zrlc=(L*(R1+R2)*omega+L*R1*R2*C*(omega).^2)./(R1+R2+(L+R1*R2*C)*omega+L*R2*C*(omega).^2);
Rrlc=(rho*c-Zrlc)./(rho*c+Zrlc);