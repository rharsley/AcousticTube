%%!Get the impulse response of an acoustic tube.
%%!@param boreInfo Structure that contains bore information: 
%% boreInfo.r contains samples of cross sectional area that correspond 
%% to positions in boreInfo.x. boreInfo.temp corresponds to temperature for
% simulations
%%@param duration Duration of impulse response in seconds
%%@param FS Sample rate of simulation in Hz
%%@return The impulse response of the acoustic tube
function  impulseResponse = TubeImpulseResponse(boreInfo, duration, FS)

    k=1/FS;
    K = floor(duration/k);

    impulseResponse = zeros(K, 1);
    [rho, c, eta, nu, gamma] = ThermConstants(boreInfo.temp);

    r0=boreInfo.r(1);rN=boreInfo.r(end);
    L=boreInfo.x(end);

    lambda=0.99;h=c*k/lambda;
    N=floor(L/h);h=L/N;lambda=c*k/h;

    xbar = [0:N]*h;
    x=[0.5:N-0.5]*h;
    S0=pi*r0^2;SN=pi*rN^2;
    r=interp1(boreInfo.x, boreInfo.r, x)';
    S=pi*r.^2;
    Sbar=[S0;0.5*(S(2:N)+S(1:N-1));SN];
    rbar=sqrt(Sbar/pi);
    
    %Visctherm
    M=16;%Filter order
    filtername=['MagFilter_0p01'];
    Chat=(gamma-1)/(rho*c^2);
    R0=zeros(N,1);Rq=zeros(N,M);Lq=zeros(N,M);
    G0=zeros(N+1,1);Gq=zeros(N+1,M);Cq=zeros(N+1,M);
    run(filtername)
    for nn=1:N+1
        if nn<N+1
            R0(nn)=(r0/r(nn))^2*exp(a0);
        end
        G0(nn)=(r0/rbar(nn))^2*exp(a0)*(gamma-1)/(rho^2*c^2*nu^2);
        for mm=1:M
            if nn<N+1
                Rq(nn,mm)=(r0/r(nn))^2*exp(a(mm));
                Lq(nn,mm)=exp(a(mm)-b(mm));
            end
            Gq(nn,mm)=(r0/rbar(nn))^2*exp(a(mm))*(gamma-1)/(rho*c^2*nu^2);
            Cq(nn,mm)=exp(a(mm)-b(mm))*(gamma-1)/(rho^2*c^2);
        end
    end

    Gbarq=2*Cq.*Gq./(2*Cq+k*Gq);
    Gbar=G0+sum(Gbarq')';
    E=rho*c^2*k*Chat./(2*Chat+k*Gbar);

    Ap=sparse(1:N+1,1:N+1,(1-E.*Gbar)./(1+E.*Gbar),2*(N+1),2*(N+1))+sparse(N+2:2*(N+1),1:N+1,1,2*(N+1),2*(N+1));
    Bv=-rho*c^2*k./(h*Sbar.*(1+E.*Gbar));
    Bv=sparse(1:N,1:N,Bv(1:N).*S,2*(N+1),2*N)-sparse(2:N+1,1:N,Bv(2:N+1).*S,2*(N+1),2*N);
    Ap0=sparse(1:N+1,1:N+1,2*E.*Gbar./(1+E.*Gbar),2*(N+1),2*(N+1));
    apq=2*E.*Gbarq./(1+E.*Gbar);

    ep=sparse(1:N+1,1:N+1,(2*Chat-k*Gbar)./(2*Chat+k*Gbar),2*(N+1),2*(N+1))+sparse(N+2:2*(N+1),1:N+1,1,2*(N+1),2*(N+1));
    nup=k*Gbar./(2*Chat+k*Gbar);
    nup=sparse(1:N+1,1:N+1,nup,2*(N+1),2*(N+1))+sparse(1:N+1,N+2:2*(N+1),nup,2*(N+1),2*(N+1));
    vp=-2*k*Gbarq./(2*Chat+k*Gbar);

    tp=(2*Cq-k*Gq)./(2*Cq+k*Gq);
    xip=k*Gq./(2*Cq+k*Gq);

    Apq=sparse(2*(N+1),M*(N+1));
    Vp=sparse(2*(N+1),M*(N+1));
    Tp=sparse(M*(N+1),M*(N+1));
    Xp=sparse(M*(N+1),2*(N+1));

    Rbarq=2*Lq.*Rq./(2*Lq+k*Rq);Rbar=R0+sum(Rbarq')';

    Av=sparse(1:N,1:N,(2*rho-k*Rbar)./(2*rho+k*Rbar),2*N,2*N)+sparse(N+1:2*N,1:N,1,2*N,2*N);
    Bp=-2*k./(h*(2*rho+k*Rbar));
    Bp=sparse(1:N,2:N+1,Bp,2*N,2*(N+1))-sparse(1:N,1:N,Bp,2*N,2*(N+1));
    avq=2*k*Rbarq./(2*rho+k*Rbar);

    tv=(2*Lq-k*Rq)./(2*Lq+k*Rq);
    xv=k*Rq./(2*Lq+k*Rq);

    Avq=sparse(2*N,M*N);
    Tv=sparse(M*N,M*N);
    Xv=sparse(M*N,2*N);
    for mm=1:M
        Apq=Apq+sparse(1:N+1,(mm-1)*(N+1)+1:mm*(N+1),apq(:,mm),2*(N+1),M*(N+1));
        Vp=Vp+sparse(1:N+1,(mm-1)*(N+1)+1:mm*(N+1),vp(:,mm),2*(N+1),M*(N+1));
        Tp=Tp+sparse((mm-1)*(N+1)+1:mm*(N+1),(mm-1)*(N+1)+1:mm*(N+1),tp(:,mm),M*(N+1),M*(N+1));
        Xp=Xp+sparse((mm-1)*(N+1)+1:mm*(N+1),1:N+1,xip(:,mm),M*(N+1),2*(N+1))+...
            sparse((mm-1)*(N+1)+1:mm*(N+1),N+2:2*(N+1),xip(:,mm),M*(N+1),2*(N+1));
        
        Avq=Avq+sparse(1:N,(mm-1)*N+1:mm*N,avq(:,mm),2*N,M*N);
        Tv=Tv+sparse((mm-1)*N+1:mm*N,(mm-1)*N+1:mm*N,tv(:,mm),M*N,M*N);
        Xv=Xv+sparse((mm-1)*N+1:mm*N,1:N,xv(:,mm),M*N,2*N)+sparse((mm-1)*N+1:mm*N,N+1:2*N,xv(:,mm),M*N,2*N);
    end
    %Set for impulse response
    Bv(1,1)=2*Bv(1,1);
    Uvec=sparse(1,1,-Bv(1,1)/S(1),2*(N+1),1);
    
    %Radiation parameters
    R1=rho*c;Lr=0.613*rho*rN;R2=0.505*rho*c;Cr=1.111*rN/(rho*c^2);
    G2=1/R2;

    Hr=0.5*k*(1+R1*G2);
    xi=k/(R1*Cr+Hr);
    Er=(rho*c^2*k/h)*(0.5*k/Lr+(0.5*G2+Cr/k)*xi);
    tau=(R1*Cr-Hr)/(R1*Cr+Hr);
    eps=-2*rho*c^2*k/(h*(1+Er));
    nu=eps*((0.5*G2)*(tau+1)+(Cr/k)*(tau-1));
    beta=2*rho*c^2*k/(SN*h*(1+Er));
    alpha=(1-Er)/(1+Er);
    
    
    Uin=zeros(K,1);Uin(1)=1;

    p=zeros(2*(N+1),1);p0=p;pq=zeros(M*(N+1),1);
    v=zeros(2*N,1);vq=zeros(M*N,1);
    vr=0;pr=0;
    
    prog = 0.1 * K;
    prog0 = prog;
    display("Starting loop");
    for nn=1:K
        p=Ap*p+Bv*v+Ap0*p0+Apq*pq+Uvec*Uin(nn);%n -> Uin n-1/2
        p0=ep*p0+nup*p+Vp*pq;%n
        pq=Tp*pq+Xp*(p-p0);%n
        
        p(N+1)=alpha*p(2*(N+1))+beta*S(N)*v(N)+eps*vr+nu*pr;
        vr=vr+0.5*k*(p(N+1)+p(2*(N+1)))/Lr;
        pr=tau*pr+0.5*xi*(p(N+1)+p(2*(N+1)));
        
        v=Av*v+Bp*p+Avq*vq;%n+1/2
        vq=Tv*vq+Xv*v;%n+1/2
       
         impulseResponse(nn) = p(1);
    %Uncomment if want to plot output      
    %     plot([0:N]*h, p(1:N+1)); hold on;
    %      plot(nn*h*[1], 0,'*', 'Markersize',15); hold off; drawnow;
         if(nn > prog) 
           display([num2str(100*nn/K), "%"]);
           prog += prog0;
          endif
    end
    display("Complete")
