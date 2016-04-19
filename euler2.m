function euler2
% Discretization details
xmin=-50;
xmax=50;
tmax=25;
% number of cells control
CF=1;
X=CF*100+1;
% Constant courant
courant=2;
dx=(xmax-xmin)/(X-1);
dt=dx/courant;
N=(tmax/dt)+1;
dt=(tmax)/(N-1);
N=round(N);
% gamma
g=1.4;
% initial primitive variables w=[density velocity pressure]
w=zeros(3,N,X);
% Initial Conditions
for j=1:round(X/2)
    w(:,1,j)=[1 0 2];
end
for j=round(X/2):X
    w(:,1,j)=[1 0 1];
end

for n=2:N
    %     boundary
    w(:,n-1,1)=w(:,n-1,2);
    for i=2:(X-1)
        %         downwind flux
        f2=flux(w(:,n-1,i),w(:,n-1,i+1),dt,dx,g);
        %         upwind flux
        f1=flux(w(:,n-1,i-1),w(:,n-1,i),dt,dx,g);
        
        u1=w2u(w(:,n-1,i),g);
        u2=u1-(dt/dx)*(f2-f1);
        
        %         update primitive variables
        w(:,n,i)=u2w(u2,g);
        
    end
    %     boundary
    w(:,n,X)=w(:,n,X-1);
end
exact=linspace(xmin,xmax,1000);
for j=1:1000
    wexact(:,j)=rmannsol(w(:,1,1),w(:,1,X),g,(xmin+j*((xmax-xmin)/1000))/(tmax));
end

subplot(1,3,1)

plot(xmin:dx:xmax,squeeze(w(1,n,:))','o',exact,(wexact(1,:))')
axis([xmin xmax 0 2])
title(strcat('2nd order Richtmeyer, dt/dx=',num2str(dt/dx)))
ylabel('Density')
xlabel('x')
subplot(1,3,2)

plot(xmin:dx:xmax,squeeze(w(2,n,:))','o',exact,(wexact(2,:))')
axis([xmin xmax 0 2])
ylabel('Velocity')
xlabel('x')
subplot(1,3,3)

plot(xmin:dx:xmax,squeeze(w(3,n,:))','o',exact,(wexact(3,:))')
axis([xmin xmax 0 2])
ylabel('Pressure')
xlabel('x')
legend('numerical solution','exact solution')

    function [uout]=w2u(win,g)
        uout=zeros(1,3);
        uout(1)=win(1);
        uout(2)=win(2)*win(1);
        uout(3)=win(3)/(g-1)+0.5*win(1)*win(2)^2;
    end
    function [wout]=u2w(uin,g)
        wout(1)=uin(1);
        wout(2)=uin(2)/uin(1);
        wout(3)=(g-1)*(uin(3)-0.5*wout(1)*wout(2)^2);
    end
    function [fout]=w2f(win,g)
        fout(1)=win(1)*win(2);
        fout(2)=win(1)*win(2)^2+win(3);
        fout(3)=((win(3)/(g-1)+0.5*win(1)*win(2)^2)+win(3))*win(2);
    end
    function [fout]=u2f(uin,g)
        fout(1)=uin(2);
        fout(2)=0.5*(3-g)*(uin(2)^2)/uin(1)+(g-1)*uin(3);
        fout(3)=0.5*(1-g)*(uin(2)^3)/(uin(1)^2)+g*uin(3)*uin(2)/uin(1);
    end
% calculate flux using Richtmeyer method
    function [fluxo]=flux(w1,w2,dt,dx,g)
        ul=w2u(w1,g);
        ur=w2u(w2,g);
        fl=w2f(w1,g);
        fr=w2f(w2,g);
        fluxo=0.5*(ur+ul)-0.5*(dt/dx)*(fr-fl);
        fluxo=u2f(fluxo,g);
    end
end