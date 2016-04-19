function muscl
% Discretization details
xmin=-50;
xmax=50;
tmax=25;
% number of cells control, choose number + ghost cells + 1 to start indexes
% from 1
X=100+4+1;
% estimate maximum number of time steps
N=2*X;
% Constant dx
dx=(xmax-xmin)/(X-4-1);
t=0;
n=2;
% set maximum number of loop iterations before program timeout
timeout=1e4;
% gamma
g=1.4;
% initial primitive variables w=[density velocity pressure]
w=zeros(3,N,X);
% Initial Conditions
for j=1:(X-1)/2
    w(:,1,j)=[1 0 1];
end
for j=round(X/2):X
    w(:,1,j)=[1 0 2];
end

while(t<tmax)
    %     boundary
    w(:,n-1,1)=w(:,n-1,3);
    w(:,n-1,2)=w(:,n-1,4);
    dt=0.5;
    %     function to set constant courant number (CFL=1) not required for this
    %     task
    %     dt=setdt(squeeze(w(:,n-1,:)),dx,1,g);
    
    %     check last time step
    if (t+dt)>tmax
        dt=tmax-t;
        t=t+dt;
    else
        t=t+dt;
    end
    
    for i=3:(X-2)
        %         downwind prediction
        [wpred2,dw2]=pred(w(:,n-1,i),w(:,n-1,i+1),w(:,n-1,i+2),dt,dx,g);
        
        %
        [wpred1,dw1]=pred(w(:,n-1,i-1),w(:,n-1,i),w(:,n-1,i+1),dt,dx,g);
        
        %         upwind prediction
        [wpred0,dw0]=pred(w(:,n-1,i-2),w(:,n-1,i-1),w(:,n-1,i),dt,dx,g);
        
        w2L=0.5*(w(:,n-1,i)+wpred1')+0.5*dw1;
        w2R=0.5*(w(:,n-1,i+1)+wpred2')-0.5*dw2;
        f2=w2f(rmannsol(w2L,w2R,g,0),g);
        
        
        w0L=0.5*(w(:,n-1,i-1)+wpred0')+0.5*dw0;
        w0R=0.5*(w(:,n-1,i)+wpred1')-0.5*dw1;
        f0=w2f(rmannsol(w0L,w0R,g,0),g);
        
        
        u1=w2u(w(:,n-1,i),g);
        u2=u1-(dt/dx)*(f2-f0);
        %         update primitive variables
        w(:,n,i)=u2w(u2,g);
        
    end
    %     boundary
    w(:,n,X)=w(:,n,X-2);
    w(:,n,X-1)=w(:,n,X-3);
    
    %     live plotting
    %%%%%%%%%%%%%%%%%%%%
    %         subplot(1,3,1)
    %
    %         plot(xmin:dx:xmax,squeeze(w(1,n,3:X-2))','o')
    %         axis([xmin xmax 0 2])
    %         title(strcat('2nd order Richtmeyer, dt/dx=',num2str(dt/dx)))
    %         ylabel('Density')
    %         xlabel('x')
    %         subplot(1,3,2)
    %
    %         plot(xmin:dx:xmax,squeeze(w(2,n,3:X-2))','o')
    %         axis([xmin xmax -2 2])
    %         ylabel('Velocity')
    %         xlabel('x')
    %         subplot(1,3,3)
    %
    %         plot(xmin:dx:xmax,squeeze(w(3,n,3:X-2))','o')
    %         axis([xmin xmax 0 2])
    %         ylabel('Pressure')
    %         xlabel('x')
    %         legend('numerical solution','exact solution')
    %         pause(0.02)
    %%%%%%%%%%%%%%%%%%%%%
    
    %     advance timestep
    n=n+1;
    %     check for timeout
    if n>timeout
        t
        error('timeout')
    end
end
% calculate exact solution
exact=linspace(xmin,xmax,1000);
for j=1:1000
    wexact(:,j)=rmannsol(w(:,1,3),w(:,1,X-2),g,(xmin+j*((xmax-xmin)/1000))/(tmax));
end
% plot results
subplot(1,3,1)
hold all
plot(xmin:dx:xmax,squeeze(w(1,n-1,3:X-2))','o',exact,(wexact(1,:))')

axis([xmin xmax 0 2])
title(strcat('MUSCL scheme, constant dt/dx=',num2str(dt/dx),',intermediate (beta=1.5)'))
ylabel('Density')
xlabel('x')

subplot(1,3,2)
hold all
plot(xmin:dx:xmax,squeeze(w(2,n-1,3:X-2))','o',exact,(wexact(2,:))')
axis([xmin xmax -2 2])
ylabel('Velocity')
xlabel('x')

subplot(1,3,3)
hold all
plot(xmin:dx:xmax,squeeze(w(3,n-1,3:X-2))','o',exact,(wexact(3,:))')
axis([xmin xmax 0 2])
ylabel('Pressure')
xlabel('x')
legend('numerical solution','exact solution')

% conversions between primitive conserved and flux variables
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

% calculate predictor step of Muscl-Hancock
    function [wp,dw]=pred(w0,w1,w2,dt,dx,g)
        %         choose slope limiter beta is superbee for b=2, minmod for
        %         b=1,1<b<2 for intermediate dissipation
        dw=beta(w1-w0,w2-w1,1.5);
        % dw=0;
        %         dw=0.5*((w2'-w0'));
        w2l=w1+0.5*dw';
        w0r=w1-0.5*dw';
        
        up=w2u(w1,g)-(dt/dx)*(w2f(w2l,g)-w2f(w0r,g));
        dw=dw';
        wp=u2w(up,g);
        
    end
% beta slope limiter
    function [lim]=beta(w1,w2,b)
        lim=zeros(1,3);
        for q=1:3
            if w2(q)>=0
                lim(q)=max([0 min([b*w1(q),w2(q)]) min([w1(q),b*w2(q)])]);
            else
                lim(q)=min([0 max([b*w1(q),w2(q)]) max([w1(q),b*w2(q)])]);
            end
        end
    end
% setdt function for constant CFL
    function [dt]=setdt(win,dx,CFL,gamma)
        dt=1000;
        for count=1:size(win,2)
            c=(gamma*win(3,count)/win(1,count))^0.5;
            speed=abs(win(2,count))+c;
            if speed>(CFL*dx/dt)
                dt=CFL*dx/speed;
            end
        end
    end
end