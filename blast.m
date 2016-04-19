function blast
% Discretization details
xmin=0;
% boundary between smaller dx and larger dx at xmed, also the boundary
% between the initial conditions
xmed=0.1;
xmax=1.5;
tmax=0.5;
% number of cells control, choose number + ghost cells + 1 since indexes
% start from 1
X=100+2;
XHR=20+2+1;
% check to make sure smaller dt on left side
if ((XHR-3)/(X-2))<((xmed-xmin)/(xmax-xmed))
    error('XHR needs to have higher resolution, increase # of cells in HR region')
end
% estimate number of time steps
N=2*(X+XHR);
% Constant dx
dx2=(xmax-xmed)/(X-2);
dx1=(xmed-xmin)/(XHR-2-1);
% initialize time variables
t=0;
t1=0;
n=2;
timeout=1e4;
% gamma
g=1.4;
% initial primitive variables w=[density velocity pressure xposition]
% xposition is not a primitive variable but it is convenient to store it in
% the same array
w=zeros(4,N,X+XHR);
tn=zeros(1,N);
% Initial Conditions
for j=1:XHR
    w(1:3,1,j)=[100 0 1000];
    w(4,1,j)=xmin+(j-3)*dx1;
end
for j=(XHR+1):(X+XHR)
    w(1:3,1,j)=[1 0 1];
    w(4,1,j)=((j-XHR)*dx2)+xmed;
end

while(t<tmax)
    %     boundary (solid wall since left bounday is symmetry surface)
    w(1:3,n-1,1)=w(1:3,n-1,4);
    w(1:3,n-1,2)=w(1:3,n-1,3);
    w(2,n-1,1:2)=-w(2,n-1,1:2);
    dt1=setdt(squeeze(w(1:3,n-1,3:XHR-1)),dx1,1,g);
    dt2=setdt(squeeze(w(1:3,n-1,XHR:X+XHR-2)),dx2,1,g);
    %     check last time step
    if (t+dt2)>tmax
        dt2=tmax-t;
        t=t+dt2;
    else
        t=t+dt2;
    end
    n2=0;
    % First compute the region with smaller dx since dt is proportionally
    % smaller. Then move on to the second region.
    % this is just a temporary array to store until we reach t+dt2
    wold=w(:,n-1,:);
    %         inner loop to run until lower dx is at same time as higher dx
    %         region
    while(t1<t)
        %         similar to global loop
        if (t1+dt1)>t
            dt1=t-t1;
            t1=t1+dt1;
        else
            t1=t1+dt1;
        end
        tn(n)=t;
        for i=3:(XHR-1)
            %         downwind predicted
            [wpred2,dw2]=pred(wold(:,i),wold(:,i+1),wold(:,i+2),dt1,g);
            
            %         current predicted
            [wpred1,dw1]=pred(wold(:,i-1),wold(:,i),wold(:,i+1),dt1,g);
            
            %         upwind predicted
            [wpred0,dw0]=pred(wold(:,i-2),wold(:,i-1),wold(:,i),dt1,g);
            
            w2L=0.5*(wold(1:3,i)+wpred1(1:3)')...
                +(0.5*abs(wold(4,i+1)-wold(4,i)))*dw1;
            w2R=0.5*(wold(1:3,i+1)+wpred2(1:3)')...
                -(0.5*abs(wold(4,i+1)-wold(4,i)))*dw2;
            f2=w2f(rmannsol(w2L,w2R,g,0),g);
            
            
            w0L=0.5*(wold(1:3,i-1)+wpred0(1:3)')...
                +(0.5*abs(wold(4,i)-wold(4,i-1)))*dw0;
            w0R=0.5*(wold(1:3,i)+wpred1(1:3)')...
                -(0.5*abs(wold(4,i)-wold(4,i-1)))*dw1;
            f0=w2f(rmannsol(w0L,w0R,g,0),g);
            
            dx=(0.5*abs(wold(4,i+1)-wold(4,i)))...
                -(-0.5*abs(wold(4,i)-wold(4,i-1)));
            u1=w2u(wold(1:3,i),g);
            u2=u1-(dt1/dx)*(f2-f0);
            u2=u2w(u2,g);
            %         update primitive variables
            w(1:3,n,i)=u2(1:3);
            
        end
        % left side is solid wall        
        w(1:3,n,1)=w(1:3,n,4);
        w(1:3,n,2)=w(1:3,n,3);
        w(2,n,1:2)=-w(2,n,1:2);
        wold(1:3,1:XHR-1)=w(1:3,n,1:XHR-1);

        % boundary between different resolutions, same as adjacent cells
        wold(1:3,XHR)=wold(1:3,XHR-1);
        wold(1:3,XHR+1)=wold(1:3,XHR);
        wold(4,:)=w(4,n-1,:);
        n2=n2+1;
        
        if n2>timeout
            error('timeout')
        end
    end
    % now compute second region with higher dx
    for i=XHR:(X+XHR-2)
        %         downwind predicted
        [wpred2,dw2]=pred(w(:,n-1,i),w(:,n-1,i+1),w(:,n-1,i+2),dt2,g);
        
        %         current predicted
        [wpred1,dw1]=pred(w(:,n-1,i-1),w(:,n-1,i),w(:,n-1,i+1),dt2,g);
        
        %         upwind predicted
        [wpred0,dw0]=pred(w(:,n-1,i-2),w(:,n-1,i-1),w(:,n-1,i),dt2,g);
        
        w2L=0.5*(w(1:3,n-1,i)+wpred1(1:3)')...
            +(0.5*abs(w(4,n-1,i+1)-w(4,n-1,i)))*dw1;
        w2R=0.5*(w(1:3,n-1,i+1)+wpred2(1:3)')...
            -(0.5*abs(w(4,n-1,i+1)-w(4,n-1,i)))*dw2;
        f2=w2f(rmannsol(w2L,w2R,g,0),g);
        
        
        w0L=0.5*(w(1:3,n-1,i-1)+wpred0(1:3)')...
            +(0.5*abs(w(4,n-1,i)-w(4,n-1,i-1)))*dw0;
        w0R=0.5*(w(1:3,n-1,i)+wpred1(1:3)')...
            -(0.5*abs(w(4,n-1,i)-w(4,n-1,i-1)))*dw1;
        f0=w2f(rmannsol(w0L,w0R,g,0),g);
        
        dx=(0.5*abs(w(4,n-1,i+1)-w(4,n-1,i)))...
            -(-0.5*abs(w(4,n-1,i)-w(4,n-1,i-1)));
        u1=w2u(w(1:3,n-1,i),g);
        u2=u1-(dt2/dx)*(f2-f0);
        u2=u2w(u2,g);
        %         update primitive variables
        w(1:3,n,i)=u2(1:3);
    end
    %     boundary (solid wall)
    w(1:3,n,X+XHR)=w(1:3,n,X+XHR-3);
    w(1:3,n,X+XHR-1)=w(1:3,n,X+XHR-2);
    w(2,n,X+XHR-1:X+XHR)=-w(2,n,X+XHR-1:X+XHR);
    w(4,n,:)=w(4,n-1,:);
    %     live plotting
    %%%%%%%%%%%%%%%%%%%%
%     subplot(1,3,1)
    
%     plot(xmin:dx1:xmed,squeeze(w(1,n-1,3:(XHR)))',xmed:dx2:xmax,squeeze(w(1,n-1,(XHR):(X+XHR-2)))')
%     title(strcat('MUSCL scheme, wall interaction, t=',num2str(t)))
%     ylabel('Density')
%     xlabel('x')
%     subplot(1,3,2)
%     
%     plot(xmin:dx1:xmed,squeeze(w(2,n-1,3:(XHR)))',xmed:dx2:xmax,squeeze(w(2,n-1,(XHR):(X+XHR-2)))')
%     ylabel('Velocity')
%     xlabel('x')
%     subplot(1,3,3)
%     
%     plot(xmin:dx1:xmed,squeeze(w(3,n-1,3:(XHR)))',xmed:dx2:xmax,squeeze(w(3,n-1,(XHR):(X+XHR-2)))')
%     ylabel('Pressure')
%     xlabel('x')
%     legend('numerical solution','exact solution')
%     pause(0.02)
    %%%%%%%%%%%%%%%%%%%%%
    
    %     advance timestep
    n=n+1;
    %     check for timeout
    if n>timeout
        t
        error('timeout')
    end
end

subplot(1,3,1)
hold all

plot(xmin:dx1:xmed,squeeze(w(1,n-1,3:(XHR)))',xmed:dx2:xmax,squeeze(w(1,n-1,(XHR):(X+XHR-2)))')

title(strcat('MUSCL scheme, wall interaction, t=',num2str(t),''))
ylabel('Density')
xlabel('x')
subplot(1,3,2)
hold all

plot(xmin:dx1:xmed,squeeze(w(2,n-1,3:(XHR)))',xmed:dx2:xmax,squeeze(w(2,n-1,(XHR):(X+XHR-2)))')

ylabel('Velocity')
xlabel('x')
subplot(1,3,3)
hold all

plot(xmin:dx1:xmed,squeeze(w(3,n-1,3:(XHR)))',xmed:dx2:xmax,squeeze(w(3,n-1,(XHR):(X+XHR-2)))')

ylabel('Pressure')
xlabel('x')

% conversions between primitive, flux and conserved variables
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
    function [wp,dw]=pred(w0,w1,w2,dt,g)
        %         choose slope limiter beta is superbee for b=2, minmod for
        %         b=1,intermediate dissipation for 1<b<2
        dw=beta(w1-w0,w2-w1,1.9);
        w2l=w1(1:3)+0.5*(dw')*abs(w2(4)-w1(4));
        w0r=w1(1:3)-0.5*(dw')*abs(w1(4)-w0(4));
        fpl=w2f(w2l,g);
        fpr=w2f(w0r,g);
        fpr(4)=w1(4)-0.5*abs(w1(4)-w0(4));
        fpl(4)=w1(4)+0.5*abs(w2(4)-w1(4));
        
        up=w2u(w1(1:3),g)-(dt/abs(fpl(4)-fpr(4)))*(fpl(1:3)-fpr(1:3));
        
        dw=dw';
        wp=u2w(up,g);
        wp(4)=w1(4);
        
    end
    function [lim]=beta(w1,w2,b)
        lim=zeros(1,3);
        for q=1:3
            w1(q)=w1(q)/abs(w1(4));
            w2(q)=w2(q)/abs(w2(4));
            if w2(q)>=0
                lim(q)=max([0 min([b*w1(q),w2(q)]) min([w1(q),b*w2(q)])]);
            else
                lim(q)=min([0 max([b*w1(q),w2(q)]) max([w1(q),b*w2(q)])]);
            end
            
        end
    end

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