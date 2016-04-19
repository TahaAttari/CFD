function linear2ndOrd
% mesh refining
min=3;
max=12;
step=3;
for cell=min:step:max
    %     plotting
    subplot((max/step),1,cell/step)
    %     discretization and parameter setting
    a=2;
    xmin=0;
    xmax=300;
    tmax=100;
    %     cell size control
    CF=cell/step;
    N=CF*20+2;
    dt=tmax/(N-2);
    %     constant courant number (of 1)
    dx=a*dt;
    X=(xmax/dx)+1;
    dx=(xmax-xmin)/(X-1);
    X=round(X);
    %     inital disturbance
    disturb=[0 20];
    
    u=zeros(N,X);
    %     the desired distrubance type is uncommented
    for j=1:(disturb(2)/dx)
        %         triangle
        u(1,j)=0.075*j*dx+0.5;
        %         smooth
        %                 u(1,j)=1-0.5*cos(pi*dx*j/10);
    end
    u(1,round(disturb(2)/dx):X)=0.5;
    disturb=disturb+tmax*a;
%     exact solution
    exact=linspace(xmin,xmax,1000);
    uexact(1:round(disturb(1)/(xmax/1000)))=0.5;
    uexact(round(disturb(2)/(xmax/1000)):1000)=0.5;
    for j=(round((disturb(1)/(xmax/1000)))):round((disturb(2)/(xmax/1000)))
        uexact(j)=0.075*(j*(xmax/1000)-disturb(1))+0.5;
%         uexact(j)=1-0.5*cos(pi*((xmax/1000)*j-disturb(1))/10);
    end
    
    for n=2:N
        %     boundary
        u(n-1,1)=0.5;
        for i=2:(X-1)
            %         downwind flux
            f2=flux(u(n-1,i),u(n-1,i+1),dt,dx,a);
            %         upwind flux
            f1=flux(u(n-1,i-1),u(n-1,i),dt,dx,a);
            
            u1=u(n-1,i);
            u(n,i)=u1-(dt/dx)*(f2-f1);
            
        end
        %     boundary
        u(n,X)=0.5;
    end
    if size(u,2)~=size((xmin:dx:xmax),2)
        X=X-1;
    end
    hold on
    plot(xmin:dx:xmax,u(n,1:X)','x')
    plot(exact,uexact)
    xlabel('x')
    ylabel('u')
    title(strcat('2nd Order Triangle number of steps in x: X=',num2str(X-1)))
    legend('numerical solution','exact solution')
    hold off
%     error calulator
    for i=2:X
        err(i)=u(n,i)-uexact(round((i-1)*dx/(xmax/1000)));
        err(i)=abs(err(i));
    end
%     l2 norm of error
    l2(cell/step)=norm(err);
%     inverse of dx
    dxinv(cell/step)=(1/dx);
end
figure
% fit the l2 norm
afit=fit(log10(dxinv)',log10(l2)','poly1')
% plot the norm
plot(afit,log10(dxinv),log10(l2))
xlabel('log( 1/dx )')
ylabel('log( |error| )')
title('2nd Order Global Convergence Analysis')

    function [fluxo]=flux(w1,w2,dt,dx,a)
        
        fluxo=0.5*(w1+w2)-0.5*a*(dt/dx)*(w2-w1);
        fluxo=a*fluxo;
    end

end