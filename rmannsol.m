function [ wout ] = rmannsol(wl,wr,gamma,xt)
%rmannsol function takes primitive variables across a discontinuity to find
%the distribution of those variables after a time t
%   w=[density,velocity,pressure]
%   wr=variables on right of discontinuity
%   wl=variables on left of discontinuity
e=5;
% initializing the error
acc=1e-3;
% desired accuracy of the P* solution
right_speed=(gamma*wr(3)/wr(1))^0.5;
cl=(gamma*wl(3)/wl(1))^0.5;
% inital sound speed on both sides
p0=wl(3)*wr(1)*right_speed+wr(3)*wl(1)*cl+...
    (wl(2)-wr(2))*wr(1)*wl(1)*right_speed*cl;
p0=p0/(wr(1)*right_speed+wl(1)*cl);
% guess the initial pressure using the acoustic wave approximation
if p0<0
    p0=1e-3;
end
%in case the guess is negative fix the value
nit=0;
% initialize the counter
while e>acc
%     decide the type of wave (exp fan or shock)
    if p0<=wr(3)
%         right-facing expansion fan
        [fpr,dfpr]=efan(wr(1),p0,wr(3),gamma,right_speed);
    else
%         right-facing shock wave
        [fpr,dfpr]=shock(wr(1),p0,wr(3),gamma);
    end
    if p0<=wl(3)        
%         left-facing expansion fan
        [fpl,dfpl]=efan(wl(1),p0,wl(3),gamma,cl);
    else
%         left-facing shock wave
        [fpl,dfpl]=shock(wl(1),p0,wl(3),gamma);
    end
%     Use Newton-Raphson method to find P*
    pstar=p0-(fpr+fpl+(wr(2)-wl(2)))/(dfpr+dfpl);
%     update error
    e=abs((pstar-p0)/pstar);
%     update p* guess
    p0=pstar;
    nit=nit+1;
%     make sure the program does not infinite loop
    if nit>1e3
        error('timeout')
    end
end
% find final values of the functions to find u*
if pstar<=wr(3)
    [fpr,dfpr]=efan(wr(1),pstar,wr(3),gamma,right_speed);
else
    [fpr,dfpr]=shock(wr(1),pstar,wr(3),gamma);
end
if pstar<=wl(3)
    [fpl,dfpl]=efan(wl(1),pstar,wl(3),gamma,cl);
else
    [fpl,dfpl]=shock(wl(1),pstar,wl(3),gamma);
end

ustar=0.5*(wl(2)+wr(2))+0.5*(fpr-fpl);

% if x/t is in front of or behind the contact surface

if xt>=ustar;
    if pstar<=wr(3)
%         W*Right calculation for region behind right-facing expansion fan
        [rhosr,csr]=kfan(wr(1),pstar,wr(3),gamma,right_speed);
%         positions of head and tail of fan
        shr=wr(2)+right_speed;
        str=ustar+csr;
        if xt<=shr
            if xt<=str
                wout=[rhosr,ustar,pstar];
            else
%                 W*Right calculations for x/t inside fan
                wout=[0,0,0];
                wout(1)=wr(1)*...
                    (2/(gamma+1)-...
                    (gamma-1)*(wr(2)-xt)/...
                    ((gamma+1)*right_speed))^...
                    (2/(gamma-1));
                wout(2)=(2/(gamma+1))*(-right_speed+(gamma-1)*wr(2)/2+xt);
                wout(3)=wr(3)*...
                    ((2/(gamma+1))-...
                    (gamma-1)*(wr(2)-xt)/...
                    ((gamma+1)*right_speed))^...
                    ((2*gamma)/(gamma-1));
            end
        else
            wout=wr;
        end
    else
%         W*Right calculations for shock wave
        rhosr=kshock(wr(1),pstar,wr(3),gamma);
        sr=wr(2)+ac(pstar,wr(3),right_speed,gamma);
        if xt>sr
            wout=wr;
        else
            wout=[rhosr,ustar,pstar];
        end
    end
else
    if pstar<=wl(3)
%         W*Left calculations for region behind left-facing fan
        [rhosl,csl]=kfan(wl(1),pstar,wl(3),gamma,cl);
%         x/t of head and tail
        shl=wl(2)-cl;
        stl=ustar-csl;
        if xt>shl
            if xt>stl
                wout=[rhosl,ustar,pstar];
            else
%                 W*Left calculations for region inside fan
                wout=[0,0,0];
                wout(1)=wl(1)*...
                    (2/(gamma+1)+...
                    (gamma-1)*(wl(2)-xt)/...
                    ((gamma+1)*cl))^...
                    (2/(gamma-1));
                wout(2)=(2/(gamma+1))*(cl+(gamma-1)*wl(2)/2+xt);
                wout(3)=wl(3)*...
                    ((2/(gamma+1))+...
                    (gamma-1)*(wl(2)-xt)/...
                    ((gamma+1)*cl))^...
                    ((2*gamma)/(gamma-1));
            end
        else
            wout=wl;
        end
    else
%         W*Left for left-facing shock wave
        rhosl=kshock(wl(1),pstar,wl(3),gamma);
        sl=wl(2)-ac(pstar,wl(3),cl,gamma);
        if xt<sl
            wout=wl;
        else
            wout=[rhosl,ustar,pstar];
        end
    end
end


    function [fp,dfp]=shock(rhok,p,pk,gamma)
%         calculates f(P*) behind a shock wave
        fp=(p-pk)*((2/(rhok*(gamma+1)))/(p+pk*(gamma-1)/(gamma+1)))^0.5;
        dfp=(1-(p-pk)/(2*(p+pk*(gamma-1)/(gamma+1))))*...
            ((2/(rhok*(gamma+1)))/(p+pk*(gamma-1)/(gamma+1)))^0.5;
    end
    function [fp,dfp]=efan(rhok,p,pk,gamma,ck)
%         calculates f(P*) behind expansion fan
        fp=((p/pk)^((gamma-1)/(2*gamma))-1)*...
            (2*ck)/(gamma-1);
        dfp=(1/(rhok*ck))*(p/pk)^(-1*(gamma+1)/(2*gamma));
    end
    function rhosk=kshock(rhok,pstar,pk,gamma)
%         Calculates the parameters behind a shock
        rhosk=rhok*...
            ((pstar/pk)+(gamma-1)/(gamma+1))...
            /((gamma-1)*pstar/((gamma+1)*pk)+1);
    end
    function [rhosk,csk]=kfan(rhok,p,pk,gamma,ck)
%         Calculates the parameters behind an expansion fan
        pr=p/pk;
        rhosk=rhok*pr^(1/gamma);
        csk=ck*pr^((gamma-1)/(2*gamma));
    end
    function fps=ac(p,pk,ck,gamma)
%         Calculates the velocity of a shock wave in the fixed reference
%         frame
        fps=ck*((gamma+1)*p/(2*gamma*pk)+(gamma-1)/(2*gamma))^0.5;
        
    end
end

