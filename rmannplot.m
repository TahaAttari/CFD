% RSTest 1
% wl=[5.99924,19.5975,460.894];
% wr=[5.99242,-6.19633,49.0950];
% RSTest 2
% wl=[1,0.75,1];
% wr=[0.125,0,0.1];
% RSTest 3
% wl=[1,0,0.01];
% wr=[1,0,100];

wl=[1,-19.59745,1000];
wr=[1,-19.59745,0.01];

g=1.4;
x0=0.3;
t=0.012;
w=zeros(1000,4);
for i=1:1000
    x=-0.5+(1/1000)*(i-0.5);
    w(i,1)=x;
    x=x-x0;
    w(i,2:4)=rmannsol(wl,wr,g,x/t); 
end
subplot(3,1,1)
plot(w(:,1),w(:,2))
title('Density')
xlabel('x')
subplot(3,1,2)
plot(w(:,1),w(:,3))
title('Velocity')
xlabel('x')
subplot(3,1,3)
plot(w(:,1),w(:,4))
title('Pressure')
xlabel('x')