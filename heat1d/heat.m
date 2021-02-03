clear
close all
clc

l=1;dx=0.025;dt=0.1;alpha=1;
x=0:dx:l;
t=0:dt:0.3;
T(1,:)=100*sin(pi*x/l);pl=zeros(1,length(t));
pl(1)=plot(x,T(1,:));
hold on
for p=1:length(t)-1
for m=2:length(x)-1
T(p+1,m)=alpha*dt/dx^2*(T(p,m+1)+T(p,m-1))+(1-2*alpha*dt/dx^2)*T(p,m);
end
pl(p+1)=plot(x,T(p+1,:));
end
legend([pl(1),pl(end)],'t=0','t=0.3')
xlabel('x')
ylabel('T')
fprintf('the final answere is:\n')
disp(T)