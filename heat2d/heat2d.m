clear
close all
clc
alpha=1/4;
dt=0.01;
dx=0.1;
t=1;
x=0:dx:1;y=x;t=0:dt:t;
T=zeros(length(y),length(x),length(t));
[x,y]=meshgrid(x,y);
%ic and bc

for i=1:length(t)
T(size(x,2),:,i)=100*(sin(1.*y(:,end)*pi/2));
T(:,size(y,1),i)=100*(sin(x(end,:).*1*pi/2));
end
T(:,:,1)=100*(sin(x.*y*pi/2));
%pde solve
surf(x,y,T(:,:,1))
hold on
for p=1:length(t)-1
    for m=2:length(x)-1
        for n=2:length(y)-1
            T(n,m,p+1)=(alpha*dt/dx^2)*(T(n,m+1,p)+T(n,m-1,p)+T(n+1,m,p)+T(n-1,m,p))+(1-4*alpha*dt/dx^2)*T(n,m,p);
        end
    end
    surf(x,y,T(:,:,p+1))
end
xlabel('x')
ylabel('y')
zlabel('T')
disp(T)