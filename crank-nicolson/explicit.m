clear
close all
clc

xl=1; %xlast in meter
x0=0; %x0 in meter
yl=1; %ylast in meter
y0=0; %y0 in meter
Tl=1000; %Tlast in secs
dt=1; %time step
dx=0.1; %x step
dy=0.1; %y step
a=2.3e-5; %alpha in m^2/s
k=200; %(W/m k)
h=100; %W/(m2 K)
tinf=400; %T(inf)
Ax=a*dt/dx^2; %lambda for x
Ay=a*dt/dy^2; %lambda for y
lx=(xl-x0)/dx+1; %length(x) with discretization
ly=(yl-y0)/dy+1; %length(y) with discretization
lt=Tl/dt+1; %length(t) with discretization
t=ones(lx,ly,lt)*273; %prealocating
f0=h*tinf; %f0 is a factor that can be considered as h*tinf
for n=1:lt-1
    for j=2:ly-1
        for i=2:lx-1
            %Boundry and initial conditions
            t(i,2,n)=(f0*dy/k)-(h*dy/k*t(i,1,n))+(t(i,1,n));
            t(:,:,1)=273;
            t(:,end,:)=273;
            t(1,:,:)=273;
            t(end,:,:)=273;
            % calculation
            t(i,j,n+1)=Ax*(t(i+1,j,n)-2*t(i,j,n)+t(i-1,j,n))+Ay*(t(i,j+1,n)-2*t(i,j,n)+t(i,j-1,n))+t(i,j,n);
        end
    end
end
%extending x,y and time for surfing
x=x0:dx:xl;
y=y0:dy:yl;
T=0:dt:Tl;
[x,y]=meshgrid(x,y);
%surfing for the last time node(t=Tl)
Tend=t(:,:,end);
surf(x,y,Tend)
save explicit.dat Tend -ascii