close all
clear
clc

%%
a=2.3e-5; %alpha in m^2/s
k=200; %(W/m k)
h=100; %W/(m2 K)
tinf=400; %T(inf)
Tl=1000; %Tlast in secs
xl=1; %xlast in meter
x0=0; %x0 in meter
yl=1; %ylast in meter
y0=0; %y0 in meter
dx=0.1; %x step
dy=0.1; %y step
dt=1; %time step
lt=Tl/dt+1; %length(t) with discretization
lx=(xl-x0)/dx+1; %length(x) with discretization
ly=(yl-y0)/dy+1; %length(y) with discretization
Td=a*Tl/(dx^2+dy^2); %dimensionless time, at/(x^2+y^2);

%%
x=x0:dx:xl;
y=y0:dy:yl;
t=0:dt:Tl;
% u=0:0.00001:Tl;u(1)=eps;
T=zeros(lx,ly);
for i=2:lx
    for j=1:ly
        j21=integral(@(u) exp((y(j)./sqrt(4*a*u)+h*(sqrt(a*u))/k).^2).*erfc(y(j)./sqrt(4*a*u)+h*(sqrt(a*u))/k).*exp(-(x(i).^2+y(j).^2)./(4*a*u))./u.^(1.5),0,Tl);
%         j21=trapz(u,j21);
        i5=0.5*integral(@(u) erfc(x(i)./sqrt(4*a*u)).*exp(-y(j)^2./(4*a.*u))./u.^(1.5),0,Tl);
%         i5=0.5*trapz(u,i5);
        t=t(end);
        T(i,j)=tinf*erf(x(i)./sqrt(4.*a.*t)).*(erfc(y(j)./sqrt(4*a.*t))-exp(-y(j).^2/(4.*a.*t)).*exp((y(j)./sqrt(4*a*t)+h*sqrt(a*t)/k).^2).*erfc(y(j)./sqrt(4*a*t)+h*sqrt(a*t)/k))...
        -tinf.*x(i)./sqrt(pi)./sqrt(4*a).*(j21-2*i5);
    end
end
[x,y]=meshgrid(x,y);
surf(y,x,T)
save analytic_int.dat T -ascii







