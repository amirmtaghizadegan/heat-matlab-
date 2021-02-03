clear
close all
clc

ro=50;
cp=4;
xl=.9; %xlast in meter
x0=0; %x0 in meter
rl=0.1; %rlast in meter
r0=0; %r0 in meter
dx=0.001; %x step
dr=0.005; %y step
u=10; %u in m/s
k=6.25; %(W/m k)
lambda=k*dx/ro/u/cp/dr; %lambda
lx=(xl-x0)/dx+1; %length(x) with discretization
lr=(rl-r0)/dr+1; %length(r) with discretization
r=r0:dr:rl; %creating r vector
x=x0:dx:xl; %creating x vector
t=ones(lx,lr)*273; %prealocating
for i=1:lx-1
    for k=2:lr-1
        t(1,:)=800+273;
        t(:,end)=20+273;
        t(i,1)=t(i,2);
        t(i+1,k)= lambda*(1/r(k)*(t(i,k+1)-t(i,k))+1/dr*(t(i,k+1)-2*t(i,k)+t(i,k-1)))+t(i,k);
    end
end
t(end,1)=t(end,2);

t=[fliplr(t(:,2:end)),t];
r=-rl:dr:rl;
disp(t)
contourf(x,r,t')
contourf(x,r,t')