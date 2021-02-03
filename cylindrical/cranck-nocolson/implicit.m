clear
close all
clc

ro=50;
cp=4;
xl=0.9; %xlast in meter
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
t=zeros(lx,lr); %prealocating
%boundry conditions
t(1,:)=800+273;
t(:,end)=20+273;
for i=1:lx-1
    ii=0;%equation counter
    A=zeros(lr);B=zeros(lr,1);
    for k=2:lr-1
        ii=ii+1;
        A(ii,k+1)=-lambda*(1/r(k)+1/dr);
        A(ii,k-1)=-lambda*(1/dr);
        A(ii,k)=lambda*(1/r(k)+2/dr)+1;
        
        t(i,1)=t(i,2);
        t(1,:)=800+273;
        t(:,end)=20+273;
        
        B(ii)=t(i,k);
    end
    %--symmetry condition
    ii=ii+1;
    A(ii,1)=1;
    B(ii)=t(i,2);
    %--
    ii=ii+1;
    A(ii,end)=1;
    B(ii)=t(i,end);
    t(i+1,:)=A\B;
end

t=[fliplr(t(:,2:end)),t];
r=-rl:dr:rl;
disp(t)
contourf(x,r,t')