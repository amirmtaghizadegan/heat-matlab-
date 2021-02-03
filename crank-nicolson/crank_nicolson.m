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
lambda=a*dt/dx^2; %lambda
lx=(xl-x0)/dx+1; %length(x) with discretization
ly=(yl-y0)/dy+1; %length(y) with discretization
lt=Tl/dt+1; %length(t) with discretization
t=ones(lx,ly,lt)*273; %prealocating
f0=h*tinf; %f0 is a factor that can be considered as h*tinf

for n=1:lt-1
    A=zeros(lx^2); %prealocating A
    id=reshape(1:(lx^2),[lx ly]); %creating an id matrix for determing constant positions in matrix
    [s1,s2]=size(id);
    ii=0; %equation counter
    B=zeros(1,lx^2); %prealcating B
    for j=2:ly-1
        for i=2:lx-1
            ii=ii+1;
            A(ii,id(i+1,j))=-lambda*0.5;
            A(ii,id(i,j+1))=-lambda*0.5;
            A(ii,id(i-1,j))=-lambda*0.5;
            A(ii,id(i,j-1))=-lambda*0.5;
            A(ii,id(i,j))=2*lambda+1;
            t(i,2,n)=(f0*dy/k)-(h*dy/k*t(i,1,n))+(t(i,1,n));
            t(:,:,1)=273;
            t(:,end,:)=273;
            t(1,:,:)=273;
            t(end,:,:)=273;
            B(ii)=t(i,j,n)+0.5*lambda*(t(i+1,j,n)+t(i,j+1,n)+t(i-1,j,n)+t(i,j-1,n)-4*t(i,j,n));
        end
    end
    %% Alocating initial and boundry conditions to A and B
    T=t(:,:,n);
    idd=id([1,end],:);idd=idd';
    id([1,end],:)=[];
    for jj=1:numel(idd)
        ii=ii+1;
        A(ii,idd(jj))=1;
        B(ii)=T(idd(jj));
    end
    idd=id(:,[1,end]);id([1,end],:)=[];
    for jj=1:numel(idd)
        ii=ii+1;
        A(ii,idd(jj))=1;
        B(ii)=T(idd(jj));
    end
    %% calculating Temprature for time(n+1)
    T=n_majhool(A,B');
    t(:,:,n+1)=reshape(T,[lx ly]);
end
% extending x,y and time for surfing
x=x0:dx:xl;
y=y0:dy:yl;
T=0:dt:Tl;
[x,y]=meshgrid(x,y);
%surfing T for last time node
Tend=t(:,:,end);
surf(x,y,Tend)
save cranck.dat Tend -ascii