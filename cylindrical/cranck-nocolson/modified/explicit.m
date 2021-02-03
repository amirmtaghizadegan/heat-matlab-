clear
close all
clc

K=6.25;Ru=50;Cp=4; 
dx=0.0001; 
dr=0.001; 
v=10; 
gamma=K*dx/Ru/v/Cp/dr; %gamma
r=0:dr:0.1;
x=0:dx:0.9;
T=ones(length(x),length(r))*273;
T(1,:)=800+273;
for i=1:length(x)-1
    for j=2:length(r)-1
        T(i,1)=T(i,2);
        T(i+1,j)= gamma.*(1./r(j).*(T(i,j+1)-T(i,j))+1./dr.*(T(i,j+1)-2.*T(i,j)+T(i,j-1)))+T(i,j);
    end
end
T(:,end)=20+273;
T(end,1)=T(end,2);

[R,X]=meshgrid(r,x);
mesh(R,X,T)

figure
plot(x,T)

figure
T=[fliplr(T(:,2:end)),T];
r=-0.1:dr:0.1;
contourf(x,r,T')