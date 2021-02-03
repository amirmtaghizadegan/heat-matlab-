clear
close all
clc
K=6.25;Ru=50;Cp=4; 
dx=0.001; 
dr=0.01; 
v=10; 
r=0:dr:0.1;
x=0:dx:0.9;
T=ones(length(x),length(r))*273;
T(1,:)=800+273;
% gamma=K*dx/Ru/v/Cp/dr; %gamma
%boundry conditions
T(1,:)=800+273;
T(:,end)=20+273;
for i=1:length(x)-1
    k=0;%equation counter
    zarayeb=zeros(length(r));javab=zeros(length(r),1);
    for j=2:length(r)-1
        k=k+1;
        gamma=dx/160/(1-(r(k)^2/0.9^2))/dr; %gamma
        zarayeb(k,j+1)=-0.5*gamma*(1/r(j)+1/dr);
        zarayeb(k,j-1)=-0.5*gamma*(1/dr);
        zarayeb(k,j)=0.5*gamma*(1/r(j)+2/dr)+1;
        T(i,1)=T(i,2);
        T(1,:)=800+273;
        T(:,end)=20+273;
        javab(k)=T(i,j)+0.5*gamma*((1/r(j)+1/dr)*T(i,j+1)+(1/dr)*T(i,j-1)-(1/r(j)+2/dr)*T(i,j));
    end
    k=k+1;
    zarayeb(k,1)=1;
    javab(k)=T(i,2);
    k=k+1;
    zarayeb(k,end)=1;
    javab(k)=T(i,end);
    T(i+1,:)=n_majhool(zarayeb,javab);
end

[R,X]=meshgrid(r,x);
mesh(R,X,T)
figure
plot(x,T)

figure
% T=[fliplr(T(:,2:end)),T];
% r=-0.1:dr:0.1;
contourf(x,r,T')