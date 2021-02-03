clear
close all


R=0.1;L=0.9;
z=0:0.001:L;
r=0:0.01:R;
T0=800;
K=6.25;Ru=50;Cp=4;v=10;
A=Ru*v*Cp/K;
x=besselzero(0,100,1); %bessel j0 first 5 roots
lambda=x/R;
t=zeros(length(z),length(r));
for i=1:length(z)
    for j=1:length(r)
        An=2/R./lambda*T0./besselj(1,lambda*R);
        T=An.*besselj(0,lambda.*r(j)).*exp(-lambda.^2*z(i)/A);
        t(i,j)=sum(T);
    end
end
plot(z,t)
figure
contourf(z,r,t')