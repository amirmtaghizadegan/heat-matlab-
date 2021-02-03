clear
close all
clc

load cranck.dat
load implicit.dat
load explicit.dat
load analytic.dat
analytic=rot90(analytic)+273;

analytic(:,1)=[];
implicit(:,1)=[];
explicit(:,1)=[];
cranck(:,1)=[];

analytic(:,end)=[];
implicit(:,end)=[];
explicit(:,end)=[];
cranck(:,end)=[];

analytic(1,:)=[];
implicit(1,:)=[];
explicit(1,:)=[];
cranck(1,:)=[];

analytic(end,:)=[];
implicit(end,:)=[];
explicit(end,:)=[];
cranck(end,:)=[];

Error_explicit=abs(analytic-explicit);
sError_explicit=sum(sum(Error_explicit));
Error_implicit=abs(analytic-implicit);
sError_implicit=sum(sum(Error_implicit));
Error_cranck=abs(analytic-cranck);
sError_cranck=sum(sum(Error_cranck));

fprintf('Error_explicit = %d \n',sError_explicit)
fprintf('Error_implicit = %d \n',sError_implicit)
fprintf('Error_cranck = %d \n',sError_cranck)