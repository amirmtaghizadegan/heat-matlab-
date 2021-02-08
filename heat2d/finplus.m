clear
close all
clc

Tfin = inf;
fin.l = 25e-3; %m
while Tfin(end) > 25
fin.l = fin.l + 1e-2;

plate.k = 43; %w/m/k
plate.ro = 0.473; %j/g/k
plate.cp = 0.234; %j/g/k
qdot = 500; %W/m2
xl = 0.4; %xlast in meter
x0 = 0; %x0 in meter
yl = 0.4; %rlast in meter
y0 = 0; %r0 in meter
dx = 0.02; %x step
dy = 0.02; %y step
lx = int16((xl - x0) / dx + 1); %length(x) with discretization
ly = int16((yl - y0) / dy + 1); %length(y) with discretization
y = y0:dy:yl; %creating r vector
x = x0:dx:xl; %creating x vector
lambda = plate.k/dx^2;

%fin setup
fin.h = 15; %w/m2/k
Tinf = 20; %°C
fin.k = 419; %W/m/k
fin.ro = 10525; %kg/m3
fin.cp = 0.234; %j/g.k
fin.thick = 2e-3; %m

dz = 0.01; %z step
z0 = 0; %m
zl = fin.l; %m
z = z0:dz:zl; %creating z vector
lz = int16((zl - z0) / dz + 1); %length(z) with discretization
fin.A = fin.thick^2;
fin.p = fin.thick * fin.l;
%boundry conditions
% T(1, :) = 100;
% T(:, 1) = 80;

A = zeros(lx^2); %preallocating A
id=reshape(1:(lx^2),[lx, ly]); %creating an id matrix for determing constant positions in matrix
[s1,s2]=size(id);
B=zeros(lx^2, 1); %preallocating B
id_center = [lx/2, ly/2];
%entering equations
ii = 0;%equation counter
for i = 2:lx-1
    for k = 2:ly-1
        ii = ii + 1;
        A(ii, id(i+1, k)) = 1;
        A(ii, id(i, k+1)) = 1;
        A(ii, id(i-1, k)) = 1;
        A(ii, id(i, k-1)) = 1;
        A(ii, id(i, k)) = -4;
        if i == 11
            if k == 11
                A(ii, id(i+1, k)) = 1;
                A(ii, id(i, k+1)) = 1;
                A(ii, id(i-1, k)) = 1;
                A(ii, id(i, k-1)) = 1;
                A(ii, id(i, k)) = -4+1/dz;
                %A(ii,tfin) = -1/dz
                ii_center = ii;
            end
        end
        
        
        B(ii) = qdot/lambda;
    end
end
for i = 1:lx
    %boundry conditions
    % T(:, 1) = 80;
    ii = ii + 1;
    A(ii, id(i, 1)) = 1;
    B(ii) = 80;
    % d(T)/dx = 0 at x = xlast
    ii = ii + 1;
    A(ii, id(i, end)) = 1;
    A(ii, id(i, end-1)) = -1;
    B(ii) = 0;
end
for i = 1:ly
    %boundry conditions
    % T(1, :) = 100;
    ii = ii + 1;
    A(ii, id(1, i)) = 1;
    B(ii) = 100;
    % d(T)/dx = 0 at x = xlast
    ii = ii + 1;
    A(ii, id(end, i)) = 1;
    A(ii, id(end-1 , i)) = -1;
    B(ii) = 0;
end
%now making the fin equations
B2 = zeros(lz, 1);
%concatinating with the previous A and B
A(end+lz, end+lz) = 0;
B = [B; B2];
last = lx^2;
for m = 2:lz-1 
    ii = ii + 1;
    A(ii, last+m) = -2-(fin.h*fin.p/fin.k/fin.A)*dz^2;
    A(ii, last + m+1) = 1;
    A(ii, last + m-1) = 1;
    B(ii) = -(fin.h*fin.p/fin.k/fin.A)*dz^2 * Tinf;
end
% fin boundry conditions
ii = ii + 1;
A(ii, last + 1) = -1;
A(ii, id(id_center)) = 1;
B(ii) = 0;

ii = ii + 1;
A(ii, end) = 1-fin.h/fin.k*dz;
A(ii, end-1) = -1;
B(ii) = -fin.h/fin.k*dz*Tinf;
%%%%%%%%%
%coupling the two equations
A(ii_center, last+1) = -1/dz;
%%%%%%%%%%%%%%
T = A\B;
Tfin = T(end-lz+1:end);
T(end-lz+1:end) = [];
T = reshape(T, lx, ly);
[X, Y] = meshgrid(x, y);
end
surf(X, Y, T)
view(0, -90)
colormap cool
colorbar
xlabel('X')
ylabel('Y')
title('steady state')
