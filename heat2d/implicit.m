clear
close all
clc

xl = 0.4; %xlast in meter
x0 = 0; %x0 in meter
yl = 0.4; %rlast in meter
y0 = 0; %r0 in meter
dx = 0.02; %x step
dy = 0.02; %y step
lx = int16((xl - x0) / dx + 1); %length(x) with discretization
ly = int16((yl - y0) / dy + 1); %length(r) with discretization
y = y0:dy:yl; %creating r vector
x = x0:dx:xl; %creating x vector


%boundry conditions
% T(1, :) = 100;
% T(:, 1) = 80;

A = zeros(lx^2); %preallocating A
id=reshape(1:(lx^2),[lx, ly]); %creating an id matrix for determing constant positions in matrix
[s1,s2]=size(id);
B=zeros(lx^2, 1); %preallocating B

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
        
        B(ii) = 0;
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
T = A\B;
T = reshape(T, lx, ly);
[X, Y] = meshgrid(x, y);

surf(X, Y, T)
view(0, -90)
colormap cool
colorbar
xlabel('X')
ylabel('Y')
title('steady state')