clear
close all
clc

%% constants and allocations
dx = 0.1;
dt = 0.1;
alpha = 2.5e-2;
beta = 0.1;
l = 1;
t = 20;
x = 0:dx:l;
t = 0:dt:t;
Ta = 25;

T = 25 * ones(length(x),length(t));
T(1, :) = 100;
error = 1;n = 1;
%% solving for the first time
for i = 2:length(x) - 1
        dTdt = alpha / dx^2 * (T(i+1, n) - 2 * T(i, n) + T(i-1, n)) ...
            - beta * (T(i, n) - Ta);
        
        T(i, n+1) = dt * dTdt + T(i, n);
end
while error > 1e-4
    n = n + 1;
    for i = 2:length(x)-1
        %method of line ->
        dTdt = alpha / dx^2 * (T(i+1, n) - 2 * T(i, n) + T(i-1, n)) ...
            - beta * (T(i, n) - Ta);
        %solving with euler method ->
        T(i, n+1) = dt * dTdt + T(i, n);
    end
    error = (sum(T(:, n)) - sum(T(:, n-1))) / sum(T(:, n));
end
T(:,n:end) = [];
fprintf("elapsed time to reach steady state is %f \n\n", t(n));
hold on
plot(x,T(:,1))
title("T vs X")
plot(x,T(:,end))
legend("t = 0","t = tlast")
