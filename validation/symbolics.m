close all; clear all; clc;

syms z T Y U a3 a2 a1 a0 b2 b1 b0


% % bilinear trafo
% z = (1 + s * T/2) / (1 - s * T/2);
% s = 2/T * (z-1) / (z+1);

% backward euler
% z = 1 / (1 - T * s);
s = (z-1) / (z * T);



G = (b2 * s^2 + b1 * s + b0) / (a3 * s^3 + a2 * s^2 + a1 * s + a0)

% G = 1 / z;
Y = U * G;


y = iztrans(Y)
% 
% t = linspace(0, 1, 100);
% u = cos(pi * t);
% y = zeros(100,1);
% for i = 2:100
%     y(i) =  (u(i) - u(i-1)) / (t(i)-t(i-1));
% end
% 
% hold on
% plot(t, y)
% plot(t, u, "--k")
% hold off