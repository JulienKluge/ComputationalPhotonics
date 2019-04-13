clc
clear all
close all

gauss1D = @(x) exp(-(x * x));
gauss2D = @(x, y) exp(-(x * x + y * y));

E1x = @(x, y) -y;
E1y = @(x, y) x;
E2x = @(x, y) -y;
E2y = @(x, y) -x;
E2 = @(x, y) -[y, x];

h = @(x) exp(1 - 1 ./ x) ./ x;

%% 1a
xMax = 2;
dx = 2 * xMax / 100;
x = -xMax:dx:xMax;
y = arrayfun(gauss1D, x);
figure(1)
plot(x, y)
saveas(gcf, 'A1_a.png');

%% 1b
gMax = 2;
dg = 2 * gMax / 50;
g = -gMax:dg:gMax;
[x, y] = meshgrid(g, g);
z = arrayfun(gauss2D, x, y);
figure(2);
subplot(2, 2, 1);
mesh(x, y, z)
subplot(2, 2, 2);
surface(x, y, z)
subplot(2, 2, 3);
pcolor(x, y, z)
subplot(2, 2, 4);
contour(x, y, z)
saveas(gcf, 'A1_b.png');

%% 1c
gMax = 1;
dg = 2 * gMax / 20;
g = -gMax:dg:gMax;
[x, y] = meshgrid(g, g);
xp1 = arrayfun(E1x, x, y);
yp1 = arrayfun(E1y, x, y);
xp2 = arrayfun(E2x, x, y);
yp2 = arrayfun(E2y, x, y);

figure(3);
set(gcf, 'position', [10, 100, 1000, 400]);
subplot(1, 2, 1);
quiver(x, y, xp1, yp1)
subplot(1, 2, 2);
quiver(x, y, xp2, yp2)
saveas(gcf, 'A1_c.png');

%% 1d
dx = 0.01;
x = 0:dx:1;
y = h(x);
figure(4);
plot(x, y)
saveas(gcf, 'A1_d.png');
int1 = integral(h, 0, 1);
int2 = integral(h, eps, 1);
h0 = h(0);
heps = h(eps);
fprintf('h(0) = %f || h(eps) = %f', h0, heps)
