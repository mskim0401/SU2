clear all
clc
% close all
format long
Re = 2.0*10^6;
rho = 1.2886;
l = 0.60;
mu = 1.853*10^-5
V = Re*mu/rho/l
AoA = 17.23; %deg
Vx = V*cosd(AoA)
Vy = V*sind(AoA)