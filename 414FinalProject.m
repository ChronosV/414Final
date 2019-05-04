%Clear things up
clear all; clc; 
clf;
%Get a TF
s = tf('s');

%Create constants
Jg = 6.2e-6;    %Gearbox inertia
Js = 1.4e-7;    %Angular sensor inertia
Gv = 5;         %Voltage amplifier gain
Ks = 10;    %touch sensor, in V/m
Ra = 1e-2;      %Radius of ball, in meters
Rh = 6e-3;      %Height difference between center of ball and top of rail, in meters
Jb = 4.52185e-7; %Ball inertia
G = 9.81;       %Acceleration of gravity
N = 10;         %Chosen Gear Ratio
