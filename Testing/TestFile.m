% Script to test  TubeImpulseResponse function. Reference impedances are
% calculated using Planar Transmission Matrix Method terminated with 
% Levine and Schwinger Planar radiation impedance. Visctherm losses are
% included using the Zwikker and Kosten loss model. 1000 segments are 
% along the length of the tube using the TMM.
clear all; close all;
%Add the path where functions are stored
addpath("../ImpulseResponse")
%Duration of impulse response
duration = 10;
%Sample rate
SR = 1e5;
%Cylinder, 1m, 0.01 cm radius bore, at temp 20C
%Define bore profile
load("Cylinder_TMM.mat");
impedanceRef = abs(Z);
fRef = freq;
cylinder.x = bore(:,1);
cylinder.r = bore(:,2);
cylinder.temp = temp;

irCyl = TubeImpulseResponse(cylinder, duration, SR);
impedanceCyl = abs(fft(irCyl));
n = length(irCyl);
figure(1)
subplot(2,1,1),semilogy([0:n-1]*SR/n, impedanceCyl,fRef, impedanceRef,'LineWidth',2);
legend("Code","Reference");
xlim([0 10000])
xlabel("Frequency (Hz)","Fontsize",12)
title("Cylinder input impedance magnitude","Fontsize",15)

%Straight cone, 1m, 0.01 to 0.1 cm radius bore, at temp 20C
%Define bore profile
load("Cone_TMM.mat");
impedanceRef = abs(Z);
fRef = freq;
cone.x = bore(:,1);
cone.r = bore(:,2);
cone.temp = temp;
%
irCone = TubeImpulseResponse(cone, duration, SR);
impedanceCone = abs(fft(irCone));
n = length(irCone);
subplot(2,1,2),
semilogy([0:n-1]*SR/n, impedanceCone,fRef, impedanceRef,'LineWidth',2);
legend("Code","Reference");
xlim([0 10000])
xlabel("Frequency (Hz)","Fontsize",12)
title("Cone input impedance magnitude","Fontsize",15)
