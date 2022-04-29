clear all; close all;
%Add the path where functions are stored
addpath("ImpulseResponse")
%Duration of impulse response
duration = 1;
%Sample rate
FS = 1e5;
%Cylinder, 1m, 0.005 m radius bore
%Define bore profile
cylinder.x = [0, 1000];
cylinder.S = pi*0.005^2*[1,1];
cylinder.temp = 19.85;

irCyl = TubeImpulseResponse(cylinder, duration, FS);
impedanceCyl = abs(fft(irCyl));
n = length(irCyl);
figure(1)
subplot(2,1,1),plot([0:n-1]/(n*FS), irCyl,'LineWidth',2);
xlabel("Time(s)","Fontsize",12)
ylabel("Pressure (Pa)","Fontsize",12)
title("Impulse response","Fontsize",15)
subplot(2,1,2),semilogy([0:n-1]*FS/n, impedanceCyl,'LineWidth',2)
xlim([0 10000])
xlabel("Frequency (Hz)","Fontsize",12)
title("Input impedance magnitude","Fontsize",15)
