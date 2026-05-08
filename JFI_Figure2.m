clc
clearvars
close all

%% Signal 1
Fs = 1000;
t = 0:1/Fs:1;
signal1 = sin(2*pi*40*t) + 0.1*sin(2*pi*30*t);
figure('DefaultAxesFontSize',14);
subplot(5,1,1)
plot(t,signal1,'b','LineWidth',1)
title('(b)','FontWeight','Normal')

%% Signal 2
Fs = 1000;
t = 0:1/Fs:1;
N = length(t);
signal_comps2 = zeros(N,2);
signal_comps2(:,1) = 0.9*sin(2*pi*(160 - (30/2)*t).*t);
signal_comps2(:,2) = 0.7*sin(2*pi*(240 - (30/2)*t).*t);
signal2 =  sum(signal_comps2,2);
subplot(5,1,2)
plot(t,signal2,'b','LineWidth',1)
title('(b)','FontWeight','Normal')

%% Signal 3
Fs = 1000;
t = 0:1/Fs:1;
N = length(t);
signal_comps3 = zeros(N,2);
signal_comps3(:,1) = 0.85*sin(2*pi*(190 + (80/2)*t).*t);
signal_comps3(:,2) = 0.6*sin(2*pi*(350 - (180/2)*t).*t);
signal3 =  sum(signal_comps3,2);
subplot(5,1,3)
plot(t,signal3,'b','LineWidth',1)
ylabel('Amplitude')
title('(c)','FontWeight','Normal')

%% Signal 4
%--------------------------------------------------------------------------
% Please download the cello G5 note signal from the link provided below.
% Link: https://github.com/dfourer/ASTRES_toolbox (new link - active)
% Old link: https://github.com/dfourer/ASTREStoolbox
%--------------------------------------------------------------------------
[signal4,Fs] = audioread("cello_short.wav");
N = length(signal4);
t = 0:1/Fs:(N-1)/Fs;
subplot(5,1,4)
plot(t,signal4,'b','LineWidth',1)
xlim([t(1) t(end)])
title('(d)','FontWeight','Normal')

%% Signal 5
%--------------------------------------------------------------------------
% Please download the bat signal from the link provided below.
% Link: https://www.ece.rice.edu/dsp/software/bat.shtml
%--------------------------------------------------------------------------
load("batsignal.mat")
signal5 = Bat_signal;
Fs = 1/7*1e6 ;
N = length(signal5);
t = 0:1/Fs:(N-1)/Fs;
subplot(5,1,5)
plot(t,signal5,'b','LineWidth',1)
xlim([t(1) t(end)])
xlabel('Time (s)')
title('(e)','FontWeight','Normal')