clc
clearvars
close all

%% Loading the real-life bat signal
%--------------------------------------------------------------------------
% Please download the bat signal from the link provided below.
% Link: https://www.ece.rice.edu/dsp/software/bat.shtml
%--------------------------------------------------------------------------
load("batsignal.mat")
signal = Bat_signal;
Fs = 1/7*1e6;   % Sampling interval is 7 microseconds.

%% Defining the window length, and overlap length 
N = length(signal);
t = 0:1/Fs:(N-1)/Fs;
window_length = 79;
overlap_length = window_length - 1;
sample_shift = window_length - overlap_length;
Nframes = floor((N - window_length)/sample_shift) + 1;

%% MDL-based proposed method
f = (0:1:Fs/2);
V1 = zeros(length(f),Nframes);
time_stamp = zeros(1,Nframes);
for fr = 1:Nframes
    sig = signal((fr-1)*sample_shift+1:(fr-1)*sample_shift+window_length);
    time_stamp(fr) = (1/Fs)*(floor(((fr-1)*sample_shift+1+(fr-1)*...
        sample_shift+window_length)/2)+1);
    K = (window_length + 1)/2;
    X = hankel(sig(1:K),sig(K:end));
    [U,V] = eig(X);
    N_opt = MDL_sig_ev(U,V,sig);
    V2 = abs(diag(V));
    comp = zeros(window_length,floor(K/2));
    V3 = zeros(1,floor(K/2));
    for i = 1:N_opt
        comp(:,i) = DiagonalAveraging(V(i,i)*U(:,i).*U(:,i)'...
            + V(end-i+1,end-i+1)*U(:,end-i+1).*U(:,end-i+1)');
        V3(1,i) = 0.5*(V2(i) + V2(end-i+1));
    end
    meanf = meanfreq(comp,Fs);
    bw = powerbw(comp,Fs);
    for i = 1:N_opt
        ind = f >= meanf(i)-0.5*bw(i) & f <= meanf(i)+0.5*bw(i);
        V1(ind,fr) = V1(ind,fr) + 2*V3(1,i)/K;
    end
end

%% Time-frequency distribution plot and Renyi entropy
fsize = 14;
figure('DefaultAxesFontSize',fsize);
RE_MDL = renyi_entropy(V1.^2,3) % Printing Renyi entropy of the TFD
imagesc(time_stamp,f,V1.^2)
colorbar
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% STP-based proposed method
f = (0:1:Fs/2);
V1 = zeros(length(f),Nframes);
time_stamp = zeros(1,Nframes);
for fr = 1:Nframes
    sig = signal((fr-1)*sample_shift+1:(fr-1)*sample_shift+window_length);
    time_stamp(fr) = (1/Fs)*(floor(((fr-1)*sample_shift+1+(fr-1)*...
        sample_shift+window_length)/2)+1);
    K = (window_length + 1)/2;
    X = hankel(sig(1:K),sig(K:end));
    [U,V] = eig(X);
    V2 = abs(diag(V));
    comp = zeros(window_length,floor(K/2));
    V3 = zeros(1,floor(K/2));
    for i = 1:floor(K/2)
        if V2(i) >= (0.1)*max(V2) || V2(end-i+1) >= (0.1)*max(V2)
            comp(:,i) = DiagonalAveraging(V(i,i)*U(:,i).*U(:,i)'...
                + V(end-i+1,end-i+1)*U(:,end-i+1).*U(:,end-i+1)');
            V3(1,i) = 0.5*(V2(i) + V2(end-i+1));
        end
    end
    meanf = meanfreq(comp,Fs);
    bw = powerbw(comp,Fs);
    for i = 1:floor(K/2)
        ind = f >= meanf(i)-0.5*bw(i) & f <= meanf(i)+0.5*bw(i);
        V1(ind,fr) = V1(ind,fr) + 2*V3(1,i)/K;
    end
end

%% Time-frequency distribution plot and Renyi entropy
figure('DefaultAxesFontSize',fsize);
RE_STP = renyi_entropy(V1.^2,3) % Printing Renyi entropy
imagesc(time_stamp,f,V1.^2)
colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% STFT and FSST
%--------------------------------------------------------------------------
% Dowload the STFT and FSST code from the link below
% Link: https://github.com/phamduonghung/FSSTn
%--------------------------------------------------------------------------
gamma = 0;   % Threshold
sigma = 0.12;   % Window parameter
N = 1024;
signal1 = [signal; zeros(624,1)];   % zero padding to make signal length=N
ft = 1:floor(N/2)+1;  % Frequency bin
bt = 1:N;  % Time bin
[stft,fsst,fsst2] = sstn(signal1,gamma,sigma,ft,bt);

%% STFT and FSST - Time-frequency distribution plot and Renyi entropy
figure('DefaultAxesFontSize',fsize);
t2 = (bt-1)/Fs;
f2 = (ft-1)/N*Fs;
TFD_stft = abs(stft(:,1:400)).^2;
RE_stft = renyi_entropy(TFD_stft,3) % Printing Renyi entropy for STFT
imagesc(t2(1:400),f2,TFD_stft);
colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

figure('DefaultAxesFontSize',fsize);
TFD_fsst = abs(fsst(:,1:400)).^2;
RE_fsst = renyi_entropy(TFD_fsst,3) % Printing Renyi entropy for FSST
imagesc(t2(1:400),f2,TFD_fsst);
colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% CWT-based Scalogram
[cwt_coeffs,frq] = cwt(signal,"amor",Fs);
TFD_cwt = abs(cwt_coeffs).^2;
figure('DefaultAxesFontSize',fsize);
RE_cwt = renyi_entropy(TFD_cwt,3)   % Printing Renyi entropy for scalogram
pcolor(t,frq,TFD_cwt);
shading interp
colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% VMD and HSA-based TFD
subbands_vmd = vmd(signal);
[TFD_vmd,f,t] = hht(subbands_vmd,Fs);
TFD_vmd = full(TFD_vmd);
RE_vmd = renyi_entropy(TFD_vmd,3)   % Printing Renyi entropy for VMD-HSA
figure('DefaultAxesFontSize',fsize);
imagesc(t,f,TFD_vmd);
colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')