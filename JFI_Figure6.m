clc
clearvars
close all

%% Defining synthetic signal
Fs = 1000;  % Sampling rate
t = 0:1/Fs:1;   % Sampling instances
N = length(t);
signal_comps = zeros(N,2);
signal_comps(:,1) = 0.85*sin(2*pi*(190 + (80/2)*t).*t);
signal_comps(:,2) = 0.6*sin(2*pi*(350 - (180/2)*t).*t);
signal =  sum(signal_comps,2);
N = length(signal);

%% Defining the window length, and overlap length 
window_length = 125;
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
[X,f1,t] = hht(signal_comps,Fs);
imagesc(t,f1,full(X))   % HSA-based reference TFD
colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

figure('DefaultAxesFontSize',fsize);
TFD_mdl = (V1).^2;
RE_mdl = renyi_entropy(TFD_mdl,3)
imagesc(time_stamp,f,TFD_mdl)
colorbar;
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
fsize = 14;
figure('DefaultAxesFontSize',fsize);
TFD_stp = abs(V1).^2;
RE_stp = renyi_entropy(TFD_stp,3)
imagesc(time_stamp,f,TFD_stp)
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
signal1 = [signal; zeros(23,1)];
ft = 1:floor(N/2)+1;  % Frequency bin
bt = 1:N;  % Time bin
[stft,fsst,fsst2] = sstn(signal1,gamma,sigma,ft,bt);
figure('DefaultAxesFontSize',fsize);
t2 = (bt-1)/Fs;
f2 = (ft-1)/N*Fs;
TFD_stft = abs(stft).^2;
RE_stft = renyi_entropy(TFD_stft,3)
imagesc(t2(1:1001),f2,TFD_stft);
cb = colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

figure('DefaultAxesFontSize',fsize);
TFD_fsst = abs(fsst(:,1:1001)).^2;
RE_fsst = renyi_entropy(TFD_fsst,3)
imagesc(t2(1:1001),f2,TFD_fsst);
cb = colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% CWT-based Scalogram
fsize = 14;
[cwt_coeffs,frq] = cwt(signal,"amor",Fs);
figure('DefaultAxesFontSize',fsize);
TFD_cwt = abs(cwt_coeffs).^2;
RE_cwt = renyi_entropy(TFD_cwt,3)
pcolor(t,frq,TFD_cwt);
shading interp
cb = colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% VMD and HSA-based TFD
subbands_vmd = vmd(signal);
[tfr,f,t] = hht(subbands_vmd,Fs);
TFD_vmd = full(tfr);
figure('DefaultAxesFontSize',fsize);

RE_vmd = renyi_entropy(TFD_vmd,3)
imagesc(t,f,TFD_vmd);
cb = colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')