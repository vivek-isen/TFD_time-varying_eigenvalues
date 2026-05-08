clc
clearvars
close all

%% Loading the real-life cello G5 note signal
%--------------------------------------------------------------------------
% Please download the cello G5 note signal from the link provided below.
% Link: https://github.com/dfourer/ASTRES_toolbox (new link - active)
% Old link: https://github.com/dfourer/ASTREStoolbox
%--------------------------------------------------------------------------
[signal,Fs] = audioread("cello_short.wav");
signal = signal(1:2048,1);
N = length(signal);
t = 0:1/Fs:(N-1)/Fs;

%% Defining the window length, and overlap length 
window_length = 231;
overlap_length = window_length - 1;
sample_shift = window_length - overlap_length;
Nframes = floor((N - window_length)/sample_shift) + 1;

%% MDL-based proposed method
f = (0:1:Fs/2);
V1 = zeros(length(f),Nframes);
time_stamp = zeros(1,Nframes);
for fr = 1:Nframes
    sig = signal((fr-1)*sample_shift+1:(fr-1)*sample_shift+window_length);
    time_stamp(fr) = (1/Fs)*(floor(((fr-1)*sample_shift+1+(fr-1)*sample_shift+window_length)/2)+1);
    K = (window_length + 1)/2;
    X = hankel(sig(1:K),sig(K:end));
    [U,V] = eig(X);
    N_opt = MDL_sig_ev(U,V,sig);
    V2 = abs(diag(V));
    comp = zeros(window_length,N_opt);
    V3 = zeros(1,N_opt);
    for i = 1:N_opt
        comp(:,i) = DiagonalAveraging(V(i,i)*U(:,i).*U(:,i)'...
                + V(end-i+1,end-i+1)*U(:,end-i+1).*U(:,end-i+1)');
        V3(1,i) = 0.5*(V2(i) + V2(end-i+1));
    end
    meanf = meanfreq(comp,Fs);
    energy_sig(fr) = sum(sig.^2);
    error_mdl(fr) = sum((sig - sum(comp,2)).^2);
    bw = powerbw(comp,Fs);
    for i = 1:N_opt
        ind = f >= meanf(i)-0.5*bw(i) & f <= meanf(i)+0.5*bw(i);
        V1(ind,fr) = V1(ind,fr) + 2*V3(1,i)/K;
    end
end

%% Framewise ESR plot
fsize = 14;
figure('DefaultAxesFontSize',fsize);
ESR = error_mdl./energy_sig;
plot(ESR,'b','LineWidth',1)
xlim([1 Nframes])
xticks([1 500 1000 1500])
xlabel('Frame number')
ylabel('ESR')

%% STP-based proposed method
f = (0:1:Fs/2);
V1 = zeros(length(f),Nframes);
time_stamp = zeros(1,Nframes);
for fr = 1:Nframes
    sig = signal((fr-1)*sample_shift+1:(fr-1)*sample_shift+window_length);
    time_stamp(fr) = (1/Fs)*(floor(((fr-1)*sample_shift+1+(fr-1)*sample_shift+window_length)/2)+1);
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
    error_stp(fr) = sum((sig - sum(comp,2)).^2);
    bw = powerbw(comp,Fs);
    for i = 1:floor(K/2)
        ind = f >= meanf(i)-0.5*bw(i) & f <= meanf(i)+0.5*bw(i);
        V1(ind,fr) = V1(ind,fr) + 2*V3(1,i)/K;
    end
end

%% Framewise ESR plot
fsize = 14;
figure('DefaultAxesFontSize',fsize);
ESR = error_stp./energy_sig;
plot(ESR,'b','LineWidth',1)
xlim([1 Nframes])
xticks([1 500 1000 1500])
xlabel('Frame number')
ylabel('ESR')