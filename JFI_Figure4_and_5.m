clc
clearvars
close all

%% Generating synthetic signal
Fs = 1000;  % Sampling rate
t = 0:1/Fs:1;   % Sampling instances
N = length(t);
signal_comps = zeros(N,2);
signal_comps(:,1) = 0.9*sin(2*pi*(160 - (30/2)*t).*t);
signal_comps(:,2) = 0.7*sin(2*pi*(240 - (30/2)*t).*t);
signal =  sum(signal_comps,2);

%% HSA-based reference TFD
fsize = 14;
figure('DefaultAxesFontSize',fsize);
[X,f1,t] = hht(signal_comps,Fs);
RE = renyi_entropy(full(X),3)
imagesc(t,f1,full(X))
cb = colorbar;
xticks([0 0.2 0.4 0.6 0.8 1])
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% Obtaining TFD for various window lengths
W = [71 91 111 131 151];
for k1 = 1:length(W)
window_length = W(k1);
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

%% STP-based proposed method
% f = (0:1:Fs/2);
% V1 = zeros(length(f),Nframes);
% time_stamp = zeros(1,Nframes);
% for fr = 1:Nframes
%     sig = signal((fr-1)*sample_shift+1:(fr-1)*sample_shift+window_length);
%     time_stamp(fr) = (1/Fs)*(floor(((fr-1)*sample_shift+1+(fr-1)*sample_shift+window_length)/2)+1);
%     K = (window_length + 1)/2;
%     X = hankel(sig(1:K),sig(K:end));
%     [U,V] = eig(X);
%     V2 = abs(diag(V));
%     comp = zeros(window_length,floor(K/2));
%     V3 = zeros(1,floor(K/2));
%     for i = 1:floor(K/2)
%         if V2(i) >= (0.1)*max(V2) || V2(end-i+1) >= (0.1)*max(V2)
%             comp(:,i) = DiagonalAveraging(V(i,i)*U(:,i).*U(:,i)'...
%                 + V(end-i+1,end-i+1)*U(:,end-i+1).*U(:,end-i+1)');
%             V3(1,i) = 0.5*(V2(i) + V2(end-i+1));
%         end
%     end
%     meanf = meanfreq(comp,Fs);
%     bw = powerbw(comp,Fs);
%     for i = 1:floor(K/2)
%         ind = f >= meanf(i)-0.5*bw(i) & f <= meanf(i)+0.5*bw(i);
%         V1(ind,fr) = V1(ind,fr) + 2*V3(1,i)/K;
%     end
% end

%% Time-frequency distribution plot and Renyi entropy
figure('DefaultAxesFontSize',fsize);
V1 = (abs(V1)).^2;
RE = renyi_entropy(V1,3)
imagesc(time_stamp,f,V1)
colorbar;
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
end