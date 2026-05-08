clc
clearvars
close all

%% Generating signal
Fs = 1000;
t = 0:1/Fs:1;
N = length(t);
window_length = 125;
K = (window_length + 1)/2;
overlap_length = window_length - 1;
sample_shift = window_length - overlap_length;
Nframes = floor((N - window_length)/sample_shift) + 1;
A = -4:0.1:0;

num_comps = zeros(length(A),1);
for i = 1:length(A)
    signal = sin(2*pi*40*t) + 10^(A(i))*sin(2*pi*30*t);
    for fr = 17
        sig = signal((fr-1)*sample_shift+1:(fr-1)*sample_shift+window_length);
        X = hankel(sig(1:K),sig(K:end));
        [U,V] = eig(X);
        V2 = abs(diag(V));
        for j = 1:floor(K/2)
            if V2(j) >= (0.1)*max(V2) || V2(end-j+1) >= (0.1)*max(V2)
                num_comps(i) = num_comps(i) + 1;
            end
        end
    end
end

fsize = 14;
figure('DefaultAxesFontSize',fsize);
semilogx(10.^(A),num_comps,'b','LineWidth',1)
xlim([10^(A(1)),10^(A(end))])
xticks([10^(-4), 10^(-3), 10^(-2), 10^(-1), 10^(0)])
xlabel('$A$', Interpreter='latex')
ylabel('$N_{\textrm{pair}}$', Interpreter='latex')

num_comps = zeros(length(A),1);
for i = 1:length(A)
    signal = sin(2*pi*40*t) + 10^(A(i))*sin(2*pi*30*t);
    for fr = 17
        sig = signal((fr-1)*sample_shift+1:(fr-1)*sample_shift+window_length);
        X = hankel(sig(1:K),sig(K:end));
        [U,V] = eig(X);
        V2 = abs(diag(V));
        for j = 1:floor(K/2)
            if V2(j) >= (0.01)*max(V2) || V2(end-j+1) >= (0.01)*max(V2)
                num_comps(i) = num_comps(i) + 1;
            end
        end
    end
end

figure('DefaultAxesFontSize',fsize);
semilogx(10.^(A),num_comps,'b','LineWidth',1)
xlim([10^(A(1)),10^(A(end))])
xticks([10^(-4), 10^(-3), 10^(-2), 10^(-1), 10^(0)])
xlabel('$A$', Interpreter='latex')
ylabel('$N_{\textrm{pair}}$', Interpreter='latex')


%% Plot of MDL-based threshold

% num_comps = zeros(length(A),1);
% for i = 1:length(A)
%     signal = sin(2*pi*40*t) + exp(A(i))*sin(2*pi*30*t);
%     signal = signal';
%     X = hankel(signal(1:501),signal(501:end));
%     [U,V] = eig(X);
%     num_comps(i) = MDL_sig_ev(U,V,signal);
% end
%
% subplot(1,2,2)
% semilogx(exp(A),num_comps)