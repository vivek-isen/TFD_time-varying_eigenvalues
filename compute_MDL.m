function MDL = compute_MDL(signal,est_signal,p)
N = length(signal);
norma = max(abs(signal));
% error = signal - est_signal;
% mean squared error loss
error = signal-est_signal;
lossf = sum(abs(error))/norma;
% log-cosh loss
% lossf = sum(log(cosh(error/norma)));
MDL = p*log(N)/2 + lossf;
end