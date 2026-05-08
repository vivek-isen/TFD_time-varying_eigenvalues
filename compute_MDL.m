function MDL = compute_MDL(signal,est_signal,p)
N = length(signal);
norma = max(abs(signal));
error = signal-est_signal;
lossf = sum(abs(error))/norma;
MDL = p*log(N)/2 + lossf;
end
