function N_opt = MDL_sig_ev(U,V,signal)
V = diag(V);
N = length(signal);
N1 = floor(length(V)/2);
comps = zeros(N,N1);
MDL = zeros(N1,1);
for i = 1:N1
    X1 = V(i)*U(:,i)*U(:,i)' + V(end-i+1)*U(:,end-i+1)*U(:,end-i+1)';
    comps(:,i) = DiagonalAveraging(X1);
    MDL(i,:) = compute_MDL(signal,sum(comps,2),i);
end
[~,N_opt] = min(MDL);
end