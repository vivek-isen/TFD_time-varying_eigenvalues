# Representation of time-varying eigenvalues in time-frequency plane
1. We derived the relationship between the amplitude of a discrete-time sinusoidal signal ($x[n] = A\sin(2\pi fn/f_s + \phi)$, for $n = 0,1,\ldots,N-1$) and magnitude of eigenvalue pair ($\lambda_1,\lambda_2$) of the sinusoidal signal obtained from eigenvalue decomposition of Hankel matrix (EVDHM).

$$\frac{|\lambda_1| + |\lambda_2|}{2} \approx \frac{AK}{2}~~~\implies~~~A \approx \frac{|\lambda_1| + |\lambda_2|}{K},~\text{where } K = \frac{N+1}{2}$$

2. A non-stationary signal is segmented into short-duration frames. The EVDHM is performed to each frame and the significant eigenvalue pairs & corresponding components are computed. For $i^\text{th}$ frame, mean frequency ($\text{MNF}_l^{(i)}$), 3 dB bandwidth ($\text{BW}_l^{(i)}$), and amplitude parameter ($A_l^{(i)}$) of $l^\text{th}$ decomposed component are computed. The time-frequency distribution (TFD) of $l^\text{th}$ component of $i^\text{th}$ frame is computed as

$$X_l[i,\zeta] = \left(A_l^{(i)}\right)^2\delta[\zeta - \zeta_l],\forall~\text{MNF}_l^{(i)} - \frac{\text{BW}_l^{(i)}}{2} \leq \zeta_l \leq \text{MNF}_l^{(i)} + \frac{\text{BW}_l^{(i)}}{2}$$

The TFD of all the significant components are added together to get the TFD of the signal. For deeper insight please read [1].

# References
[1]. V.K. Singh and R.B. Pachori, "Eigenvalues-based time-frequency analysis," _Journal of the Franklin Institute_, vol. 363, no. 7, p. 108561, 2026.
