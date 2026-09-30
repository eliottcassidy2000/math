"""The Collatz frequency recursion as a weighted circulant digraph on Z/L_n (L_n = 2 3^(n-1)) composed with the
Gauss-sum diagonal omega_n(k) = e(2^k/3^n):  (T_n v)(k) = sum_{a=1}^{A} 2^-a omega_n(k-a) v(k-a).
(1) autocorrelation of omega_n: R(d) = sum_k omega(k) conj(omega(k+d)) = L (d=0), -L/2 (d = +-L/3), 0 otherwise (n >= 2);
(2) mean square of the one-step Gauss sums over the units: (1/L) sum_k |G_n(2^k)|^2 = sum_a 4^-a = 1/3 (A < 2 3^(n-2));
(3) the geometric circulant's spectrum lies on the circle |w - 1/3| = 2/3;
(4) eigenvalues of T_n (n <= NMAX): spectral radius, |det|^(1/L) (= 1/2 by Jensen for the untruncated kernel), numerical range."""
import numpy as np, math, sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
A = 40
for n in range(2, NMAX+1):
    L = 2*3**(n-1); mod = 3**n
    r = np.array([pow(2, k, mod) for k in range(L)], dtype=np.float64)
    om = np.exp(2j*np.pi*r/mod)
    # (1) autocorrelation via FFT
    F = np.fft.fft(om); R = np.fft.ifft(np.abs(F)**2).real
    R0 = R[0]; R3 = R[L//3]; Rmax_other = max(abs(R[d]) for d in range(1, L) if d not in (L//3, 2*L//3))
    # spectrum flatness
    prim = [xi for xi in range(L) if xi % 3]; imp = [xi for xi in range(L) if xi % 3 == 0]
    fl = np.abs(F[prim]) / 3**(n/2); imp_max = np.abs(F[imp]).max()
    # (2) one-step Gauss sums G(2^k) = sum_a 2^-a omega(k - a)
    w = 2.0**(-np.arange(1, A+1))
    G = np.zeros(L, dtype=complex)
    for a in range(1, A+1): G += w[a-1] * np.roll(om, a)      # G[k] = sum_a w_a om[k-a]
    ms = np.mean(np.abs(G)**2)
    # (4) T_n = C D with C circulant (kernel w at shifts 1..A), D = diag(om)
    C = np.zeros((L, L)); 
    for a in range(1, A+1): C += w[a-1] * np.roll(np.eye(L), a, axis=0)   # (C v)(k) = sum_a w_a v(k-a)
    T = C * om[None, :]          # (T v)(k) = sum_j C[k,j] om[j] v[j]
    ev = np.linalg.eigvals(T)
    rho = np.abs(ev).max(); gm = np.exp(np.mean(np.log(np.abs(ev))))
    # numerical radius (approx): max |<Tv,v>| over random unit v and over the eigen/singular directions
    U, s, Vh = np.linalg.svd(T)
    print(f"n={n} L={L}: R(0)={R0:.3f}=L? {abs(R0-L)<1e-6}; R(L/3)={R3:.3f} (=-L/2={-L/2}); max|R| elsewhere={Rmax_other:.2e}; |omega_hat|/3^(n/2) on primitive xi in [{fl.min():.6f},{fl.max():.6f}], on 3|xi max {imp_max:.2e}")
    print(f"      mean |G_n(2^k)|^2 over the cycle = {ms:.9f} (1/3 = {1/3:.9f}, sum_a 4^-a to A = {sum(4.0**-a for a in range(1,A+1)):.9f});  spectral radius rho(T_n) = {rho:.5f}, geometric mean |eig| = {gm:.5f} (Jensen 1/2), sigma_max = {s[0]:.5f}, sigma_min = {s[-1]:.5f}")
print("geometric circulant spectrum: w(theta) = 1/(2 e^{i theta} - 1): check |w - 1/3| = 2/3 at theta = 0.7:", abs(1/(2*np.exp(0.7j)-1) - 1/3))
