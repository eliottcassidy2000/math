"""Exact spectrum of the untruncated level operator T_n = (S/2)(I - S/2)^-1 D on C^L, L = 2 3^(n-1), D = diag(e(2^k/3^n)):
claim: eigenvalues are the roots of (2^L - 1) lambda^L + lambda^(L/2) - 1 = 0, i.e. lambda^(L/2) = mu_+- with
(2^L - 1) mu^2 + mu - 1 = 0; eigenvectors w(k-1) = 2 lambda/(omega_k + lambda) w(k)."""
import numpy as np, sys
for n in range(2, 7):
    L = 2*3**(n-1); mod = 3**n
    om = np.exp(2j*np.pi*np.array([pow(2,k,mod) for k in range(L)])/mod)
    S = np.roll(np.eye(L), 1, axis=0)              # (S v)(k) = v(k-1)
    C = (S/2) @ np.linalg.inv(np.eye(L) - S/2)     # untruncated geometric circulant
    T = C * om[None, :]
    ev = np.sort_complex(np.linalg.eigvals(T))
    # predicted
    mu = np.roots([2**L - 1, 1, -1])
    pred = np.concatenate([np.array([m])**(2/L) * np.exp(2j*np.pi*np.arange(L//2)/(L//2)) for m in mu])  # principal root times (L/2)-th roots of unity
    # compare multisets by sorting moduli/args
    a = np.sort(np.abs(ev)); b = np.sort(np.abs(pred))
    # also check the polynomial vanishes at the numerical eigenvalues
    P = lambda lam: (2**L-1)*lam**L + lam**(L/2) - 1
    resid = max(abs(P(l))/(2**L) for l in ev)   # scaled residual
    # eigenvector check for one eigenvalue: w(k-1) = 2 lam/(om_k+lam) w(k) -> build w and test T? (w is the auxiliary vector; v = (S/2) w / lam)
    lam = ev[0]; w = np.ones(L, dtype=complex)
    for k in range(L-1, 0, -1): w[k-1] = 2*lam/(om[k]+lam)*w[k]
    v = (S @ w)/(2*lam); res_vec = np.linalg.norm(T @ v - lam*v)/np.linalg.norm(v)
    print(f"n={n} L={L}: |eig| in [{a.min():.6f},{a.max():.6f}] vs predicted moduli {b.min():.6f},{b.max():.6f}; max|eig - nearest predicted| = {max(min(abs(e-p) for p in pred) for e in ev):.2e}; scaled residual of the char. polynomial {resid:.1e}; eigenvector check {res_vec:.1e}; mu = {mu}")
