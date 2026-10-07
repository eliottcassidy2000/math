#!/usr/bin/env python3
"""Merit factors / sup norms of the rows of the bordered Paley skew-Hadamard matrix of order p+1 (p = 3 mod 4),
rows in the cyclic (Z_p) vertex order: row i = (border, chi(j - i) with chi(0) := +1).
Høholdt-Jensen: rotation fraction r gives 1/F -> 2/3 - 4|r| + 8 r^2 (|r| <= 1/2), max 6 at r = 1/4.
Also: Barker 7 vs Legendre 7, and the tower's H_8 vs the Paley H_8 (both are skew Hadamard of order 8)."""
import numpy as np, sympy
from tower_merit import autocorr_int, supnorm, tower

def legendre(p):
    chi = np.zeros(p, dtype=np.int8)
    qr = set((x * x) % p for x in range(1, p))
    for j in range(1, p):
        chi[j] = 1 if j in qr else -1
    chi[0] = 1
    return chi

def paley_bordered(p):
    chi = legendre(p)
    N = p + 1
    H = np.ones((N, N), dtype=np.int8)
    H[1:, 0] = -1
    for i in range(p):
        H[i + 1, 1:] = np.roll(chi, i)  # entry j: chi(j - i)
    return H

def F_of(X):
    C = autocorr_int(X)
    N = X.shape[1]
    return N * N / (2.0 * (C[:, 1:] ** 2).sum(axis=1))

if __name__ == "__main__":
    b7 = np.array([[1, 1, 1, -1, -1, 1, -1]], dtype=np.int8)
    print("Barker 7 F =", F_of(b7)[0], " C_u:", autocorr_int(b7)[0][1:])
    l7 = legendre(7)
    rots = np.array([np.roll(l7, -s) for s in range(7)], dtype=np.int8)
    revs = rots[:, ::-1]
    allv = np.concatenate([rots, revs, -rots, -revs])
    print("Barker 7 is a (reversed/negated) rotation of Legendre-7 (chi(0)=+1):",
          any(np.array_equal(v, b7[0]) for v in allv))
    for p in [7, 11, 19, 23, 31, 43, 127, 251, 499, 1019, 2039, 4091, 8191]:
        if not sympy.isprime(p) or p % 4 != 3:
            continue
        H = paley_bordered(p); N = p + 1
        G = H.astype(np.int64) @ H.T.astype(np.int64)
        assert np.array_equal(G, N * np.eye(N, dtype=np.int64)), p
        assert np.array_equal(H.astype(int) + H.T.astype(int), 2 * np.eye(N, dtype=int)), p
        F = F_of(H[1:])
        # also plain Legendre rotations (length p) best
        chi = legendre(p)
        R = np.array([np.roll(chi, -s) for s in range(p)], dtype=np.int8)
        FL = F_of(R)
        ub, mx, mn = supnorm(H[1:], 8)
        print(f"p={p:5d} N={N:5d} bordered-Paley rows: Fmax {F.max():.4f} at row {int(F.argmax())+1} (r={((int(F.argmax()))/p):.3f}), Fmean {F.mean():.4f}, "
              f"best sup/sqrtN {(mx/np.sqrt(N)).min():.3f} | Legendre rotations Fmax {FL.max():.4f} at shift {int(FL.argmax())} (r={int(FL.argmax())/p:.3f})")
    H8 = tower(3)
    print("tower H_8 rows F:", np.round(F_of(H8), 4))
    print("Paley-7 bordered H_8 rows F:", np.round(F_of(paley_bordered(7)), 4))
