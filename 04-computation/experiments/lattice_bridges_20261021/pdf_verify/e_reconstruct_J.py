#!/usr/bin/env python3
"""Try to reconstruct paper 4's Hermitian factor J: O_K (K = Q(zeta35)) with h_J(x,y) = Tr_{K/F}(lambda x ybar),
lambda = (u * zeta^30 * delta / sqrt(-7))^{-1}, delta = P'(zeta) (P = min poly of zeta over F = Q(sqrt(-7))),
reading the undefined 'u = -u4 u8 u16' as u_k = zeta^k + zeta^-k (real cyclotomic units).
Checks: reality and total positivity of lambda, integrality of Tr_{K/Q}(lambda x ybar) on Z[zeta], det, divisibility by 7,
and writes the 24x24 Gram (divided by 7 if possible) for PARI."""
import itertools, sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nint, fabs
mp.dps = 80
n = 35
units = [k for k in range(1, n) if __import__('math').gcd(k, n) == 1]
H = [k for k in units if pow(k, 3, 7) == 1]              # Gal(K/F)
def leg7(k): return 1 if pow(k % 7, 3, 7) == 1 else -1
z = lambda k: exp(2j*pi*mpf(k)/n)
def sig_delta(k):
    r = mpc(1)
    for h in H:
        if h != 1: r *= (z(k) - z(k*h))
    return r
def sig_s7(k): return leg7(k) * 1j*sqrt(7)
def sig_u(k, which):
    r = mpc(1)
    for a in which: r *= (z(k*a) + z(-k*a))
    return r
def lam(k, sign=-1, which=(4, 8, 16)):
    return 1/(sign*sig_u(k, which) * z(30*k) * sig_delta(k) / sig_s7(k))
# reality: sigma_{-k}(lambda) = conj(sigma_k(lambda))
real_ok = all(fabs((lam(k) - lam(n-k).conjugate())) < mpf(10)**-40 for k in units)
imag_ok = all(fabs(lam(k).imag) < mpf(10)**-40 for k in units)
print("lambda in K^+ (sigma_k(lambda) real for all k):", imag_ok)
signs = {k: (1 if lam(k).real > 0 else -1) for k in units}
print("signs of sigma_k(lambda) for k in reps 1..17:", [signs[k] for k in units if k < n/2])
print("totally positive:", all(s > 0 for s in signs.values()))
def gram(lamf):
    G = [[None]*24 for _ in range(24)]
    for a in range(24):
        for b in range(24):
            s = sum(lamf(k) * z(k*(a-b)) for k in units)
            G[a][b] = s
    return G
if all(s > 0 for s in signs.values()):
    G = gram(lam)
    maxerr = max(fabs(G[a][b] - nint(G[a][b].real)) for a in range(24) for b in range(24))
    Gi = [[int(nint(G[a][b].real)) for b in range(24)] for a in range(24)]
    print("max distance from integers:", mp.nstr(maxerr, 3))
    div7 = all(x % 7 == 0 for r in Gi for x in r)
    print("all entries divisible by 7:", div7)
    import sympy as sp
    d = sp.Matrix(Gi).det()
    print("det =", d, "= 7^24:", d == 7**24)
    with open('J_gram_div7.gp', 'w') as f:
        M = [[x//7 for x in r] for r in Gi] if div7 else Gi
        f.write("GJ = [" + ";".join(",".join(str(x) for x in r) for r in M) + "];\n")
    print("wrote J_gram_div7.gp (divided by 7)" if div7 else "wrote J_gram_div7.gp (NOT divided)")
