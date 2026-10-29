#!/usr/bin/env python3
"""Analyse the P19 two-run rows from a13_p19.c (path given as argv[1]): B_0, det B_1 in Z[zeta_19] (CRT),
invariance under sigma_4, ninth-power law, a + b*eta, content, quadratic norm (trial division to 10^7)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a3_runobs import bareiss_det, detB1_coeffs, sigma, norm_Qzeta, eta_coords, quad_norm, split_valuations, kron, partial_factor
from math import gcd
q = 19
rows = {}
for line in open(sys.argv[1]):
    s2, t1, t2, val = map(int, line.split())
    rows[(s2, t1, t2)] = val
R = [(0, s2) for s2 in range(1, q)]
def M(r, t):  # r = (0, s2)
    return rows.get((r[1], t[0], t[1]), 0)
Ms = [[[M(r, ((rp[0] + s) % q, (rp[1] + s) % q)) for s in range(q)] for rp in R] for r in R]
B0 = [[sum(x) for x in row] for row in Ms]
dB0 = bareiss_det(B0)
print("det B_0 =", dB0, "=", partial_factor(dB0))
c, used = detB1_coeffs(Ms, q)
QR = sorted({(x * x) % q for x in range(1, q)})
print("CRT primes used:", used, "; fixed by sigma_4:", sigma(c, 4, q) == c, "; fixed by all of QR:", all(sigma(c, m, q) == c for m in QR),
      "; moved by sigma_2:", sigma(c, 2, q) != c)
a, b = eta_coords(c, q, set(QR))
nq = quad_norm(a, b, q)
NB = norm_Qzeta(c, q)
print("N(det B_1) == quadratic norm^9:", NB == nq ** 9, "; quadratic norm digits:", len(str(nq)))
print("content gcd(a,b) =", partial_factor(gcd(a, b)))
pf = partial_factor(nq)
print("quadratic norm partial factorisation:", pf)
for p, e in sorted((p, e) for p, e in pf.items() if isinstance(p, int)):
    kr = kron(-q, p)
    if kr == 1:
        print(f"   p={p}: split, exponent {e}, (v_P, v_Pbar) = {split_valuations(a, b, p, q, e)}")
    else:
        print(f"   p={p}: {'inert' if kr == -1 else 'ramified'}, exponent {e}")
