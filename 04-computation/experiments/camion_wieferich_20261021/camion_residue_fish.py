# Fishing test for spark 16 (refined target): is the Camion residue rho_C(p) a fixed small-rational combination
#   rho_C = alpha*q_2 + beta*w_p + gamma*B_(p-3) + delta   (mod p, for all primes)?
# Coefficients alpha, beta, gamma, delta range over rationals n/d with |n| <= 4, d in {1,2,3,4,6,8,12}.  Uses camion_residue_2000.out.
# mac-mini-2026-10-21-fanopsp.
import re, numpy as np
from fractions import Fraction
rows = []
for line in open('camion_residue_2000.out'):
    m = re.match(r'p=(\d+) eps=([+-]1) rho_C=(\d+) q_2=(\d+) q_C=(\d+) w_p=(\d+) B_\(p-3\)=(\d+)', line)
    if m: rows.append(tuple(int(x) for x in m.groups()))
vals = sorted({Fraction(n, d) for n in range(-4, 5) for d in (1, 2, 3, 4, 6, 8, 12)})
print(len(rows), "primes;", len(vals), "candidate values per coefficient;", len(vals) ** 4, "relations")
def res(fr, p): return fr.numerator * pow(fr.denominator, -1, p) % p
cands = None
for (p, eps, rho, q2, qc, wp, b) in rows:
    if p < 50: continue
    v = np.array([res(f, p) for f in vals], dtype=np.int64)
    if cands is None:
        A, B, G, Dl = np.meshgrid(np.arange(len(vals)), np.arange(len(vals)), np.arange(len(vals)), np.arange(len(vals)), indexing='ij')
        A, B, G, Dl = A.ravel(), B.ravel(), G.ravel(), Dl.ravel()
    lhs = (v[A] * q2 + v[B] * wp + v[G] * b + v[Dl]) % p
    keep = lhs == rho % p
    A, B, G, Dl = A[keep], B[keep], G[keep], Dl[keep]
    cands = len(A)
    if cands == 0: print("all relations eliminated by p =", p); break
if cands: print("surviving relations:", [(str(vals[a]), str(vals[b_]), str(vals[g]), str(vals[d])) for a, b_, g, d in zip(A[:10], B[:10], G[:10], Dl[:10])])
# also: rho_C against eps (a constant depending on the class of p mod 8)?
for e in (1, -1):
    sub = [r for r in rows if r[1] == e]
    print(f"eps={e:+d}: rho_C = 0 for {sum(1 for r in sub if r[2] == 0)} of {len(sub)} primes; rho_C = 1 for {sum(1 for r in sub if r[2] == 1)}")
