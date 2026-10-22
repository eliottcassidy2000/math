#!/usr/bin/env python3
"""Build paper 4's unglued lattice L_raw = E8 (x)_{O_F} J for reconstructed J's (totally positive twists found by
e_search_J.py), as a 96x96 integer Gram matrix (trace form), and look for short vectors:
 (1) export for PARI LLL; (2) targeted rank-2 tensors x1(x)y1 + x2(x)y2 with x_i E8-roots, y_i minimal vectors of J."""
import itertools, math, sys, json
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nint, fabs
import numpy as np
mp.dps = 60
n = 35
units = [k for k in range(1, n) if math.gcd(k, n) == 1]
reps = [k for k in units if k < n/2]
Hgal = [k for k in units if pow(k, 3, 7) == 1]
def leg7(k): return 1 if pow(k % 7, 3, 7) == 1 else -1
z = lambda k: exp(2j*pi*mpf(k)/n)
def sig_delta(k):
    r = mpc(1)
    for h in Hgal:
        if h != 1: r *= (z(k) - z(k*h))
    return r
lam0 = {k: 1/(z(30*k) * sig_delta(k) / (leg7(k)*1j*sqrt(7))) for k in units}
gens = {('r', a): {k: (z(k*a) + z(-k*a)).real for k in units} for a in range(1, 18)}
gens[('m1',)] = {k: mpf(-1) for k in units}
names = list(gens)
def signvec(vals): return np.array([0 if vals[k] > 0 else 1 for k in reps], dtype=np.int64)
S = np.array([signvec(gens[g]) for g in names]); t = signvec({k: lam0[k].real for k in units})
def rref_mod2(A):
    A = A.copy() % 2; rows, cols = A.shape; piv = []; r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i, c]), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]
        for i in range(rows):
            if i != r and A[i, c]: A[i] ^= A[r]
        piv.append(c); r += 1
    return A, piv
R, piv = rref_mod2(np.concatenate([S.T, t.reshape(-1, 1)], axis=1))
x0 = np.zeros(len(names), dtype=np.int64)
for i, c in enumerate(piv): x0[c] = R[i, -1]
free = [c for c in range(len(names)) if c not in piv]
kernel = []
for f in free:
    v = np.zeros(len(names), dtype=np.int64); v[f] = 1
    for i, c in enumerate(piv): v[c] = R[i, f]
    kernel.append(v)
which = [int(a) for a in sys.argv[1:]] or [0]
s7 = sqrt(7)
# O_F-basis of O_K: zeta^(7a+5b)
E = [7*a + 5*b for a in range(4) for b in range(3)]
def to_F(v):
    # v = (p + q sqrt(-7))/2 ; return (p, q) as Fractions with small denominators
    p = Fr(int(nint(2*v.real*7)), 7) if fabs(2*v.real*7 - nint(2*v.real*7)) < 1e-25 else None
    q = Fr(int(nint(2*v.imag/s7*7)), 7) if fabs(2*v.imag/s7*7 - nint(2*v.imag/s7*7)) < 1e-25 else None
    return p, q
# H_E from paper 4 as (p,q) pairs: value = (p + q sqrt(-7))/2
# H_E = [[2,3, -3i/sqrt7, -2i/sqrt7],[3,5,-5i/sqrt7,-3i/sqrt7],[3i/sqrt7,5i/sqrt7,1,0],[2i/sqrt7,3i/sqrt7,0,1]]
# i/sqrt7 = sqrt(-7)/7 -> (p,q) = (0, 2/7)
def F(p, q): return (Fr(p), Fr(q))
HE = [[F(4,0), F(6,0), F(0,Fr(-6,7)), F(0,Fr(-4,7))],
      [F(6,0), F(10,0), F(0,Fr(-10,7)), F(0,Fr(-6,7))],
      [F(0,Fr(6,7)), F(0,Fr(10,7)), F(2,0), F(0,0)],
      [F(0,Fr(4,7)), F(0,Fr(6,7)), F(0,0), F(2,0)]]
def fmul(a, b):   # ((p1+q1 s)/2)((p2+q2 s)/2), s^2=-7  -> (p,q) of result
    p1, q1 = a; p2, q2 = b
    return ((p1*p2 - 7*q1*q2)/2, (p1*q2 + q1*p2)/2)
def fre(a): return a[0]/2           # real part
w = (Fr(1), Fr(1))                  # omega = (1 + sqrt(-7))/2
wb = (Fr(1), Fr(-1))
for ci in which:
    x = x0.copy()
    bits = [(ci >> j) & 1 for j in range(len(kernel))]
    for b, v in zip(bits, kernel): x = (x + b*v) % 2
    lam = {}
    for k in units:
        val = lam0[k].real
        for g, xi in zip(names, x):
            if xi: val *= gens[g][k]
        lam[k] = val
    HJ = [[None]*12 for _ in range(12)]
    for i in range(12):
        for j in range(12):
            v = sum(lam[k]*z(k*(E[i]-E[j])) for k in Hgal)
            HJ[i][j] = to_F(v)
            assert None not in HJ[i][j], (i, j, v)
    integral = all(((p + q) % 2 == 0) and p.denominator == 1 and q.denominator == 1 for row in HJ for (p, q) in row)
    print(f"candidate {ci}: H_J entries in O_F: {integral}; diag h_J(b_i,b_i) = {[int(HJ[i][i][0]/2) for i in range(12)]}")
    # Kronecker product and trace form
    Hraw = [[fmul(HE[i//12][j//12], HJ[i % 12][j % 12]) for j in range(48)] for i in range(48)]
    G = np.zeros((96, 96), dtype=object)
    for f_ in range(48):
        for g_ in range(48):
            h = Hraw[f_][g_]
            G[f_, g_] = 2*fre(h)
            G[48+f_, g_] = 2*fre(fmul(w, h))
            G[f_, 48+g_] = 2*fre(fmul(wb, h))
            G[48+f_, 48+g_] = 2*fre(fmul(fmul(w, wb), h))
    assert all(x.denominator == 1 for x in G.flatten())
    Gi = np.array([[int(x) for x in row] for row in G], dtype=object)
    print("  G_raw symmetric:", (Gi == Gi.T).all(), " even:", all(Gi[i, i] % 2 == 0 for i in range(96)), " min diag:", min(Gi[i, i] for i in range(96)))
    with open(f'Lraw_cand{ci}.gp', 'w') as fh:
        fh.write("G = [" + ";".join(",".join(str(int(v)) for v in row) for row in Gi) + "];\n")
    with open(f'HJ_cand{ci}.json', 'w') as fh:
        json.dump([[[str(p), str(q)] for (p, q) in row] for row in HJ], fh)
print("done")
