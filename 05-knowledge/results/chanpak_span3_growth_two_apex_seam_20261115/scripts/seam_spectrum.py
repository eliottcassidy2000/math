#!/usr/bin/env python3
"""seam_spectrum.py: the value set of the signed-seam family H = (2^k - b)(2^r + a) 2^t + C(k, r, a, b) (PRODUCT.md (*)),
by forward enumeration over k, r <= KMAX, t <= TMAX, all masks; reports the density of odd N covered below NMAX."""
import sys, numpy as np
KMAX, TMAX, NMAX = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
def rev(z, w): return int(format(z, f'0{w}b')[::-1], 2) if w > 0 else 0
def J(z):
    bits = [i for i in range(z.bit_length()) if z >> i & 1]; s = 0
    for x in range(len(bits)):
        for y in range(x):
            i, j = bits[x], bits[y]; s += 2 ** (i - j - 1) + (1 if i == j + 1 else 0)
    return s
cov = np.zeros(NMAX, dtype=np.uint8); first = {}
Jt = {}
def Jc(z):
    if z not in Jt: Jt[z] = J(z)
    return Jt[z]
for k in range(2, KMAX + 1):
    for r in range(2, KMAX + 1):
        bs = np.arange(2 ** (k - 1)); as_ = np.arange(2 ** (r - 1))
        revb = np.array([rev(int(b), k) for b in bs]); reva = np.array([rev(int(a), r) for a in as_])
        Jb = np.array([Jc(int(b)) for b in bs]); Ja = np.array([Jc(int(a)) for a in as_])
        D = (2 ** k - bs)[:, None]; E = (2 ** r + as_)[None, :]
        P = (D * E).astype(np.int64); C = (1 + revb[:, None] - reva[None, :] - Jb[:, None] - Ja[None, :]).astype(np.int64)
        for t in range(1, TMAX + 1):
            Nv = (P << t) + C
            sel = Nv[(Nv > 0) & (Nv < NMAX)]
            cov[sel] = 1
print(f"k, r <= {KMAX}, t <= {TMAX}: odd N < {NMAX} covered: {int(cov[1::2].sum())} of {NMAX // 2}")
for lo, hi in ((1, 2**8), (2**8, 2**10), (2**10, 2**12), (2**12, 2**14), (2**14, 2**16), (2**16, 2**18), (2**18, 2**20)):
    if hi > NMAX: break
    odd = np.arange(lo | 1, hi, 2); c = cov[odd].sum()
    print(f"  [{lo}, {hi}): covered {c} of {len(odd)} ({100 * c / len(odd):.1f}%)")
miss = [n for n in range(1, min(NMAX, 400), 2) if not cov[n]]
print("misses below 400:", miss)
