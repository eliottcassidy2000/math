#!/usr/bin/env python3
"""The pairwise sector co-emptiness matrix M (the S2 carrier) -- what is its structure?

HYP-3084: the gK8/Delsarte concentration bound is led by S2 = sum_{i<j} meas{sectors i,j both empty}.
This script forms the 7x7 matrix  M[i][j] = meas{tau : sectors i AND j both empty}  at the binding
config (consec_k) and asks: is it CIRCULANT on Z/7 (the apex-7 / Gauss-sum / sqrt(-7) face), and what
is its spectrum?  M[i][i] = meas{sector i empty} (the S1 marginals); off-diagonal = the covariance.

Two candidate structures to test:
  (a) Clebsch/biplane:  N^T N = 4I + 2J  (the 2-(16,6,2) design Gram, the H4 / cut-space face)
  (b) Gauss-sum circulant on Z/7:  eigenvalues built from the QR/Fano characters (the H2 / apex-7 face)
A circulant on Z/7 is diagonalized by the additive characters; its eigenvalues are the DFT of one row,
and the QR-symmetric part yields the Gauss sum sqrt(-7).  Knowing which structure M has tells us which
face of 14=2.7 certifies the bound.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def sector_of(p): return int((p % 1) * 7)

def pairwise_matrix(E):
    """M[i][j] = meas{tau: sectors i,j both empty}, i,j in 0..6.  Exact rationals."""
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b)
    M = [[F(0)] * 7 for _ in range(7)]
    for t in range(len(b) - 1):
        x0, x1 = b[t], b[t + 1]
        if x1 <= x0: continue
        w = x1 - x0
        covered = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        empty = [s for s in range(7) if s not in covered]
        for i in empty:
            for j in empty:
                M[i][j] += w
    return M

def is_circulant(M, mod=7, restrict=None):
    """Check M[i][j] depends only on (i-j) mod 7.  restrict: subset of indices to test."""
    idx = restrict if restrict is not None else range(7)
    by_diff = {}
    ok = True
    for i in idx:
        for j in idx:
            d = (i - j) % mod
            if d in by_diff and by_diff[d] != M[i][j]: ok = False
            by_diff.setdefault(d, M[i][j])
    return ok, by_diff

def eigvals_float(M):
    import numpy as np
    A = np.array([[float(x) for x in row] for row in M])
    return sorted(np.linalg.eigvalsh((A + A.T) / 2).tolist())

print("=" * 78)
print(" PAIRWISE SECTOR CO-EMPTINESS MATRIX M -- structure & spectrum at the binding config")
print("=" * 78)
for k in (8, 9, 10):
    consk = tuple(range(k))
    M = pairwise_matrix(consk)
    # sector 0 always hit (0 in E) -> row/col 0 should be ~0
    row0 = [float(M[0][j]) for j in range(7)]
    inner = list(range(1, 7))
    circ_inner, diffs = is_circulant(M, restrict=inner)
    circ_all, diffs_all = is_circulant(M, restrict=range(7))
    S1 = sum(M[i][i] for i in range(7))           # = E[N], the 1st moment
    S2 = sum(M[i][j] for i in range(7) for j in range(7) if i < j)  # the pairwise total
    print(f"\n--- k={k}  consec_{k} ---")
    print(f"  sector-0 row (should be ~0, pinned by e=0): max|M[0][j]| = {max(abs(x) for x in row0):.2e}")
    print(f"  diagonal (S1 marginals, meas{{sector i empty}}): "
          f"[{', '.join(f'{float(M[i][i]):.4f}' for i in range(7))}]")
    print(f"  S1 = sum diag = {float(S1):.4f}   S2 = sum_{{i<j}} = {float(S2):.4f}")
    print(f"  CIRCULANT on inner 6 sectors (M[i][j]=f(i-j mod 7))? {circ_inner}")
    if circ_inner:
        # show the circulant generator over the 6 inner sectors (differences 1..6, plus self)
        gen = {d: float(v) for d, v in sorted(diffs.items())}
        print(f"     circulant row by (i-j) mod 7: {{ {', '.join(f'{d}:{v:.4f}' for d,v in gen.items())} }}")
        # QR vs NQR symmetry: QR(7)={1,2,4}, NQR={3,5,6}.  Gauss-sum face <=> M[d] same on each class
        QR, NQR = {1, 2, 4}, {3, 5, 6}
        qr_vals = set(round(diffs.get(d, F(0)), 12) for d in QR if d in diffs)
        nqr_vals = set(round(diffs.get(d, F(0)), 12) for d in NQR if d in diffs)
        print(f"     QR{{1,2,4}}-class off-diag values: {[float(v) for v in qr_vals]}")
        print(f"     NQR{{3,5,6}}-class off-diag values: {[float(v) for v in nqr_vals]}")
        print(f"     -> QR/NQR split (Gauss-sum/Fano face)? "
              f"{'YES (QR!=NQR, the sqrt(-7) signature)' if qr_vals != nqr_vals else 'NO (QR==NQR, plain circulant)'}")
    ev = eigvals_float(M)
    print(f"  spectrum of M (symmetrized): [{', '.join(f'{v:.4f}' for v in ev)}]")
    # rank / multiplicities hint at design (4I+2J -> eigenvalues {4+2*6, 4} = {16,4} on 16 pts; here 7x7)
print("\nINTERPRETATION KEY:")
print(" - circulant on inner 6 + QR/NQR split  => the apex-7 / Gauss-sum / Fano (H2) face certifies S2")
print(" - eigenvalues of shape {a, b (mult m)}  => a design/SRG (Clebsch 4I+2J, H4) face certifies S2")
print(" - the consec extremum of S2 is then read from the top eigenvector / the QR-class sign")
