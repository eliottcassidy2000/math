#!/usr/bin/env python3
"""
invariant_gaps_klein_S399.py  (klein-2026-07-21-S399, ZOO ATLAS gap-fill V.1)

Compare the two "centrality rankings" of a tournament that the invariant zoo
flagged as never-compared:

  * Kendall-Wei / Perron : RIGHT dominant eigenvector of A (strength propagates
    from beaten opponents; source-heavy).  NB the LEFT eigenvector (A^T) is
    weakness-heavy and gives a spurious near-reversal -- same direction trap as
    THM-1750; use the RIGHT eigenvector.
  * arborescence vector a_r : (r,r) cofactor of the in-Laplacian = out-trees
    rooted at r = stationary distribution of the who-dominates-me walk (THM-1750);
    also source-heavy.

Result: strongly correlated (mean Kendall-tau ~0.94) but NOT identical
(exact-match 174/353 at n=7).  Two distinct-but-aligned centralities; the
divergence set is a new object.  See ZOO ATLAS Part V.1.
"""
import numpy as np
from fractions import Fraction as Fr

# reuse the iso-class harness (BFS + score-refined canon) from the WOWII script
exec(open('04-computation/directed_wowii_klein_S397.py').read().split('# ---- collect invariants')[0])


def perron_rank(om, n):
    A = np.array([[1.0 if beats(om, i, j) else 0.0 for j in range(n)] for i in range(n)])
    w, V = np.linalg.eig(A)            # RIGHT eigenvector = Kendall-Wei STRENGTH (source-heavy)
    k = np.argmax(w.real)
    return tuple(np.abs(V[:, k].real))


def arb_vec(om, n):
    """a_r = (r,r) cofactor of the in-Laplacian D_in - A (out-trees from r)."""
    L = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and beats(om, i, j):
                L[j][j] += 1
                L[i][j] -= 1
    av = []
    for r in range(n):
        idx = [i for i in range(n) if i != r]
        M = [[L[i][j] for j in idx] for i in idx]
        d = Fr(1); sz = len(idx)
        for c in range(sz):
            p = next((x for x in range(c, sz) if M[x][c] != 0), None)
            if p is None:
                d = Fr(0); break
            if p != c:
                M[c], M[p] = M[p], M[c]; d = -d
            d *= M[c][c]; inv = 1 / M[c][c]
            for x in range(c + 1, sz):
                f = M[x][c] * inv
                for y in range(c, sz):
                    M[x][y] -= f * M[c][y]
        av.append(float(d))
    return av


def order(v):
    return tuple(np.argsort(-np.array(v)))


if __name__ == '__main__':
    print("=" * 76)
    print("Kendall-Wei RIGHT-Perron (strength) vs arborescence vector, source-heavy both")
    print("=" * 76)
    for n in (4, 5, 6, 7):
        cls = classes(n); agree = 0; tot = 0; tausum = 0.0; ex = None
        for om in cls:
            if not strong(om, n):
                continue
            tot += 1
            pv = perron_rank(om, n); av = arb_vec(om, n)
            po, ao = order(pv), order(av)
            if po == ao:
                agree += 1
            elif ex is None:
                ex = (po, ao)
            conc = disc = 0
            for i in range(n):
                for j in range(i + 1, n):
                    s1 = pv[i] - pv[j]; s2 = av[i] - av[j]
                    if s1 * s2 > 0:
                        conc += 1
                    elif s1 * s2 < 0:
                        disc += 1
            tausum += (conc - disc) / (conc + disc) if conc + disc else 1
        print(f" n={n}: strong classes {tot}; ranking EXACT-match {agree}/{tot}; "
              f"mean Kendall-tau {tausum/tot:.3f}" + (f"; first mismatch {ex}" if ex else ""))
    print("""
 VERDICT: with both source-heavy, the Perron (spectral) and arborescence
 (combinatorial) centralities are strongly CORRELATED but NOT identical. Two
 distinct-but-aligned centralities; the divergence set is a new object.
""")
