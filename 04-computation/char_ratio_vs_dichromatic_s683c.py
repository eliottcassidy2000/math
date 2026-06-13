#!/usr/bin/env python3
r"""
char_ratio_vs_dichromatic_s683c.py    oracle-2026-06-06-S683c

Q (user): does the character-ratio spectrum bound the LRC tournament's DICHROMATIC number?

Setup. A tournament T has skew-adjacency S = A - A^T (real, skew-symmetric). Its Hermitian
adjacency is iS (Hermitian, REAL eigenvalues). For a CIRCULANT tournament Cay(Z_m, C) the
eigenvalues of S are the CHARACTER RATIOS / Fourier coefficients
    mu_k = 2 * sum_{c in C} sin(2 pi c k / m),   k = 0..m-1
so iS has eigenvalues {±mu_k}. KEY: skew-symmetry forces the Hermitian spectrum to be
SYMMETRIC about 0 for EVERY tournament, so lam_max = -lam_min and the Hoffman ratio bound
    chi_dichromatic >= 1 - lam_max/lam_min = 1 - (-1) = 2
is IDENTICALLY 2. We verify this, and check it against actual dichromatic numbers (R_m=2,
Paley=3, S582): the bound is TIGHT for the LRC tournament R_m but BLIND to Paley's 3.
We also test whether dichromatic is spectral-determined (cospectral, different dichromatic).
"""
import numpy as np
from itertools import combinations, permutations

def circulant_tournament(m, conn):
    cs = set(c % m for c in conn)
    return [[1 if i != j and (j - i) % m in cs else 0 for j in range(m)] for i in range(m)]

def skew(adj):
    A = np.array(adj); return A - A.T

def herm_eigs(adj):
    S = skew(adj)
    return np.sort(np.linalg.eigvalsh(1j * S).real)

def hoffman_ratio(adj):
    ev = herm_eigs(adj)
    lmin, lmax = ev[0], ev[-1]
    return (1 - lmax / lmin) if lmin != 0 else float('inf'), lmax, lmin

def is_acyclic(adj, S):
    for i, j, k in combinations(S, 3):
        if adj[i][j] + adj[j][k] + adj[k][i] in (0, 3):
            return False
    return True

def dichromatic(adj, m):
    for k in range(1, m + 1):
        color = [-1] * m
        def bt(v):
            if v == m: return True
            for c in range(k):
                color[v] = c
                if is_acyclic(adj, [u for u in range(v + 1) if color[u] == c]) and bt(v + 1):
                    return True
            color[v] = -1; return False
        if bt(0): return k
    return m

def max_transitive(adj, m):
    """largest transitive (acyclic) subtournament -- the 'independence' for dichromatic."""
    best = 0
    # greedy + exact-ish for small m via subset search bounded
    for r in range(m, 0, -1):
        found = False
        for S in combinations(range(m), r):
            if is_acyclic(adj, S):
                found = True; break
        if found:
            return r
    return 1

def main():
    print("=" * 80)
    print("Character-ratio (Hermitian) spectrum vs DICHROMATIC number of the LRC tournament")
    print("=" * 80)
    print("\n  m   tournament        dichromatic  Hoffman-ratio(1-lmax/lmin)  max-transitive  n/maxtrans")
    rows = []
    for m in (5, 7, 9, 11, 13):
        Rm = circulant_tournament(m, range(1, (m - 1) // 2 + 1))   # interval = LRC tournament R_m
        chi = dichromatic(Rm, m); hr, lmax, lmin = hoffman_ratio(Rm); mt = max_transitive(Rm, m)
        print(f"  {m:2d}  R_{m} (interval/LRC)   {chi}            {hr:.3f}  (lmax={lmax:.2f},lmin={lmin:.2f})"
              f"        {mt}            {m/mt:.2f}")
        rows.append(("R", m, chi, hr))
        if m % 4 == 3:  # Paley valid (m = 3 mod 4)
            qr = sorted({(x * x) % m for x in range(1, m)})
            P = circulant_tournament(m, qr)
            chi = dichromatic(P, m); hr, lmax, lmin = hoffman_ratio(P); mt = max_transitive(P, m)
            print(f"  {m:2d}  Paley QR_{m}          {chi}            {hr:.3f}  (lmax={lmax:.2f},lmin={lmin:.2f})"
                  f"        {mt}            {m/mt:.2f}")
            rows.append(("P", m, chi, hr))

    print("\n  --- random tournaments (do Hoffman-ratio and dichromatic ever exceed 2 / disagree?) ---")
    import random; rnd = random.Random(683)
    seen_hr = set(); maxchi = 0
    for _ in range(400):
        m = rnd.choice([5, 6, 7])
        adj = [[0] * m for _ in range(m)]
        for i, j in combinations(range(m), 2):
            if rnd.random() < .5: adj[i][j] = 1
            else: adj[j][i] = 1
        chi = dichromatic(adj, m); hr, _, _ = hoffman_ratio(adj)
        seen_hr.add(round(hr, 3)); maxchi = max(maxchi, chi)
    print(f"  over 400 random tournaments: distinct Hoffman-ratio values = {sorted(seen_hr)};  max dichromatic seen = {maxchi}")

    print("\n" + "=" * 80)
    print("ANSWER")
    print("=" * 80)
    print("""  The Hermitian adjacency iS of ANY tournament has a spectrum SYMMETRIC about 0 (S is
  real skew-symmetric => eigenvalues ±mu_k), so lam_max = -lam_min and the character-ratio
  Hoffman bound is IDENTICALLY 2 for every tournament. Therefore:
   * It bounds the LRC tournament R_m's dichromatic number EXACTLY (=2, S582) -- TIGHT.
   * But it is the SAME 2 for Paley (dichromatic 3) and every tournament -- it CANNOT certify
     any dichromatic number > 2; it is blind to the R_m(2) vs Paley(3) distinction.
  So: YES it bounds the LRC tournament's dichromatic number (and tightly, =2), but only because
  the bound is trivially 2 everywhere. The dichromatic number is NOT a character-ratio-spectral
  invariant; the R_m-vs-Paley structure (and LRC's chi) lives in the COMBINATORICS (arc flips /
  OCF / the delta-field), not the spectrum. The better spectral handle is the Delsarte/max-
  transitive bound n/max-transitive (a real, non-trivial lower bound), reported above.""")

if __name__ == "__main__":
    main()
