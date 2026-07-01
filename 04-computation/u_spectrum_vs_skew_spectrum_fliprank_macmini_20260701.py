#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S89 -- does the U-spectrum see the flip-rank excess the skew-spectrum misses?

DEFINITIONS (the natural reading; owner introduced 'U-spectrum'):
 - U-spectrum   = eigenvalues of the 0/1 adjacency A (the 'unsigned'/directed tournament matrix) = kps's cpA.
 - skew-spectrum= eigenvalues of S = A - A^T (the +/-1 skew-adjacency; purely imaginary; converse-EVEN).
 - flip-rank excess = the S_n-folding excess in the covering number k(n)=kappa (0,0,0,1,3 for n=3..7, opus/klein),
   DRIVEN by high-|Aut| classes (the few-rep 'needles' a thin subcube can't cover). Per-class proxies:
   |Aut| (symmetry) and MFAS (min feedback arc set = covering-radius contribution).

Note A = (J - I + S)/2, so A and S are the SAME matrix data, but their SPECTRA differ: the U-spectrum couples
S to the all-ones vector J (the score/row-sum direction, converse-ODD); the skew-spectrum does not.

TEST (n=5,6 exhaustive): (a) is the U-spectrum strictly finer than the skew-spectrum? (b) are there
skew-cospectral classes with DIFFERENT |Aut| (skew blind to the covering symmetry) that the U-spectrum
SEPARATES? (c) does each spectrum's degeneracy track |Aut|?
"""
import itertools
import numpy as np
from math import comb
from collections import defaultdict

def all_reps(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    perms = list(itertools.permutations(range(n)))
    seen = {}
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        A = [[0]*n for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b: A[i][j] = 1
            else: A[j][i] = 1
        key = min(tuple(A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)
        if key not in seen: seen[key] = A
    return list(seen.values()), perms

def uspec(A):
    ev = np.linalg.eigvals(np.array(A, float))
    return tuple(sorted(np.round(ev.real, 4))) + tuple(sorted(np.round(ev.imag, 4)))
def skewspec(A):
    n = len(A); S = np.array(A, float) - np.array(A, float).T
    ev = np.linalg.eigvals(S)  # purely imaginary
    return tuple(sorted(np.round(ev.imag, 4)))
def aut(A, perms, n):
    return sum(1 for p in perms if all(A[i][j] == A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j))
def mfas(A, n):
    C2 = comb(n, 2); best = 0
    for p in itertools.permutations(range(n)):
        f = sum(1 for a in range(n) for b in range(a+1, n) if A[p[a]][p[b]])
        best = max(best, f)
    return C2 - best

for n in [5, 6]:
    reps, perms = all_reps(n)
    rows = []
    for A in reps:
        rows.append((uspec(A), skewspec(A), aut(A, perms, n), mfas(A, n)))
    nU = len(set(r[0] for r in rows)); nS = len(set(r[1] for r in rows))
    print(f"=== n={n}: {len(reps)} iso classes ===")
    print(f"  #distinct U-spectra (adjacency) = {nU} / {len(reps)}  ({nU/len(reps):.3f})")
    print(f"  #distinct skew-spectra          = {nS} / {len(reps)}  ({nS/len(reps):.3f})")
    # (b) skew-cospectral groups: does |Aut| vary within them? does U resolve?
    byskew = defaultdict(list)
    for r in rows: byskew[r[1]].append(r)
    skew_blind_aut = 0; u_resolves = 0; groups_autvary = 0
    for sp, grp in byskew.items():
        auts = set(r[2] for r in grp)
        if len(grp) > 1 and len(auts) > 1:
            groups_autvary += 1
            skew_blind_aut += 1  # skew groups together classes of different |Aut|
            if len(set(r[0] for r in grp)) == len(grp): u_resolves += 1
    print(f"  skew-cospectral GROUPS mixing different |Aut| (skew BLIND to symmetry): {groups_autvary}")
    print(f"    of those, U-spectrum fully separates the classes: {u_resolves}")
    # (c) does the spectrum determine |Aut|? (U-cospectral -> same |Aut|? skew-cospectral -> same |Aut|?)
    def determines_aut(idx):
        by = defaultdict(set)
        for r in rows: by[r[idx]].add(r[2])
        return all(len(s) == 1 for s in by.values())
    print(f"  U-spectrum determines |Aut|? {determines_aut(0)}   skew-spectrum determines |Aut|? {determines_aut(1)}")
    # correlation of degeneracy (max eigenvalue multiplicity) with |Aut|
    def maxmult(sp):
        from collections import Counter
        return max(Counter(sp).values())
    uc = np.corrcoef([maxmult(r[0]) for r in rows], [r[2] for r in rows])[0,1]
    sc = np.corrcoef([maxmult(r[1]) for r in rows], [r[2] for r in rows])[0,1]
    print(f"  corr(max-U-spectrum-multiplicity, |Aut|) = {uc:.3f}   corr(max-skew-mult, |Aut|) = {sc:.3f}")
    print()
