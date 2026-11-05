#!/usr/bin/env python3
"""collective.py -- the collective invariant in its local-limit form.
(1) exact echelon basis of the difference lattice L of the reachable set (level NL) and the base point;
    verification that A_m L subset L and A_m y - y in L for all letters (proof of Y_n subset y_ref + L);
(2) for n = 9..NMAX: multiplicity r_n(N) = #{|u| = n : H(u) = N} on the solid window of the level (min, median,
    fraction of N with r_n(N) >= 1), and residue equidistribution of H(u) mod q (q = 3,5,7,8,9,16) over the window
    (max deviation of class frequencies from uniform among admissible odd classes);
(3) proximality: the dominant eigenvalue of each letter and of the two-letter products is simple (numpy)."""
import sys
from collections import Counter, defaultdict
import numpy as np, flint
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 13
def load_intrep(path):
    with open(path) as f:
        r, nl = map(int, f.readline().split())
        y0 = list(map(int, f.readline().split())); hc = list(map(int, f.readline().split()))
        A = {}
        for _ in range(nl):
            mask = int(f.readline().strip()); A[mask] = [list(map(int, f.readline().split())) for _ in range(r)]
    return r, y0, hc, A
r, y0, hc, A = load_intrep('span2_intrep.txt')
Am = {m: np.array(M, dtype=np.int64) for m, M in A.items()}
masks = sorted(Am); hidx = hc.index(1)
levels = []; cur = Counter({tuple(y0): 1})
for n in range(1, NMAX + 1):
    nxt = Counter()
    for y, c in cur.items():
        v = np.array(y, dtype=np.int64)
        for m in masks:
            if n >= 2 and (m >> (n - 1)): continue
            if n == 1 and m: continue
            nxt[tuple((Am[m] @ v).tolist())] += c
    cur = nxt; levels.append(cur)
    print(f"n={n}: {len(cur)} states, {sum(cur.values())} words", flush=True)
# (1) lattice
NL = 7
Y = sorted(levels[NL - 1]); y_ref = Y[0]
diffs = [[a - b for a, b in zip(y, y_ref)] for y in Y[1:]]
Hn = flint.fmpz_mat(diffs).hnf()
rows = [[int(Hn[i, j]) for j in range(Hn.ncols())] for i in range(Hn.nrows()) if any(int(Hn[i, j]) for j in range(Hn.ncols()))]
print(f"(1) base point y_ref (level {NL}) = {list(y_ref)}")
print("(1) echelon basis of L (rows):")
for rw in rows: print("    ", rw)
def in_lattice(v, rows):
    v = list(v)
    for rw in rows:
        p = next(c for c, x in enumerate(rw) if x)
        q, rem = divmod(v[p], rw[p])
        if rem: return False
        v = [a - q * b for a, b in zip(v, rw)]
    return not any(v)
ok = all(in_lattice((Am[m] @ np.array(rw, dtype=np.int64)).tolist(), rows) for m in masks for rw in rows)
ok2 = all(in_lattice((Am[m] @ np.array(y_ref, dtype=np.int64) - np.array(y_ref)).tolist(), rows) for m in masks)
ok3 = all(in_lattice([a - b for a, b in zip(y, y_ref)], rows) for lev in levels[NL - 1:] for y in lev)
later = all(in_lattice([a - b for a, b in zip(y, y_ref)], rows) for lev in levels[NL - 1:NMAX] for y in lev)
print(f"(1) A_m L subset L for all letters: {ok}; A_m y_ref - y_ref in L: {ok2}; every state of levels {NL}..{NMAX} lies in y_ref + L: {ok3}")
# congruence description: for each coordinate, modulus of the pivot column (diagonal) when the basis is diagonal-like
print("(1) pivot moduli per column:", {next(c for c, x in enumerate(rw) if x) + 1: next(x for x in rw if x) for rw in rows})
# (2) multiplicities and residues
for n in range(9, NMAX + 1):
    lev = levels[n - 1]
    rN = Counter()
    for y, c in lev.items(): rN[y[hidx]] += c
    vals = sorted(rN)
    # solid window [49, d] of this level
    d = 47
    while d + 2 in rN: d += 2
    window = list(range(49, d + 1, 2))
    mult = [rN[N] for N in window]
    print(f"(2) n={n}: solid window [49, {d}], r_n(N): min {min(mult)}, median {sorted(mult)[len(mult)//2]}, max {max(mult)}; words {sum(rN.values())}")
    line = []
    for q in (3, 5, 7, 8, 9, 16):
        cls = Counter()
        for N in window: cls[N % q] += rN[N]
        adm = [a for a in range(q) if any((N % q) == a for N in window)]
        tot = sum(cls.values()); freqs = [cls[a] / tot for a in adm]
        line.append(f"q={q}: max|f-1/k|={max(abs(f - 1/len(adm)) for f in freqs):.3f}")
    print("      residues:", "; ".join(line))
# (3) proximality
names = {0: '--', 1: '1-', 2: '-2', 3: '12'}
print("(3) dominant eigenvalues (modulus of the two largest):")
for m in (1, 2, 3):
    ev = sorted(np.abs(np.linalg.eigvals(Am[m].astype(float))), reverse=True)
    print(f"    {names[m]}: {ev[0]:.6f}, {ev[1]:.6f}")
for m1 in (1, 2, 3):
    for m2 in (1, 2, 3):
        if m1 >= m2: continue
        ev = sorted(np.abs(np.linalg.eigvals((Am[m2] @ Am[m1]).astype(float))), reverse=True)
        print(f"    {names[m1]}.{names[m2]}: {ev[0]:.6f}, {ev[1]:.6f}")
