"""
mac-mini-2026-07-07-S51 (HYP-5097) -- THE WINDING-BIT CYLINDER MODEL and its relation
to TRANSITIVITY (owner: extend the two-circle model; how does the cocylinder relate to
transitivity).

MODEL (derived): vertices on two concentric circles (inner A at angles i/m, outer B at
angles (j + tw)/n', tw = global twist); every edge = a radially-monotone curve; a curve
from inner angle a to outer angle b with winding lift d = (b + w) - a (w in Z, binary
near optimum).  TWO such curves cross exactly #(Z cap open interval (x, x + delta))
times, x = a1 - a2, delta = d1 - d2 (solve a1 + s d1 == a2 + s d2 mod 1, s in (0,1)).
Within-circle chords live in their own face: their crossings are FORCED by the cyclic
order: C(m,4) + C(n',4) (every 4 points on a circle contribute exactly one chord
crossing).  So for K_n with split m + n':
    Q_cyl = C(m,4) + C(n',4) + CROSS(windings, twist).

THE TOURNAMENT READING (the cylinder tiling model): fix transitive chains on each
circle (the two spines -- the cylinder's analog of the book's single Redei spine);
the FREE data = the m*n' cross-arc orientations := the winding bits (w=0: A-beats-B
"aligned"; w=1: B-beats-A "wound").  Cube 2^{mn'} -> tournaments; Q_cyl is a quadratic
form on the cube.

QUESTIONS:
 Q1 (validation): does min Q_cyl over windings+twist at balanced splits reproduce the
    2-page/Guy optimum Z(n) (cylindrical drawings are conjectured optimal)?
 Q2 (THE TRANSITIVITY QUESTION): corr(Q_cyl, log H) over the winding cube -- the book
    model gave corr < 0 (S50: crossing-minimal = cycle-rich).  Does the cylinder
    INVERT (aligned/transitive-ish windings optimal)?  Where does the fully transitive
    tournament (all bits aligned) sit in the Q_cyl landscape?
"""
import numpy as np
from itertools import combinations
from math import comb
import random as rnd
rnd.seed(51)

def cross_count(a1, d1, a2, d2):
    x = a1 - a2; delta = d1 - d2
    lo, hi = min(x, x+delta), max(x, x+delta)
    import math
    return max(0, math.floor(hi) - math.floor(lo) - (1 if hi == math.floor(hi) else 0)) \
        if hi > lo else 0

def cyl_Q(m, np_, w, tw):
    """crossings among cross edges only; w[i][j] in {0,-1} winding; returns int."""
    edges = [(i, j) for i in range(m) for j in range(np_)]
    tot = 0
    for (i1, j1), (i2, j2) in combinations(edges, 2):
        a1, b1 = i1/m, (j1 + tw)/np_
        a2, b2 = i2/m, (j2 + tw)/np_
        d1 = b1 + w[i1][j1] - a1
        d2 = b2 + w[i2][j2] - a2
        tot += cross_count(a1, d1, a2, d2)
    return tot

def Z(n):
    return (n//2)*((n-1)//2)*((n-2)//2)*((n-3)//2)//4

print("=== Q1: validation -- balanced-split cylinder minimum vs Z(n) ===")
for n in (5, 6, 7, 8):
    m = n//2; np_ = n - m
    forced = comb(m, 4) + comb(np_, 4)
    best = 10**9
    # search: all-0, spiral-like, local search from random
    for trial in range(60):
        if trial == 0:
            w = [[0]*np_ for _ in range(m)]
        else:
            w = [[rnd.choice([0, -1]) for _ in range(np_)] for _ in range(m)]
        tw = rnd.random()
        cur = cyl_Q(m, np_, w, tw)
        for it in range(150):
            i, j = rnd.randrange(m), rnd.randrange(np_)
            w[i][j] = -1 - w[i][j]   # toggle 0 <-> -1
            v = cyl_Q(m, np_, w, tw)
            if v <= cur: cur = v
            else: w[i][j] = -1 - w[i][j]
        best = min(best, cur + forced)
    print(f"  n={n} (split {m}+{np_}): forced = {forced}, best Q_cyl = {best}, Z(n) = {Z(n)}"
          f"  {'MATCH' if best == Z(n) else f'delta {best - Z(n):+d}'}")

print("\n=== Q2: transitivity -- corr(Q_cyl, log H) over the winding cube ===")
def ham_count_bits(m, np_, wbits):
    """tournament: A = 0..m-1 transitive chain (i beats i' if i < i'), B likewise,
    cross: A_i beats B_j iff wbits[i][j] == 0 (aligned)."""
    n = m + np_
    adj = np.zeros((n, n), dtype=bool)
    for i in range(m):
        for i2 in range(i+1, m): adj[i, i2] = True
    for j in range(np_):
        for j2 in range(j+1, np_): adj[m+j, m+j2] = True
    for i in range(m):
        for j in range(np_):
            if wbits[i][j] == 0: adj[i, m+j] = True
            else: adj[m+j, i] = True
    full = 1 << n
    dp = np.zeros((full, n), dtype=np.int64)
    for v in range(n): dp[1 << v, v] = 1
    for S in range(full):
        for v in range(n):
            if not dp[S, v] or not (S >> v) & 1: continue
            for u in range(n):
                if (S >> u) & 1: continue
                if adj[v, u]: dp[S | (1 << u), u] += dp[S, v]
    return int(dp[full-1].sum())

for n in (6, 7):
    m = n//2; np_ = n - m
    nb = m*np_
    Qs = []; Hs = []
    trans_Q = None
    for mask in range(1 << nb):
        w = [[0]*np_ for _ in range(m)]
        for b in range(nb):
            if (mask >> b) & 1: w[b // np_][b % np_] = -1
        q = cyl_Q(m, np_, w, 0.25) + comb(m,4) + comb(np_,4)
        h = ham_count_bits(m, np_, [[(-w[i][j]) for j in range(np_)] for i in range(m)])
        Qs.append(q); Hs.append(h)
        if mask == 0: trans_Q = q
    Qs = np.array(Qs); Hs = np.array(Hs, dtype=float)
    r = np.corrcoef(Qs, np.log(Hs))[0, 1]
    print(f"  n={n} (cube 2^{nb}): corr(Q_cyl, log H) = {r:+.4f} "
          f"[book model S50 gave {-0.76 if n==6 else -0.59:+.2f}]")
    print(f"    all-aligned (transitive-ish) Q = {trans_Q}; global min Q = {int(Qs.min())} "
          f"at H = {sorted(set(int(h) for q,h in zip(Qs,Hs) if q == Qs.min()))[:6]}; "
          f"max-H tournaments' Q = {sorted(set(int(q) for q,h in zip(Qs,Hs) if h == Hs.max()))}")
    print(f"    H at all-aligned = {int(Hs[0])}; H range on cube = [{int(Hs.min())}, {int(Hs.max())}]")
