# opus-2026-07-16-S327 -- HYP-7035 / THM-895:
# (A) THE LEVEL-3 CROSSING: order-3 Bonferroni lower bound
#     uncovered >= muE - S1 + S2 - S3 (Bonferroni: truncation at -S3 is a
#     LOWER bound). Coefficient wall: coercive iff
#     2197 - 338 m' + 26 m'(m'-1) - (8/6) m'(m'-1)(m'-2)... computed exactly.
#     Referee: exact S1, S2, S3 on real radius-8 packets over prefix {1..4}
#     (m' = 8 lifted residues of [12]): verify the exact bound is positive.
# (B) triple sawtooth empirics: rho3(x_i,x_j,x_k) vs (2/13)^3.
# (C) THE TOOTHPICK FERTILITY FIELD: one-vertex extensions n = 5 -> 6:
#     fertility per 5-class vs x-position (does growth localize?).
from fractions import Fraction
from math import gcd, comb
import math, itertools, sys
sys.path.insert(0, '04-computation')

# ---------- (A) coefficient wall (exact)
print("(A) level-3 coefficient wall (units 1/2197, muE = 1):")
for mp in range(7, 14):
    val = 2197 - 338*mp + 52*comb(mp, 2) - 8*comb(mp, 3)
    print(f"   m' = {mp:2d}: 2197 - 338m' + 52C(m',2) - 8C(m',3) = {val:5d} "
          f"{'COERCIVE' if val > 0 else 'dead'}")

# ---------- exact machinery
def safe_arcs(S, delta):
    ivs = [(Fraction(0), Fraction(1))]
    for v in S:
        bands = []
        for k in range(v):
            lo = (Fraction(k) + delta)/v
            hi = (Fraction(k+1) - delta)/v
            if lo < hi: bands.append((lo, hi))
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def comb_teeth_in(x, a, b):
    w = Fraction(1, 13*x)
    out = []
    for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
        lo, hi = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def restrict(arcs, x):
    out = []
    for (a, b) in arcs: out.extend(comb_teeth_in(x, a, b))
    return out

def mu(arcs): return sum(b - a for (a, b) in arcs)

# radius-8 packets over the 4-speed prefix {1,2,3,4}: 8 lifted residues 5..12
P = [1, 2, 3, 4]
E = safe_arcs(P, Fraction(1, 14))
muE = mu(E)
print(f"\n   prefix {P}: mu(E) = {muE} ~ {float(muE):.5f}")
import random
random.seed(11)
print("   exact order-3 Bonferroni on radius-8 packets (S1, S2, S3 exact):")
for trial in range(3):
    xs = [r + 13*random.randint(3, 30) for r in range(5, 13)]
    Ei = [restrict(E, x) for x in xs]
    S1 = sum(mu(e) for e in Ei)
    S2 = Fraction(0); S3 = Fraction(0)
    inter = {}
    for i, j in itertools.combinations(range(8), 2):
        ij = []
        for (a, b) in Ei[i]: ij.extend(comb_teeth_in(xs[j], a, b))
        inter[(i, j)] = ij
        S2 += mu(ij)
    for i, j, k in itertools.combinations(range(8), 3):
        ijk = []
        for (a, b) in inter[(i, j)]: ijk.extend(comb_teeth_in(xs[k], a, b))
        S3 += mu(ijk)
    bound = muE - S1 + S2 - S3
    # actual uncovered
    U = list(E)
    for x in xs:
        out = []
        for (a, b) in U:
            cur = a
            for (lo, hi) in sorted(comb_teeth_in(x, a, b)):
                if lo > cur: out.append((cur, min(lo, b)))
                cur = max(cur, hi)
                if cur >= b: break
            if cur < b: out.append((cur, b))
        U = [iv for iv in out if iv[0] < iv[1]]
    print(f"   packet {trial+1}: BONF3 >= {float(bound):+.5f}  actual = "
          f"{float(mu(U)):.5f}  {'COERCIVE' if bound > 0 else 'not coercive'}")

# ---------- (B) triple sawtooth empirics
print("\n(B) triple overlaps rho3 vs (2/13)^3 = {:.6f} (global, sample):".format((2/13)**3))
full = [(Fraction(0), Fraction(1))]
for (a, b, c) in [(101, 137, 211), (150, 151, 152), (99, 198+1, 300), (77, 143, 169)]:
    t2 = restrict(full, a)
    t2b = []
    for (lo, hi) in t2: t2b.extend(comb_teeth_in(b, lo, hi))
    t3 = []
    for (lo, hi) in t2b: t3.extend(comb_teeth_in(c, lo, hi))
    print(f"   ({a},{b},{c}): rho3 = {float(mu(t3)):.6f}")

# ---------- (C) the toothpick fertility field
print("\n(C) THE FERTILITY FIELD n = 5 -> 6 (one-vertex extensions per class):")
from smith_diagram_of_the_metagraph_opus_S307 import build
from collections import defaultdict
B5, B6 = build(5), build(6)
# representative adjacency for each 5-class; insert vertex 6 in 2^5 ways;
# classify the 6-tournament by B6's classifier via its tiling encoding --
# build the 6-tiling from an adjacency matrix with relabeling to make the
# base path standard: use scores to find a Hamiltonian path via the classic
# insertion (tournaments always have one: greedy insertion).
def ham_path(adj, n):
    path = [0]
    for v in range(1, n):
        placed = False
        for i in range(len(path)):
            if adj[v][path[i]]:
                path.insert(i, v); placed = True; break
        if not placed: path.append(v)
        # fix consistency (insertion sort works for tournaments)
    return path
def tiling_of(adj, n, tiles, tidx):
    p = ham_path(adj, n)
    lam = [0]*n
    for pos, v in enumerate(p): lam[v] = n - pos   # path head -> label n
    t = 0
    for i, (xx, yy) in enumerate(tiles):
        u = lam.index(xx); w = lam.index(yy)
        if adj[u][w]: t |= 1 << i
    return t
tiles5 = [(x, y) for y in range(1, 4) for x in range(5, y+1, -1)]
tiles6 = [(x, y) for y in range(1, 5) for x in range(6, y+1, -1)]
tidx6 = {t: i for i, t in enumerate(tiles6)}
rep5 = {}
for t in range(1 << 6):
    c = B5['cls_of'][t]
    if c not in rep5: rep5[c] = t
fert = {}
for c5, t5 in rep5.items():
    adj = [[False]*6 for _ in range(6)]
    for k in range(2, 6): adj[k-1][k-2] = True
    for i, (xx, yy) in enumerate(tiles5):
        if (t5 >> i) & 1: adj[xx-1][yy-1] = True
        else: adj[yy-1][xx-1] = True
    seen = set()
    for mask in range(32):
        for v in range(5):
            adj[5][v] = bool((mask >> v) & 1)
            adj[v][5] = not adj[5][v]
        t6 = tiling_of(adj, 6, tiles6, tidx6)
        seen.add(B6['cls_of'][t6])
    fert[c5] = len(seen)
xs5 = B5['x_of']
byx = defaultdict(list)
for c, f in fert.items(): byx[xs5[c]].append(f)
print("   x-level -> fertility list (5-classes):")
for x in sorted(byx): print(f"      x = {x:3d}: {sorted(byx[x])}")
tot = sum(fert.values())
print(f"   total extension pairs = {tot}; distinct 6-classes = 56; "
      f"mean fertility = {tot/len(fert):.2f}")
