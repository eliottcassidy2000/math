#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S6 -- HYP-4292: THE 7-SPREAD LATTICE BRACKET CENSUS.

The (A)-residual after kps-S20's torus_split_rung (which formally kills all
<= 6-lifted couplings = support-<=6 directions, via GL2(Z) normalization):
the genuine 7-spread rank-2 lattices.

CHARACTERIZATION (mac-mini-S6, derived): 7-spread <=> the 12 pair-vectors
(u_i, v_i) fall into direction-classes each of size <= 5 (so >= 3 classes).
Within class D (primitive normal d_D), coords have (u_i,v_i) = c_i * d_D, so
    u_i t + v_i s = c_i * (d_D . (t,s)) = c_i * tau_D,   tau_D = d_D . (t,s).
=> a 7-spread torus is >= 3 COUPLED <=5-runner LR systems in transversal
   linear forms tau_D.
    M(U) = max_{(t,s)} min_D min_{i in D} || c_i tau_D ||.

QUESTION: does ANY 7-spread lattice have M(U) in (1/13, 2/25]?  If none:
(A) closes (with the support-6 kill), and via the sibling's unification the
|S| >= 7 (C)-residual inherits it.

RIGOR: f(t,s) = min_i ||u_i t + v_i s|| is (max|u_i|, max|v_i|)-Lipschitz;
grid_max <= M(U) <= grid_max + (L1 h1 + L2 h2)/2.  Brackets need only CLEAR
the window; doubled-grid + local refinement on straddles.

FAMILIES (all genuinely >= 3 classes, verified max-class <= 5):
 (1) 3 classes of 4, directions from {(1,0),(0,1),(1,1),(1,-1),(1,2),(2,1)},
     speeds AP-blocks {1,2,3,4} / {5,6,7,8} / {9,10,11,12} and permutations;
 (2) 4 classes of 3, tighter speeds;
 (3) 2 classes of 5 + 1 class of 2 (the boundary max-class = 5 case);
 (4) 6 classes of 2 (maximally spread);
 (5) random >= 3-class lattices with speeds up to 12 and up to 24;
 (6) the DELIBERATELY-TIGHT constructions: speeds within each class = a
     dilated-AP so each class is individually as tight as possible, directions
     chosen for maximal transversal coverage (the triangular-tiling threat).
"""
from fractions import Fraction as F
from math import gcd
import random, time
from itertools import combinations, product

T0 = time.time()
def log(m=""):
    print(m, flush=True)
random.seed(6)

LO, HI = F(1, 13), F(2, 25)
LOf, HIf = float(LO), float(HI)

def dist1(x):
    x = x - int(x)
    if x < 0:
        x += 1
    return min(x, 1 - x)

def build_uv(classes):
    """classes: list of (direction (p,q), [speeds]).  Returns u, v (len 12)."""
    u, v = [], []
    for (p, q), speeds in classes:
        for c in speeds:
            u.append(c * p)
            v.append(c * q)
    return u, v

def max_class_size(u, v):
    """direction-class sizes of the pair-vectors (up to sign/scale)."""
    from collections import Counter
    dirs = Counter()
    for a, b in zip(u, v):
        if a == 0 and b == 0:
            return 99          # zero pair-vector: degenerate (not proper)
        g = gcd(abs(a), abs(b))
        p, q = a // g, b // g
        if p < 0 or (p == 0 and q < 0):
            p, q = -p, -q
        dirs[(p, q)] += 1
    return max(dirs.values())

def bracket(u, v, N=360, depth=0):
    L1, L2 = max(abs(a) for a in u), max(abs(b) for b in v)
    best = 0.0
    bt = bs = 0.0
    for i1 in range(N):
        t = i1 / N
        ut = [a * t for a in u]
        for i2 in range(N):
            s = i2 / N
            m = 1.0
            for k in range(12):
                d = dist1(ut[k] + v[k] * s)
                if d < m:
                    m = d
                    if m <= best:
                        break
            if m > best:
                best = m
                bt, bs = t, s
    slack = (L1 / N + L2 / N) / 2
    lower, upper = best, best + slack
    if depth < 3 and not (lower > HIf or upper < LOf):
        # local refine near (bt,bs)
        span = 3.0 / N
        M = 240
        bb = lower
        bbt, bbs = bt, bs
        for i1 in range(M + 1):
            t = bt - span / 2 + span * i1 / M
            ut = [a * t for a in u]
            for i2 in range(M + 1):
                s = bs - span / 2 + span * i2 / M
                m = 1.0
                for k in range(12):
                    d = dist1(ut[k] + v[k] * s)
                    if d < m:
                        m = d
                        if m <= bb:
                            break
                if m > bb:
                    bb = m
                    bbt, bbs = t, s
        lower = max(lower, bb)
        slack2 = (L1 * span / M + L2 * span / M) / 2
        upper = min(upper, max(bb + slack2, best + slack))
        if not (lower > HIf or upper < LOf):
            return bracket(u, v, N * 2, depth + 1)
    return lower, upper

def verdict(lo, up):
    if lo > HIf:
        return "SAFE-ABOVE"
    if up <= LOf:
        return "BELOW-13TH"
    if up <= HIf and lo > LOf:
        return "IN-WINDOW"
    return "UNRESOLVED"

DIRS = [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)]

cases = []
# (1) 3 classes of 4
sp4 = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]
for d3 in combinations(DIRS, 3):
    cases.append(("3x4-AP", [(d3[j], sp4[j]) for j in range(3)]))
# (1b) 3 classes of 4, interleaved AP speeds
sp4i = [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]
for d3 in combinations(DIRS[:6], 3):
    cases.append(("3x4-interleave", [(d3[j], sp4i[j]) for j in range(3)]))
# (2) 4 classes of 3
sp3 = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
for d4 in combinations(DIRS[:6], 4):
    cases.append(("4x3", [(d4[j], sp3[j]) for j in range(4)]))
# (3) 2 classes of 5 + 1 of 2 (boundary)
for d3 in combinations(DIRS[:6], 3):
    cases.append(("5-5-2", [(d3[0], [1,2,3,4,5]), (d3[1], [6,7,8,9,10]), (d3[2], [11,12])]))
# (4) 6 classes of 2
sp2 = [[1,2],[3,4],[5,6],[7,8],[9,10],[11,12]]
cases.append(("6x2", [(DIRS[j], sp2[j]) for j in range(6)]))
# (6) tight: each class a dilated AP {1,2,3,4} scaled, transversal dirs
for scale in ([1,1,1],[1,2,3],[1,3,5]):
    cases.append((f"3x4-dilated{scale}", [((1,0),[scale[0]*c for c in [1,2,3,4]]),
                                          ((0,1),[scale[1]*c for c in [1,2,3,4]]),
                                          ((1,1),[scale[2]*c for c in [1,2,3,4]])]))
# (5) random >= 3-class
def rand_case(maxspeed, maxdir):
    m = random.randint(3, 6)
    sizes = [0]*m
    for _ in range(12):
        sizes[random.randrange(m)] += 1
    if max(sizes) > 5 or 0 in sizes:
        return None
    dirs = random.sample([(p,q) for p in range(0,maxdir+1) for q in range(-maxdir,maxdir+1)
                          if (p,q)!=(0,0) and gcd(abs(p),q if q>=0 else -q) <= 1 and (p>0 or q>0)], m)
    cl = []
    for j in range(m):
        cl.append((dirs[j], random.sample(range(1, maxspeed+1), sizes[j])))
    return cl
nr = 0
while nr < 400:
    c = rand_case(random.choice([12, 18, 24]), random.choice([2, 3]))
    if c:
        cases.append(("rand", c))
        nr += 1

log(f"7-spread census: {len(cases)} lattices; window (1/13, 2/25] = ({LOf:.6f}, {HIf:.6f}]")
stats = {}
flagged = []
mins = []
for name, cl in cases:
    u, v = build_uv(cl)
    if len(u) != 12:
        continue
    mc = max_class_size(u, v)
    if mc > 5:
        continue                 # not 7-spread (support-6 kill handles it)
    lo, up = bracket(u, v)
    vd = verdict(lo, up)
    stats[vd] = stats.get(vd, 0) + 1
    mins.append((lo, name, cl))
    if vd in ("IN-WINDOW", "UNRESOLVED", "BELOW-13TH"):
        flagged.append((name, cl, lo, up, vd))
        log(f"  {vd}: {name} bracket=({lo:.6f},{up:.6f})  classes={cl}")
mins.sort()
log(f"\ncensus verdicts: {stats}")
log(f"lowest 6 M-brackets seen (all should be SAFE-ABOVE):")
for lo, name, cl in mins[:6]:
    log(f"   M >= {lo:.6f}  [{name}]")
log("VERDICT: " + ("NO 7-spread lattice in the window -- (A)'s residual is clean on the census"
                   if not flagged else f"{len(flagged)} FLAGGED -- critical, see above"))
log(f"[t = {time.time()-T0:.0f}s]")
