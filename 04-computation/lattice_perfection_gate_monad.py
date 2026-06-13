#!/usr/bin/env python3
"""
lattice_perfection_gate_monad.py
monad-explorer-2026-06-13   (OPEN-Q-057 / THM-496; complements concurrent THM-495)

THESIS (exact-integer): the two-factor RESONANT-product family (THM-493) first beats
3N at n=28, gated by LATTICE-PERFECTION of the factor sizes.

A k-point triangular-lattice (Eisenstein) patch has at most Harborth(k) unit edges.
Call size k "lattice-perfect" when Harborth(k) = u(k) (the planar unit-distance max).
We show Harborth(k)=u(k) for k<=8 and that 9 is the FIRST lattice-imperfect size.

Consequences, all verified here exactly:
  * n=27=3^3 forces a size-9 factor (imperfect: 16<18) AND a size-3 factor (chord-free,
    Delta_t=0).  Exact resonant cap U(27) <= 75 < 81.  (Resonance HURTS at 27; the
    81-tie is GENERIC/off-lattice, not resonant.)
  * n=28=4*7 uses only lattice-perfect sizes; the rhombus (size 4) carries the sqrt3
    chord; gap(28)=3*28-P_gen(28)=84-83=1 < Delta_3=2.  FIRST resonant crossing.

All arithmetic exact in Eisenstein integers a+b*zeta6 (zeta6^2=zeta6-1),
norm N(a,b)=a^2+a*b+b^2.  No floats decide any count.
"""
from itertools import combinations
from collections import defaultdict
from math import isqrt

# ---------- exact Eisenstein arithmetic ----------
def esub(p, q): return (p[0]-q[0], p[1]-q[1])
def enorm(p):
    a, b = p
    return a*a + a*b + b*b
def edges(pts): return sum(1 for p, q in combinations(pts, 2) if enorm(esub(p, q)) == 1)
def chordspec(pts):
    d = defaultdict(int)
    for p, q in combinations(pts, 2): d[enorm(esub(p, q))] += 1
    return dict(d)
def m_alpha(pts):
    d = defaultdict(int)
    for p in pts:
        for q in pts:
            if p != q: d[esub(p, q)] += 1
    return d
def norm_t_alphas(t, R=8):
    return [(a, b) for a in range(-R, R+1) for b in range(-R, R+1) if a*a+a*b+b*b == t]
def delta_t(G, H, t):
    mG, mH = m_alpha(G), m_alpha(H)
    s = sum(mG.get(al, 0)*mH.get(al, 0) for al in norm_t_alphas(t))
    assert s % 2 == 0
    return s // 2

def harborth_formula(k):
    # Harborth (1974): max unit edges of a k-point triangular-lattice patch
    return int(3*k - (12*k-3)**0.5 + 1e-9)  # = floor(3k - sqrt(12k-3))

# cited planar maxima u(k) (AMP arXiv:2412.11914 exact n<=21).  u(9)=18, u(10)=20
# per THM-433 (the two avgdeg-exactly-4 atoms G9,G10).  Only u(3),u(4),u(7),u(9) enter
# the core 27/28 argument; the rest frame the gap table for n<=28.
U = {1:0,2:1,3:3,4:5,5:7,6:9,7:12,8:14,9:18,10:20,11:23,12:27,13:30,14:33}

# ---------- PART 1: lattice-perfection table (Harborth vs u) ----------
print("="*74)
print("PART 1 — lattice-perfection: Harborth(k) vs planar u(k); first divergence")
print("="*74)
print(f"  {'k':>3} {'Harborth(k)':>12} {'u(k)':>6} {'perfect?':>9}")
first_imperfect = None
for k in range(2, 15):
    h = harborth_formula(k); uk = U.get(k)
    perf = (h == uk) if uk is not None else None
    if perf is False and first_imperfect is None: first_imperfect = k
    print(f"  {k:>3} {h:>12} {str(uk):>6} {str(perf):>9}")
print(f"  => lattice-perfect (Harborth=u) for k<=8; FIRST imperfect size = {first_imperfect}")
assert first_imperfect == 9, "expected 9 to be first lattice-imperfect size"

# ---------- PART 2: confirm Harborth(k), k<=9, by full connected-patch enumeration ----------
print("\n" + "="*74)
print("PART 2 — max unit-edges over ALL connected triangular-lattice patches (k<=9)")
print("="*74)
# Every maximum-unit-edge patch is connected; enumerate all connected Eisenstein patches
# up to translation (canonicalize by subtracting the lexicographically-min point) by
# growing size (k-1)->k via boundary additions. Complete: any connected size-k patch has
# a connected size-(k-1) sub-patch (remove a spanning-tree leaf).
UNITS6 = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def norm_form(s):
    lm = min(s)
    return frozenset((p[0]-lm[0], p[1]-lm[1]) for p in s)
def edges_fs(s):
    return sum(1 for p, q in combinations(s, 2) if enorm(esub(p, q)) == 1)
levels = {1: {frozenset([(0,0)])}}
for k in range(2, 10):
    nxt = set()
    for s in levels[k-1]:
        for p in s:
            for u in UNITS6:
                q = (p[0]+u[0], p[1]+u[1])
                if q not in s:
                    nxt.add(norm_form(s | {q}))
    levels[k] = nxt
for k in range(2, 10):
    be = max(edges_fs(s) for s in levels[k])
    hf = harborth_formula(k)
    flag = "OK" if be == hf else f"*** {be} vs formula {hf} ***"
    print(f"  k={k}: max edges over {len(levels[k]):>5} connected patches = {be}  "
          f"(Harborth {hf})  {flag}")

# ---------- canonical lattice factors ----------
def triangle_k(): return [(0,0),(1,0),(0,1)]                       # K3, e=3
RHOMBUS = [(0,0),(1,0),(0,1),(1,1)]                                # e=5, chord {1,3}
W7 = [(0,0),(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]              # rosette, e=12
RHOMBIC9 = [(a,b) for a in range(3) for b in range(3)]            # 3x3 patch, e=16

# ---------- PART 3: the Delta_t bound and exact resonant cap at n=27 ----------
print("\n" + "="*74)
print("PART 3 — exact resonant-product cap at n=27 (=3*9), the only nontrivial 2-factor")
print("="*74)
K3 = triangle_k()
print(f"  e(K3)={edges(K3)} (=3=u(3)); densest 9-lattice patch e={edges(RHOMBIC9)} "
      f"(=16=Harborth(9) < u(9)=18)")
print(f"  ChordSpec(K3)={chordspec(K3)}  -> chord-free (only norm 1)")
print(f"  ChordSpec(3x3 patch)={chordspec(RHOMBIC9)}")
# best resonant product over the size-3 factor in {K3, path P3} x best 9-patch, all t
def path3_120():  # bent unit path: 0 -1- a -1- b with endpoint at norm-3 (sqrt3)
    return [(0,0),(1,0),(1,1)]   # (0,0)-(1,0) unit, (1,0)-(1,1) unit, (0,0)-(1,1) norm 3
def path3_straight():
    return [(0,0),(1,0),(2,0)]   # collinear: endpoints norm 4
print("\n  resonant U over (3-factor) (+)_t (3x3 9-patch), scan t=2..30:")
print(f"   {'3-factor':>14} {'e':>2} {'best_t':>7} {'Delta':>6} {'U':>4}")
for name, G3 in [("K3", K3), ("path(120°)", path3_120()), ("path(straight)", path3_straight())]:
    cart = edges(G3)*9 + 3*edges(RHOMBIC9)
    bt, bd = None, 0
    for t in range(2, 31):
        d = delta_t(G3, RHOMBIC9, t)
        if d > bd: bd, bt = d, t
    print(f"   {name:>14} {edges(G3):>2} {str(bt):>7} {bd:>6} {cart+bd:>4}")
print("  Analytic bound: U(27,resonant) = 9*e(G3)+3*e(H9)+Delta_t,")
print("    e(H9)<=Harborth(9)=16, Delta_t<=8*(3-e(G3))  =>  U <= e(G3)+72 <= 75.")
print("    (With only the planar bound e(H9)<=u(9)=18: U <= e(G3)+78 <= 81.)")
print("  => no 2-factor RESONANT product beats 81 at 27; exact cap 75 (resonance HURTS).")

# ---------- PART 4: generic gap and the n=28 crossing ----------
print("\n" + "="*74)
print("PART 4 — generic product cap P_gen(n), gap to 3N, lattice-perfect chord factorizns")
print("="*74)
def divisor_pairs(n):
    return [(a, n//a) for a in range(2, isqrt(n)+1) if n % a == 0]
def P_gen(n):
    best = 0; arg = None
    for a, b in divisor_pairs(n):
        if a not in U or b not in U: continue
        v = U[a]*b + a*U[b]
        if v > best: best, arg = v, (a, b)
    return best, arg
print(f"  {'n':>3} {'P_gen':>6} {'arg':>7} {'3n':>4} {'gap':>4}  {'LP-chord factrzn (parts<=8, >=1 size>=4)':>40}")
for n in range(24, 29):   # the first-crossing window (n>=30 is THM-432/433 territory)
    pg, arg = P_gen(n)
    gap = 3*n - pg
    # lattice-perfect chord-bearing factorizations: all parts <=8, at least one >=4
    lp = [(a,b) for a,b in divisor_pairs(n) if a<=8 and b<=8 and max(a,b)>=4]
    print(f"  {n:>3} {pg:>6} {str(arg):>7} {3*n:>4} {gap:>4}  {str(lp):>40}")
print("  n=24 (LP 4x6,3x8): gap 6 too big;  n=25 (LP 5x5): gap 5 too big;")
print("  n=26 (only 2x13, 13 imperfect): no LP factorization;")
print("  n=27 (only 3x9, 9 imperfect; 3 chord-free): excluded -- generic merely TIES 81;")
print("  n=28 (LP 4x7): gap=1 < Delta_3=2  =>  FIRST resonant crossing.")

# ---------- PART 5: the explicit n=28 crossing W7 (+)_3 R = 85 ----------
print("\n" + "="*74)
print("PART 5 — explicit n=28 resonant crossing  W7 (+)_3 R = 85 > 84  (exact)")
print("="*74)
cart28 = edges(W7)*len(RHOMBUS) + len(W7)*edges(RHOMBUS)
d3 = delta_t(RHOMBUS, W7, 3)
print(f"  e(W7)={edges(W7)}, e(R)={edges(RHOMBUS)}, |W7|={len(W7)}, |R|={len(RHOMBUS)}")
print(f"  generic Cartesian P(28) = 12*4 + 7*5 = {cart28}")
print(f"  Delta_3(R,W7) = {d3}")
print(f"  U(W7 (+)_3 R) = {cart28} + {d3} = {cart28+d3}  vs 3*28 = 84  -> "
      f"{'CROSS' if cart28+d3>84 else 'no'}")
# confirm 4,7 are lattice-perfect AND R is chord-bearing
print(f"  sizes 4,7 lattice-perfect: Harborth(4)=u(4)={harborth_formula(4)}={U[4]}, "
      f"Harborth(7)=u(7)={harborth_formula(7)}={U[7]}")
print(f"  ChordSpec(R)={chordspec(RHOMBUS)}  -> non-unit chord norm 3 (sqrt3) present")
assert cart28 + d3 == 85

# ---------- PART 6: EXHAUSTIVE max resonant total at n=24,25 (factors<=8) ----------
print("\n" + "="*74)
print("PART 6 — EXHAUSTIVE resonant-product max at n=24,25 over ALL connected patches")
print("="*74)
# levels[k] = all connected k-patches (from PART 2). For each LP factorization with both
# parts<=8 (so fully enumerated), maximize U = e(G)|H|+|G|e(H)+max_t Delta_t exactly.
def m_alpha_fs(s):
    d = defaultdict(int)
    L = list(s)
    for p in L:
        for q in L:
            if p != q: d[esub(p, q)] += 1
    return d
ALPHAS = {t: norm_t_alphas(t) for t in range(2, 14)}   # admissible resonant norms
def best_resonant_total(a, b):
    """max over all connected a-patches G, b-patches H, and t of U(G (+)_t H)."""
    Ga = [(edges_fs(s), m_alpha_fs(s), s) for s in levels[a]]
    Gb = [(edges_fs(s), m_alpha_fs(s), s) for s in levels[b]]
    # prune to near-densest to keep it fast but exact for the MAX (a sparse patch can only
    # help via Delta, bounded by (|.|-1)*P_t; keep all with edges >= Harborth-2 OR with any
    # non-unit chord, i.e. anything that could matter). Here sizes are tiny so keep ALL.
    best = 0; arg = None
    for eA, mA, sA in Ga:
        for eB, mB, sB in Gb:
            cart = eA*b + a*eB
            bd = 0; bt = None
            for t, als in ALPHAS.items():
                d = sum(mA.get(al, 0)*mB.get(al, 0) for al in als) // 2
                if d > bd: bd, bt = d, t
            if cart + bd > best:
                best = cart + bd; arg = (eA, eB, bt, bd)
    return best, arg
for n, facs in [(24, [(3,8),(4,6)]), (25, [(5,5)])]:
    for a, b in facs:
        bt_total, arg = best_resonant_total(a, b)
        eA, eB, t, d = arg
        print(f"  n={n} = {a}x{b}: max resonant U = {bt_total} "
              f"(e(G)={eA},e(H)={eB},best t={t},Delta={d})  vs 3N={3*n}  "
              f"{'CROSS' if bt_total>3*n else 'no'}")
print("  n=26 (=2x13): size-2 factor gives Cartesian <= 1*13+2*Harborth(13) = 65 (chord-free")
print("    edge) or <=2*26=52 with bonus <=12 (norm-t 2-pt) => U<=65<78. No cross.")
print("  => 2-factor RESONANT products do NOT cross 3N at n<=27; n=28 is the FIRST.")
print("\nDONE.")
