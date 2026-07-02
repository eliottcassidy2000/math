#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3953: closing hpartA via the c-averaged ruler. kind-pasteur-2026-07-01-S30 (pure-math session).

hpartA (Lean axiom): 0 < G2(shape) => M(S) >= 1/14, where S = P u L, L = cluster {V - u_i},
G2 = meas{x in G_P : maxgap{frac(u_i x)} > 1/7}.

THE THREE MOVES (this script stress-tests each):
 (1) EXACT c-RULER IDENTITY: at tau_j = (j+c)/V,  ||l_i tau_j|| = ||u_i tau_j - c||  EXACTLY
     (l_i = V - u_i;  l_i tau_j = j + c - u_i tau_j).  So M(S) >= 1/14 iff some ruler point lies in
       G_c := {x in G_P : forall i, ||u_i x - c|| >= 1/14},
     and the lattice count gives:  hpartA  <==  exists c: V*meas(G_c) > arcCount(G_c).
 (2) FUBINI GAP IDENTITY: Int_0^1 meas(G_c) dc = Int_{G_P} F(x) dx,  F(x) = Sum_gaps (gap - 1/7)^+.
     (The free-c measure at fixed x is exactly the total gap surplus.)  So G2 > 0 ==> sup_c meas(G_c)
     >= Int F > 0, and the TIGHT case (sum u_i small vs V * Int F) closes outright.
 (3) THE c-AVERAGE KILLS Sum(m) != 0 RELATIONS: A(U) := E_c[L^c(U)] (average inhomogeneous lonely
     measure of the co-offset set) has Fourier support on {m : Sum m_i u_i = 0 AND Sum m_i = 0}.
     Covering sets kill L^0 but CANNOT kill all targets c.  Conjecture: robust positive floor.

TESTS:
 T1: Fubini identity (numeric, fine grid, several shapes) — should match to grid error.
 T2: end-to-end: real 13-sets S = P u cluster; compare [exists c with a ruler point in G_c] against
     exact M(S) >= 1/14; verify the count lower bound V*meas - arcCount.
 T3: THE LEDGER OF A(U) = E_c[L^c(U)]: min over k-sets U (k = 2..7; exhaustive bounded + adversarial
     covering sets) vs (6/7)^k and vs homogeneous L^0 min.  KEY QUESTION: does the c-average have a
     robust floor where the homogeneous inf collapses (covering adversary)?
 T4: k = 7 criticality: A(U) for near-tiling 7-sets (consecutive, AP, {1..7}-like) — the union-bound
     death row; compare with THM-594-E's Parseval floor scale.
 T5: two-scale decorrelation: E[1_{G_B} * F_T] vs meas(G_B)*E[F_T] as V grows (rate ~ 1/V?), B bounded,
     T = top window — the level-cut induction's glue constant.
"""
import sys, itertools, random
from math import gcd, floor
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
rng = random.Random(3953)

BAND = 1.0/14.0

def frac(x): return x - floor(x)

def circ(x):
    f = frac(x)
    return min(f, 1.0 - f)

# ---------- interval engine (floats) ----------
def clip_target(arcs, v, c, band=BAND):
    """intersect arc list with {x : ||v x - c|| >= band} (teeth centered at (k+c)/v, width 2band/v)."""
    out = []
    for (lo, hi) in arcs:
        k0 = int(floor(lo*v - c)) - 1; k1 = int(floor(hi*v - c)) + 1
        # safe intervals between teeth: [(k+c+band)/v, (k+1+c-band)/v]
        for k in range(k0, k1+1):
            a = (k + c + band)/v; b = (k + 1 + c - band)/v
            l = lo if lo > a else a
            h = hi if hi < b else b
            if l < h - 1e-15: out.append((l, h))
    return out

def lonely_arcs_target(U, cs, box=(0.0, 1.0)):
    """arcs of {x in box : ||u x - c_u|| >= band for all (u, c_u) in zip(U, cs)}"""
    arcs = [box]
    for u, c in sorted(zip(U, cs)):
        arcs = clip_target(arcs, u, c)
        if not arcs: return arcs
    return arcs

def measN(arcs): return sum(h-l for l, h in arcs), len(arcs)

def G_P_arcs(P):
    return lonely_arcs_target(P, [0.0]*len(P))

# ---------- gap functional F(x) = sum (gap - 1/7)^+ of the offset phases ----------
def F_of_x(U, x, width=1.0/7.0):
    ph = sorted(frac(u*x) for u in U)
    tot = 0.0
    k = len(ph)
    for i in range(k):
        g = (ph[(i+1) % k] - ph[i]) % 1.0    # circular gap (wraps for i = k-1)
        if g > width: tot += g - width
    return tot

# ============================================================================================
print("="*100)
print(" T1: FUBINI IDENTITY  Int_c meas(G_c) dc  =  Int_{G_P} F(x) dx   (grid check)")
print("="*100)
for (P, U) in [ ((1,2,3), (0, 3, 7)), ((1,3,4), (0, 2, 9, 11)), ((2,5), (0, 1, 4, 6, 13)) ]:
    GP = G_P_arcs(P)
    # RHS: integrate F over G_P
    NGRID = 200001
    rhs = 0.0
    for (lo, hi) in GP:
        n = max(2, int((hi-lo)*NGRID))
        for i in range(n):
            x = lo + (hi-lo)*(i+0.5)/n
            rhs += F_of_x(U, x) * (hi-lo)/n
    # LHS: average over c of meas(G_c)
    NC = 700
    lhs = 0.0
    for ic in range(NC):
        c = (ic+0.5)/NC
        arcs = GP
        for u in U:
            arcs = clip_target(arcs, u, c) if u != 0 else [ (l,h) for (l,h) in arcs ] if circ(c) >= BAND else []
            if not arcs: break
        m = sum(h-l for l, h in arcs) if arcs else 0.0
        lhs += m/NC
    print(f"  P={P} U={U}:  LHS (avg_c meas G_c) = {lhs:.6f}   RHS (Int_GP F) = {rhs:.6f}   ratio = {lhs/rhs if rhs else float('nan'):.4f}")

# ============================================================================================
print("\n" + "="*100)
print(" T2: END-TO-END on real 13-sets: [exists c, ruler point in G_c] vs exact M >= 1/14")
print("="*100)
def M_exact_approx(S, grid=400000):
    """approximate M(S) = max_t min ||s t|| by fine grid + local refine (good enough to compare to 1/14)."""
    best = 0.0; bt = 0.0
    for i in range(grid):
        t = (i+0.5)/grid
        m = min(circ(s*t) for s in S)
        if m > best: best, bt = m, t
    # refine
    for _ in range(3):
        span = 2.0/grid
        for i in range(200):
            t = bt - span + 2*span*(i+0.5)/200
            m = min(circ(s*t) for s in S)
            if m > best: best, bt = m, t
        span /= 100
    return best, bt

def ruler_test(P, U, V, NC=280):
    """search c-grid: does some ruler point (j+c)/V lie in G_c? Also report the count bound."""
    GP = G_P_arcs(P)
    found = None; best_margin = -1e9; best_c = None
    for ic in range(NC):
        c = (ic+0.5)/NC
        if circ(c) < BAND:  # u=0 constraint (the V-runner itself): ||0*x - c|| = ||c|| >= band
            continue
        arcs = GP
        for u in U:
            if u == 0: continue
            arcs = clip_target(arcs, u, c)
            if not arcs: break
        if not arcs: continue
        m, N = measN(arcs)
        margin = V*m - N
        if margin > best_margin: best_margin, best_c = margin, c
        # actual ruler point check
        for (lo, hi) in arcs:
            j0 = int(floor(lo*V - c)); j1 = int(floor(hi*V - c)) + 1
            for j in range(j0, j1+1):
                x = (j + c)/V
                if lo <= x <= hi:
                    found = (c, j, x); break
            if found: break
        if found: break
    return found, best_margin, best_c

print(f"  {'P':>18} {'cluster (V-u)':>28} {'V':>7} {'ruler finds?':>12} {'V*meas-N':>10} {'M>=1/14?':>9} {'M(S)':>8}")
for (P, offs, V) in [ ((1,2,3,4,5), (0,1,2,3), 500),
                      ((1,2,3,4,5), (0,1,2,3), 5000),
                      ((1,2,3,4,5,6,7,8), (0,1,3,7), 2000),
                      ((1,2,3), (0,1,2,3,5,8,11,14,17,20), 3000),
                      ((2,3,4,5,6), (0,2,4,6,8,10,12,14), 10000) ]:
    S = tuple(P) + tuple(V - u for u in offs)
    found, marg, bc = ruler_test(P, offs, V)
    M, bt = M_exact_approx(S)
    print(f"  {str(P):>18} {str(tuple(V-u for u in offs)):>28} {V:>7} {str(bool(found)):>12} {marg:>10.2f} {str(M >= 1.0/14 - 1e-9):>9} {M:.5f}")

# ============================================================================================
print("\n" + "="*100)
print(" T3: THE c-AVERAGED INHOMOGENEOUS LONELY LEDGER  A(U) = E_c[L^c(U)]  (k = 2..7)")
print("     homogeneous inf collapses on covering sets; does the c-average keep a floor?")
print("="*100)
def A_of_U(U, NC=1200):
    """E_c[L^c(U)]: same shifted target c for every u (the c-ruler structure: targets q_u*c with q_u=1)."""
    tot = 0.0
    for ic in range(NC):
        c = (ic+0.5)/NC
        arcs = lonely_arcs_target(U, [c]*len(U))
        tot += sum(h-l for l, h in arcs)/NC
    return tot

def L0_of_U(U):
    arcs = lonely_arcs_target(U, [0.0]*len(U))
    return sum(h-l for l, h in arcs)

for k in range(2, 8):
    pool = set()
    for V0 in range(k, min(k+5, 13)):
        for C in itertools.combinations(range(1, V0+1), k):
            g = 0
            for c_ in C: g = gcd(g, c_)
            pool.add(tuple(sorted(c_//g for c_ in C)))
    for _ in range(150):
        C = tuple(sorted(rng.sample(range(1, 40), k)))
        pool.add(C)
    pool = list(pool)
    if len(pool) > 320: pool = rng.sample(pool, 320)
    best = None; worstL0 = None
    for U in pool:
        A = A_of_U(U, NC=400)
        L0 = L0_of_U(U)
        if best is None or A < best[0]: best = (A, U)
        if worstL0 is None or L0 < worstL0[0]: worstL0 = (L0, U)
    ind = (6.0/7.0)**k
    print(f"  k={k}: min A(U) = {best[0]:.6f} at U={best[1]}   [(6/7)^k = {ind:.4f}, ratio {best[0]/ind:.3f}]"
          f"   min L^0(U) = {worstL0[0]:.6f} at {worstL0[1]}")

# ============================================================================================
print("\n" + "="*100)
print(" T4: k=7 criticality — A(U) on the union-bound death row (consecutive / AP / dilates)")
print("="*100)
for U in [ tuple(range(1,8)), tuple(range(2,9)), tuple(range(3,10)), (1,2,3,4,5,6,14),
           (2,4,6,8,10,12,14), (1,2,4,8,16,32,64), (7,14,21,28,35,42,49) ]:
    g = 0
    for u in U: g = gcd(g, u)
    Up = tuple(sorted(u//g for u in U))
    A = A_of_U(Up, NC=1000); L0 = L0_of_U(Up)
    print(f"  U={str(U):>28} (prim {Up if Up!=U else '='}):  A = {A:.6f}   L^0 = {L0:.6f}"
          f"   [THM-594-E Parseval floor scale ~ 0.00545]")

# ============================================================================================
print("\n" + "="*100)
print(" T5: TWO-SCALE DECORRELATION  E[1_GB * F_T] vs meas(GB)*E[F_T]  (level-cut glue; rate in V?)")
print("="*100)
B = (1,2,3,4,5)
GB = G_P_arcs(B)
mB, NB = measN(GB)
for offs in [ (0,1,2,3,4,5,6), (0,1,3,7,12,18,25) ]:
    print(f"  B={B} meas(GB)={mB:.5f}; top offsets {offs}:")
    for V in (200, 2000, 20000):
        U = tuple(V - o for o in offs)
        # E[1_GB F_U] over x
        val = 0.0
        for (lo, hi) in GB:
            n = max(2, int((hi-lo)*120000))
            for i in range(n):
                x = lo + (hi-lo)*(i+0.5)/n
                val += F_of_x(U, x)*(hi-lo)/n
        # E[F_U] global
        EFU = 0.0; NG2 = 120000
        for i in range(NG2):
            x = (i+0.5)/NG2
            EFU += F_of_x(U, x)/NG2
        prod = mB*EFU
        print(f"    V={V:>6}: E[1_GB F] = {val:.6f}   meas(GB)*E[F] = {prod:.6f}   rel dev = {(val-prod)/prod if prod else 0:+.4f}")
print("DONE.")
