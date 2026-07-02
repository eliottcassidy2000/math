#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3953 part 2: adversarial stress of the FLOOR object (*) and the k=8..13 ledger.

 T6: ADVERSARIAL JOINT  J(B, U) = E_x[1_{G_B}(x) * F_U(x)]  — can (B, U) pairs drive it to 0?
     B = bounded small parts (incl. covering-flavored), U = offset-difference patterns tuned
     against them (mod-7 aligned, mod-14, Bohr-aligned with B's lonely arcs).
     hpartA/hfloor need J > 0 for the shapes that occur; the SHAPE hypothesis G2>0 means we only
     need J > 0 WHEN the gap event has positive measure — the adversarial question is whether
     J can vanish while the raw gap event survives (it cannot: J = Int_{G_B} F and G2>0 means
     F > 0 on a positive-measure subset of G_B — J > 0 AUTOMATICALLY. The REAL question for the
     FLOOR route is a uniform lower bound; we measure how low J gets.)
 T7: LARGE-k A-LEDGER: A(U) = E_c[L^c(U)] = E_x[F_U(x)] for k = 8..13 offsets vs the skeleton's
     rhoGlobFloorRat ledger (8152/24255=0.336, 143/420=0.340, 19/49=0.388, 3193/8820=0.362,
     10358/24255=0.427, 477/1078=0.443) and vs witnessMP = 14249/252252 = 0.0565.
     NOTE: rhoGlobFloorRat bounds meas{x in G_P: gap>1/7} for the k-cluster + admissible P;
     A(U) = E_x[F_U] over ALL x is the P-free analogue (upper-layer of the same object).
 T8: THE TWO-WINDOW JOINT: E[1_{G_B} * F_{U1} * F_{U2}] with U1 bounded-diff (narrow window) and
     U2 = V2-scale window (V2 = 50..5000) — the (*)-object of the level-cut; is it bounded below
     uniformly in V2 (scale-invariance like T5) and positive for adversarial alignments?
"""
import sys, itertools, random
from math import gcd, floor
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
rng = random.Random(3954)
BAND = 1.0/14.0

def frac(x): return x - floor(x)
def clip_target(arcs, v, c, band=BAND):
    out = []
    for (lo, hi) in arcs:
        k0 = int(floor(lo*v - c)) - 1; k1 = int(floor(hi*v - c)) + 1
        for k in range(k0, k1+1):
            a = (k + c + band)/v; b = (k + 1 + c - band)/v
            l = lo if lo > a else a
            h = hi if hi < b else b
            if l < h - 1e-15: out.append((l, h))
    return out
def G_arcs(P):
    arcs = [(0.0, 1.0)]
    for p in sorted(P): arcs = clip_target(arcs, p, 0.0)
    return arcs
def F_of_x(U, x, width=1.0/7.0):
    ph = sorted(frac(u*x) for u in U)
    tot = 0.0; k = len(ph)
    for i in range(k):
        g = (ph[(i+1) % k] - ph[i]) % 1.0
        if g > width: tot += g - width
    return tot
def J_of(B, U, ngrid=140000):
    GB = G_arcs(B)
    val = 0.0
    for (lo, hi) in GB:
        n = max(2, int((hi-lo)*ngrid))
        for i in range(n):
            x = lo + (hi-lo)*(i+0.5)/n
            val += F_of_x(U, x)*(hi-lo)/n
    return val, sum(h-l for l, h in GB)

print("="*100)
print(" T6: ADVERSARIAL JOINT  J(B,U) = Int_{G_B} F_U dx  — hunt the minimum")
print("="*100)
worst = None
Bs = [ (1,2,3,4,5), (1,2,3,5,7), (2,3,4,5,6), (1,2,3,4,5,6), (1,3,5,7,9), (2,4,6,8,10,12),
       (1,2,4,8,9,12), (7,14,3,6,12), (1,2,3,4,5,6,7,8) ]
Us = [ tuple(range(0,7)), tuple(range(0,8)), (0,7,14,21,28,35,42), (0,2,4,6,8,10,12),
       (0,1,2,3,4,5,6,7,8,9), (0,1,3,7,12,18,25), (0,14,28,42,56,70,84),
       (0,1,2,14,15,16,28), (0,5,10,15,20,25,30), (0,3,6,9,12,15,18) ]
for B in Bs:
    for U0 in Us:
        U = tuple(sorted(set(U0)))
        if len(U) < 3: continue
        J, mB = J_of(B, U, ngrid=60000)
        tag = ""
        if worst is None or J/mB < worst[0]: worst = (J/mB, J, B, U); tag = " <— new min J/meas"
        # print only notable rows
for _ in range(150):   # random adversarial search around the worst
    B = worst[2]; U = list(worst[3])
    i = rng.randrange(len(U))
    U2 = sorted(set(U[:i] + [max(0, U[i] + rng.choice([-2,-1,1,2,7,-7,14]))] + U[i+1:]))
    if len(U2) != len(U): continue
    J, mB = J_of(B, tuple(U2), ngrid=40000)
    if J/mB < worst[0]: worst = (J/mB, J, B, tuple(U2))
print(f"  minimum found: J/meas(G_B) = {worst[0]:.6f}  (J = {worst[1]:.6f})  at B={worst[2]}, U={worst[3]}")
print("  [J > 0 is AUTOMATIC given G2>0 (J = Int_GB F); the number above measures the uniform-floor room]")

print("\n" + "="*100)
print(" T7: LARGE-k A-LEDGER  A(U) = E_x[F_U]  for k = 8..13  vs rhoGlobFloorRat and witnessMP")
print("="*100)
rho = {8: 8152/24255, 9: 143/420, 10: 19/49, 11: 3193/8820, 12: 10358/24255, 13: 477/1078}
MP = 14249/252252
def A_of(U, ngrid=250000):
    tot = 0.0
    for i in range(ngrid):
        x = (i+0.5)/ngrid
        tot += F_of_x(U, x)/ngrid
    return tot
for k in range(8, 14):
    pool = set()
    for V0 in range(k, min(k+4, 15)):
        from math import comb
        if comb(V0, V0-k) <= 600:
            for C in itertools.combinations(range(1, V0+1), k):
                g = 0
                for c_ in C: g = gcd(g, c_)
                pool.add(tuple(sorted(c_//g for c_ in C)))
    for _ in range(80):
        C = tuple(sorted(rng.sample(range(1, 30), k)))
        pool.add(C)
    pool = list(pool)
    if len(pool) > 150: pool = rng.sample(pool, 150)
    best = None
    for U in pool:
        A = A_of(U, ngrid=60000)
        if best is None or A < best[0]: best = (A, U)
    print(f"  k={k:2d}: min A(U) = {best[0]:.6f} at {best[1]}   [rhoGlobFloorRat = {rho[k]:.4f},"
          f" witnessMP = {MP:.4f}]  clears MP x{best[0]/MP:.1f}")

print("\n" + "="*100)
print(" T8: TWO-WINDOW JOINT  E[1_GB * F_U1 * F_U2(V2-scale)]  — scale-invariance + positivity")
print("="*100)
B = (1,2,3,4,5)
U1 = (0,1,2,3,4,5,6)          # worst narrow window from T5
GB = G_arcs(B); mB = sum(h-l for l,h in GB)
for offs2 in [ (0,1,2,3,4,5,6), (0,1,3,7,12,18,25) ]:
    for V2 in (50, 500, 5000):
        U2 = tuple(V2 - o for o in offs2)
        val = 0.0
        for (lo, hi) in GB:
            n = max(2, int((hi-lo)*160000))
            for i in range(n):
                x = lo + (hi-lo)*(i+0.5)/n
                val += F_of_x(U1, x)*F_of_x(U2, x)*(hi-lo)/n
        print(f"  U1={U1}, U2 = {V2}-{offs2}: E[1_GB F1 F2] = {val:.6f}")
print("DONE.")
