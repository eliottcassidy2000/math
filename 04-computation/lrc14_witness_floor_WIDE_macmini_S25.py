#!/usr/bin/env python3
"""
lrc14_witness_floor_WIDE_macmini_S25.py  (mac-mini-2026-06-22-S25)

THE ANALYTIC-CLOSURE GAP: the witness floor G2(P,E) >= m_P for WIDE E
(some offset > 14).  HYP-2830 proved it for BOUNDED E (consec worst).  Here we
test the reduction "wide E equidistributes => larger maxgap => G2 >= bounded".
Reuses the EXACT G2 machinery from lrc14_witness_floor_all_bounded_claudeopus_0622s5.py.

Tests, per k=10,11,12 (worst-witness P from HYP-2830 Phase 1):
  (a) single-far wide  E = consec_{k-1} u {FAR}, FAR=15..30
  (b) dilated wide     E = d*consec_{k-1} (span > 14)
  (c) spread wide      E = consec base + 2 far
vs the BOUNDED consec G2 (HYP-2830) and the floor m_P.  If all wide G2 >= consec
bounded G2 (>= m_P), the wide case reduces to the bounded compact floor.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools

HALF = Fr(1, 14)   # G_P threshold ||p x|| >= 1/14
T1   = Fr(1, 7)    # maxgap threshold
P7   = 7
m_P  = Fr(14249, 252252)

def primitive(E):
    if len(E) < 2: return len(E) == 1 and E[0] != 0
    return reduce(gcd, [E[i]-E[0] for i in range(1, len(E))]) == 1
def phases_at(E, x): return sorted((int(e)*x) % 1 for e in E)
def maxgap(ph):
    if len(ph) <= 1: return Fr(1)
    g = max(b-a for a, b in zip(ph, ph[1:]))
    return max(g, ph[0] + 1 - ph[-1])
def in_GP(Pset, x):
    for p in Pset:
        r = (int(p)*x) % 1
        if min(r, 1-r) < HALF: return False
    return True
def grid(E, Pset):
    bp = {Fr(0), Fr(1)}
    for e in list(E)+list(Pset):
        e = int(e)
        if e == 0: continue
        for ki in range(0, e+1):
            for s in (HALF, -HALF):
                v = (Fr(ki)+s)/e
                if Fr(0) <= v <= Fr(1): bp.add(v)
        for t in range(0, P7*e+1): bp.add(Fr(t, P7*e))
    El = list(E)
    for i in range(len(El)):
        for j in range(i+1, len(El)):
            d = abs(El[i]-El[j])
            if d == 0: continue
            for m in range(0, d+1):
                for s in (T1, -T1):
                    v = Fr(m, d)+s/d
                    if Fr(0) <= v <= Fr(1): bp.add(v)
    return sorted(b for b in bp if Fr(0) <= b <= Fr(1))
def measure_G2(E, Pset):
    pts = grid(E, Pset); G2 = Fr(0)
    for a, b in zip(pts, pts[1:]):
        if b <= a: continue
        mid = (a+b)/2
        if not in_GP(Pset, mid): continue
        if maxgap(phases_at(E, mid)) > T1: G2 += b-a
    return G2

# worst-witness P per k (HYP-2830 Phase 1) + bounded consec G2
WORST_P = {10: [1,11,12], 11: [1,12], 12: [1]}
CONSEC_G2 = {10: Fr(19,49), 11: Fr(3193,8820), 12: Fr(10358,24255)}

print("WIDE-E witness floor: G2(worstP, wide E) >= m_P? and >= bounded consec G2?")
print(f"  m_P = {float(m_P):.5f}\n")
for k in (10, 11, 12):
    P = WORST_P[k]; base = list(range(k-1)); cg = CONSEC_G2[k]
    print(f"k={k}  worstP={P}  bounded-consec G2={float(cg):.5f} (={cg})")
    wides = []
    for FAR in [15, 18, 22, 30]:
        wides.append(("single-far {"+str(FAR)+"}", base + [FAR]))
    for d in [2, 3]:
        E = [d*i for i in range(k)]
        if max(E) > 14: wides.append((f"dilated x{d}", E))
    wides.append(("2-far {16,20}", base + [16, 20]))
    wides.append(("spread", [0,2,4,6,8,10,12,15,18,21][:k]))
    worst_ratio = None
    for name, E in wides:
        E = sorted(set(int(x) for x in E))
        if len(E) != k or not primitive(E): continue
        G2 = measure_G2(E, P)
        okfloor = G2 >= m_P; okcons = G2 >= cg
        r = float(G2/m_P)
        if worst_ratio is None or r < worst_ratio: worst_ratio = r
        print(f"   {name:18s} E={E} G2={float(G2):.5f}  >=m_P:{okfloor} >=consec:{okcons}  ({r:.2f}x m_P)")
    print(f"   => worst wide G2/m_P at k={k}: {worst_ratio:.2f}x\n")
print("=> if every wide G2 >= consec bounded G2 (>> m_P), the WIDE case reduces to the")
print("   bounded compact floor (HYP-2830): wide is SAFER, binding is bounded consec.")
