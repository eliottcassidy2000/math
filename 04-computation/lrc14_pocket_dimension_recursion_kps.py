#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The pocket structure of L_y = Freiman DIMENSION d, and its recursion in k. kps-2026-06-19-S12.

Pocket d = sets of Freiman dimension d (d-dim generalized AP / GAP).
  d=1: AP (excess 0)  -> L_y = L_y(consec), the UNIQUE max (pocket 1).
  d=2: 2-dim GAP      -> pocket 3.
  d=3: 3-dim GAP      -> pocket 4 (look here).
  d=k-1: dissociated  -> pocket 2, L_y -> independent limit.
Relation-lattice rank of a d-dim GAP among the k generators = k-1-d.
GOAL: (1) confirm L_y strictly DECREASES with dimension d (d>=2 all << cap);
      (2) find pocket 4 (d=3); (3) the recursion as k grows (#pockets ~ how many dims fit in k).
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd, prod
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist_p(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=hi-lo
    return p
def g_poly(k,t):
    if k==8: return Fraction((t-1)*(t-2)*(t-4)*(t-5),40)
    if k in(9,10): return Fraction(-(t-2)*(t-3)*(t-6),36)
    return Fraction((t-3)*(t-4),12)
def L_y(E):
    k=len(set(E)); p=dist_p(E); return sum(p[t]*g_poly(k,t) for t in range(7))
def excess(E):
    E=set(E); ss={a+b for a in E for b in E}; return len(ss)-(2*len(E)-1)
def primitive(E): return reduce(gcd,E)==1

def gap(dims, steps):
    """proper GAP {sum x_i*step_i : 0<=x_i<dims_i}."""
    pts=[0]
    for L,s in zip(dims,steps):
        pts=[p+x*s for p in pts for x in range(L)]
    return sorted(set(pts))

caps={8:0.38153,9:0.49426,10:0.6044,12:0.857}

if __name__ == "__main__":
    print("=== L_y by Freiman dimension d (GAP), per k ===\n")
    # dimension-d GAPs sized to ~k. choose steps so it's PROPER (no collisions) and primitive.
    configs = {
      8:  {1: [(8,),(1,)], 2:[(2,4),(1,7)], 3:[(2,2,2),(1,5,17)]},  # d=3: 2x2x2=8
      9:  {1: [(9,),(1,)], 2:[(3,3),(1,11)], 3:[(3,3,1),(1,11,0)] }, # d=3 fallback handled below
      12: {1: [(12,),(1,)], 2:[(3,4),(1,13)], 3:[(2,2,3),(1,5,23)], 4:[(2,2,3,1),(1,5,23,0)]},
    }
    for k in [8,9,12]:
        print(f"--- k={k} (cap={caps[k]}) ---")
        # d=1 (consec)
        C=list(range(k)); print(f"  d=1 AP consec: L_y={float(L_y(C)):.5f} excess={excess(C)} relrank={k-1-1}  (POCKET 1, the MAX)")
        # explicit GAPs for d=2,3,(4)
        gap_specs = {
          8:  [(2,[(2,),(4,)]), ],
          9:  [],
          12: [],
        }
        # d=2
        for (dims,steps,label) in {
          8:  [((2,4),(1,7),'d2'),((2,4),(1,9),'d2b')],
          9:  [((3,3),(1,11),'d2'),((3,3),(1,7),'d2b')],
          12: [((3,4),(1,13),'d2'),((4,3),(1,11),'d2b')],
        }[k]:
            E=gap(dims,steps)
            if len(E)==k:
                print(f"  d=2 GAP {dims} steps{steps}: L_y={float(L_y(E)):.5f} exc={excess(E)} relrank={k-1-2}  {E if k<=9 else ''}")
        # d=3
        for (dims,steps,label) in {
          8:  [((2,2,2),(1,5,17),'d3'),((2,2,2),(1,3,11),'d3b')],
          9:  [((3,3,1),(1,11,0),'deg'),],  # 9 not = 2x2x2; use 3x3 (d2) — true d=3 needs >=8
          12: [((2,2,3),(1,5,23),'d3'),((2,3,2),(1,7,29),'d3b')],
        }[k]:
            E=gap(dims,steps)
            if len(E)==k and primitive(E):
                print(f"  d=3 GAP {dims} steps{steps}: L_y={float(L_y(E)):.5f} exc={excess(E)} relrank={k-1-3}  (POCKET 4)  {E if k<=12 else ''}")
        # d=4 for k=12 (2x2x3 is d3; 2x2x... need 16 for 2^4; use 2x2x3 done. skip)
        print()

    print("=== RECURSION in k: max L_y within each excess band, all k=8..12 (bounded scan) ===")
    print("    (shows the dimension-1 vs dimension>=2 GAP across k; the 'pockets' as k grows)")
    for k in [8,9,10]:
        C=list(range(k)); Lc=L_y(C)
        # max L_y over excess>=1 (the non-AP max) within bounded spread
        box={8:16,9:15,10:14}[k]
        best_nonAP=(Fraction(0),None)
        for tail in itertools.combinations(range(1,box+1),k-1):
            E=(0,)+tail
            if not primitive(E): continue
            if excess(E)==0: continue   # skip APs (only consec among primitive)
            L=L_y(E)
            if L>best_nonAP[0]: best_nonAP=(L,E)
        gap_to_consec=Lc-best_nonAP[0]
        print(f"  k={k}: consec L_y={float(Lc):.5f}  max_(non-AP) L_y={float(best_nonAP[0]):.5f} at {best_nonAP[1]}  "
              f"GAP={float(gap_to_consec):.5f}  (non-AP max < cap {caps[k]}? {best_nonAP[0]<caps[k]})")
