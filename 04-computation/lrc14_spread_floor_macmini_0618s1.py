#!/usr/bin/env python3
"""
lrc14_spread_floor_macmini_0618s1.py  (mac-mini-2026-06-18-S1)

CRITICAL TEST after the consecutive-extremality correction (FALSE at k=7).
Question: is the EXTREMAL (min pure measure mu) cluster shape BOUNDED-SPREAD?
If yes  -> residual is a finite-dimensional/compact problem (floor > 0 decidable).
If no   -> large-spread shapes drive mu -> 0, reduction fails.

mu(E) = meas{x in [0,1): maxgap{frac(e x): e in E} > 2/7}.

TESTS:
 (A) large-spread family {0,1,...,k-2, M} for M up to 4000: does mu -> 0 or stay bounded?
 (B) broad min-mu search per k (k=7,9,11) over shapes with spread up to 4k: the true min
     and its spread (is the minimizer bounded-spread?).
 (C) is there a UNIVERSAL positive floor mu_min(k) = min over ALL k-subsets? estimate it.
"""
import sys, random, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(61811)

def maxgap(points):
    pts=sorted(p%1.0 for p in points); g=0.0
    for i in range(len(pts)-1): g=max(g,pts[i+1]-pts[i])
    return max(g,(pts[0]+1.0)-pts[-1])
def mu(E,N=200000):
    cnt=0
    for t in range(N):
        x=(t+0.5)/N
        if maxgap([(e*x)%1.0 for e in E])>2.0/7.0: cnt+=1
    return cnt/N

print("="*80)
print("(A) large-spread {0,1,...,k-2,M}: does mu stay bounded as M->inf?  (k=9)")
print("="*80)
k=9; base=list(range(k-1))
cons=mu(list(range(k)))
print(f"   consecutive {{0..{k-1}}}: mu={cons:.5f}")
for M in (10, 30, 100, 400, 1500, 4000):
    E=base+[M]; print(f"   {{0..{k-2}}}∪{{{M}}}: mu={mu(E):.5f}")

print("\n"+"="*80)
print("(B) broad MIN-mu search per k (spread up to 4k); is the minimizer bounded-spread?")
print("="*80)
for k in (7,9,11):
    cons=mu(list(range(k)))
    best=(cons,list(range(k)))
    for _ in range(300):
        sp=random.randint(k-1, 4*k)
        E=sorted(set([0,sp]+random.sample(range(1,sp), min(k-2, sp-1))))
        if len(E)!=k: continue
        m=mu(E, N=80000)
        if m<best[0]: best=(m,E)
    sp=max(best[1])-min(best[1])
    print(f"   k={k}: consecutive mu={cons:.4f}; MIN found mu={best[0]:.4f} at E={best[1]} (spread={sp})")

print("\n"+"="*80)
print("(C) estimate universal floor mu_min(k)=min over ALL k-subsets, k=4..13")
print("="*80)
for k in range(4,14):
    cons=mu(list(range(k)), N=120000)
    best=cons
    for _ in range(120):
        sp=random.randint(k-1, 3*k)
        E=sorted(set([0,sp]+random.sample(range(1,sp), min(k-2,sp-1))))
        if len(E)!=k: continue
        m=mu(E, N=60000)
        best=min(best,m)
    print(f"   k={k:2d}: mu_min ~ {best:.4f}  (consecutive={cons:.4f})")
print("\nDONE.")
