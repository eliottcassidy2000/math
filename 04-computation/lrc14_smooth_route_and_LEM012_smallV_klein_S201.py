#!/usr/bin/env python3
"""
klein-2026-07-09-S201  ("work on completely closing the covering case")

A rigorous audit of the GOOD-PERIOD leg on the EXTREMAL / RESONANT small-Vmax corner.
Two "closures" of the near-AP branch are REFUTED there, both by the SAME small-V mechanism,
and the correct home of those clusters (the density-floor leg) is confirmed.

Good period j (THM-527): j in {1,..,V-1} with maxgap{ e*j mod V : e in E } > V/7.
  E_x[maxgap]   = continuum mean over x in [0,1)  (opus-S170's smooth surrogate)
  E_grid[maxgap]= grid mean over j=1..V-1         (what max>=mean actually needs)
  best/V        = max over j of maxgap/V           (EXISTENCE: >1/7 iff a good period exists)
  mu_good       = measure{ x : maxgap(xE) > 1/7 }  (density-floor leg quantity, ~ mu_{1/7})

FINDINGS
  (1) opus-S170 smooth-MEAN route is REFUTED. It claims E_x>1/7 (+ disc<=0.006) => good period,
      "INCLUDING the tight AP {1..13}". But the tight AP {0..12} at its ruler V=13 is the EXTREMAL
      LRC(14) instance: EVERY j gives maxgap = 1/13 (13 pts equidistribute at every rational j/13),
      so E_grid = 1/13 = 0.077 < 1/7 and NO good period -- yet E_x = 0.211. disc = 0.137 (not 0.006):
      at the resonant ruler V=13 the grid lands on maxgap's equidistribution NULLS, the OPPOSITE of E_x.
      => E_x > 1/7 does NOT imply a good period.  (And V=33 below: E_grid<1/7 yet good period EXISTS,
         so max>=mean is one-way; the smooth mean is neither necessary nor sufficient at resonant V.)
  (2) LEM-012 (klein-S196, canon PROVED) is BUGGY as stated. Hyp "V > max E" is too weak: its
      Dirichlet step yields j with ||jd/V|| < 1/Q, and for {0..12}, d=1, V=13, Q=14 the only such j
      is j=13 == V == 0 (mod V) -- the EXCLUDED trivial period. Correct hypothesis: V >= Q+1.
      At V=15(=Q+1) Dirichlet gives a valid nonzero j and a good period exists. Small V (<=Q) is
      density-floor territory, NOT good-period.
  (3) The correct home of the no-good-period resonant clusters is the DENSITY-FLOOR leg:
      mu_good({0..12}) = 0.44 >= bar_13 (comfortably). No hole in the covering case -- only in the
      two good-period "closures".

LESSON (same as klein-S200 arc-count): EXISTENCE IS A MAX, NOT A MEAN/COUNT. Certify a good period by
exhibiting the best j (Dirichlet/collapse, LEM-012 with V>=Q+1 / LEM-013), or route to the density
floor -- never by an average, which the extremal/resonant cluster fools.
"""
import numpy as np
from math import gcd, ceil
from functools import reduce
INV7=1/7

def maxgap_res(res,V):
    p=np.unique(np.sort(np.array(res)%V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()
def best_period(E,V):
    E=np.array(sorted(set(E))); return max(maxgap_res(E*j,V) for j in range(1,V))/V
def Ex(E,Nx=2000000):
    E=np.array(sorted(set(E))); x=(np.arange(Nx)+0.5)/Nx
    ph=np.sort(np.mod(np.outer(x,E),1.0),axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1); return g.max(axis=1).mean()
def Egrid(E,V):  # j=1..V-1 (exclude trivial j=0)
    E=np.array(sorted(set(E))); j=np.arange(1,V); x=j/V
    ph=np.sort(np.mod(np.outer(x,E),1.0),axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1); return g.max(axis=1).mean()
def mu_good(E,Nx=2000000):
    E=np.array(sorted(set(E))); x=(np.arange(Nx)+0.5)/Nx
    ph=np.sort(np.mod(np.outer(x,E),1.0),axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1); return (g.max(axis=1)>INV7).mean()
def longest_AP(E):
    Es=sorted(set(E)); S=set(Es); b=1
    for i in range(len(Es)):
        for j in range(i+1,len(Es)):
            d=Es[j]-Es[i]; L=2; x=Es[j]+d
            while x in S: L+=1; x+=d
            b=max(b,L)
    return b

print(f"1/7 = {INV7:.5f}   (good period exists  <=>  best/V > 1/7)\n")
print(f"{'cluster':>30} {'V':>4} {'L':>3} {'E_x':>6} {'E_grid':>7} {'disc':>7} {'best/V':>7} {'good?':>6} {'mu_good':>7}")
cases=[
 ("tightAP {0..12} EXTREMAL", list(range(13)), 13),
 ("tightAP {0..12} @V=Q+1",   list(range(13)), 15),
 ("3-struct {0,3..30,10,32}", [0,3,6,9,10,12,15,18,21,24,27,30,32], 33),
 ("2-struct",                 [0,2,4,6,8,9,10,12,14,16,18,20,21], 22),
 ("7-struct dissoc (S200)",   [0,7,14,21,26,29,37,44,51,58,67,75,82], 91),
]
for nm,E,V in cases:
    E=[e-min(E) for e in sorted(set(E))]
    ex,eg,bp,mu,L=Ex(E),Egrid(E,V),best_period(E,V),mu_good(E),longest_AP(E)
    print(f"{nm:>30} {V:>4} {L:>3} {ex:>6.3f} {eg:>7.3f} {eg-ex:>+7.3f} {bp:>7.3f} {str(bp>INV7):>6} {mu:>7.3f}")

print("\nLEM-012 Dirichlet j for {0..12} (L=13,k=13 => Q=14):")
def dirichlet_j(d,V,Q):
    for j in range(1,Q+1):
        r=(j*d)%V
        if min(r,V-r)/V < 1/Q: return j,r
    return None
for V in (13,15):
    j,r=dirichlet_j(1,V,14)
    print(f"  V={V:>2}: Q=14, Dirichlet j={j}, jd mod V={r}, excluded(j==0 mod V)? {j%V==0}"
          f"   {'<-- BUG: hyp V>maxE holds but j excluded' if j%V==0 else '<-- ok (V>=Q+1)'}")

print("\nCONCLUSION: good-period leg needs V>=Q+1 (LEM-012 corrected) / a REAL good period (LEM-013);")
print("the small-V resonant corner (extremal tight AP) has NO good period and closes via mu_good>=bar.")
print("opus-S170 smooth-MEAN route: REFUTED (E_x>1/7 =/=> good period; disc=0.137 not 0.006).")
