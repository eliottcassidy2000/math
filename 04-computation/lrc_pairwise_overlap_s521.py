#!/usr/bin/env python3
"""
lrc_pairwise_overlap_s521.py   claudebox-2026-06-01-S521

Exact pairwise overlap of LRC danger arcs, and why the moment method is ruled out.
(Reflection: 07-reflections/lrc-danger-coalescent-formalization-s521.md, T1.)

B_v = {t : ||v t|| < 1/n}, measure 2/n.  Correlation ratio r(a,b) = mu(B_a ∩ B_b)*n^2/4
(=1 iff independent).  Findings (exact Fraction arithmetic):
  * scale-invariant: r(ga,gb)=r(a,b); depends only on coprime (a,b) and n.
  * MAX = doubling pair b=2a: mu=1/n, r=n/4 (grows linearly in n).
  * MIN = negation pair b ≡ -a (mod n): mu=1/(n(n-1)), r->1/2.
  * mu(lonely set)=0 EXACTLY on the LRC-extremal/tight sets (& scalar multiples).
The second moment is blind to the mu=0 extremizers (ordinary variance), so it
cannot prove LRC; the surviving route is the covering-system reformulation.
"""
from fractions import Fraction as F
from math import gcd

def dist(x):
    x = x % 1; return min(x, 1 - x)

def overlap_measure(a, b, n):
    """exact measure of {t in [0,1): ||a t||<1/n and ||b t||<1/n} via cell decomposition."""
    ev = set([F(0), F(1)])
    for v in (a, b):
        for k in range(v):
            ev.add(F(k*n - 1, n*v) % 1); ev.add(F(k*n + 1, n*v) % 1)
    pts = sorted(e for e in ev if 0 <= e <= 1)
    tot = F(0)
    for x, y in zip(pts, pts[1:]):
        if y <= x: continue
        mid = (x + y) / 2
        if dist(F(a)*mid) < F(1, n) and dist(F(b)*mid) < F(1, n): tot += y - x
    return tot

def main():
    print("Correlation ratio r(a,b) = overlap * n^2/4   (1 = independent)\n")
    for n in (5, 7, 9, 12):
        # doubling pair (1,2) and negation pair (1,n-1)
        rd = overlap_measure(1, 2, n) * n*n / 4
        rn = overlap_measure(1, n-1, n) * n*n / 4
        print(f"  n={n:>2}: doubling (1,2) r={float(rd):.3f} (=n/4={n/4:.3f}); "
              f"negation (1,{n-1}) r={float(rn):.4f} (->1/2)")
    print("\n  scale-invariance check (n=7): r(1,2) vs r(3,6):",
          float(overlap_measure(1, 2, 7)*49/4), float(overlap_measure(3, 6, 7)*49/4))
    # confirm max over coprime pairs is the doubling pair
    print("\n  max correlation ratio over coprime pairs (a<b<=2n), per n:")
    for n in (5, 6, 7):
        best = (0, None)
        for a in range(1, 2*n):
            for b in range(a+1, 2*n):
                if gcd(a, b) != 1: continue
                r = overlap_measure(a, b, n) * n*n / 4
                if r > best[0]: best = (float(r), (a, b))
        print(f"    n={n}: max r={best[0]:.3f} at pair {best[1]}  (predicted doubling, r=n/4={n/4:.3f})")

if __name__ == "__main__":
    main()
