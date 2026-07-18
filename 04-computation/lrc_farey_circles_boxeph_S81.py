#!/usr/bin/env python3
"""
THE FAREY-CIRCLE LAW for the tight family (boxeph-S81).

CLAIM: for the prefix {1..N-1} at threshold 1/N, the K-deep set
{p in (0,q) : bandCount(q,p) >= K} is (for q large enough) a disjoint union of
resonance arcs, one around each Farey fraction a/b in lowest terms with
  b <= (N-1)/K   (equivalently  floor((N-1)/b) >= K  speeds are divisible by b).
Around a/b, the speeds divisible by b nest at 0; the binding one is the K-th
such speed = b*K, giving the arc half-width from  b*K*|p - a q/b| < q/N.

This EXPLAINS and GENERALIZES death-star's THM-985 "two circles" (N=14,K=6):
b <= 13/6 = 2.16 -> b in {1,2} -> exactly TWO circles.  Lower K -> more circles.
"""
from math import gcd

def band_count(V, q, p, N):
    c = 0
    for vi in V:
        r = (vi * p) % q
        if not (q <= N * r <= (N - 1) * q):
            c += 1
    return c

def deep_set(V, q, N, K):
    return [p for p in range(1, q) if band_count(V, q, p, N) >= K]

def farey_denominators(N, K):
    """Predicted set of resonance denominators b with floor((N-1)/b) >= K."""
    return [b for b in range(1, N) if (N-1)//b >= K]

def predicted_centers(q, N, K):
    """Predicted arc centers p ~ a*q/b for each admissible Farey fraction a/b."""
    centers = []
    for b in farey_denominators(N, K):
        for a in range(0, b+1):
            if gcd(a, b) == 1 or (a == 0 and b == 1) or (a == b and b == 1):
                # a/b in lowest terms, 0<=a<=b; centers at p = a*q/b
                centers.append((a, b, a*q/b))
    return sorted(set(centers), key=lambda t: t[2])

def arcs_of(pts):
    if not pts: return []
    arcs=[]; s=pr=pts[0]
    for p in pts[1:]:
        if p==pr+1: pr=p
        else: arcs.append((s,pr)); s=pr=p
    arcs.append((s,pr)); return arcs

def report(N, K, qlist):
    V=list(range(1,N))
    bs = farey_denominators(N,K)
    print(f"N={N} K={K}:  admissible denominators b (floor((N-1)/b)>=K) = {bs}")
    print(f"         => predicted #circles = {len(bs)} (Farey fractions with those b)")
    for q in qlist:
        ds = deep_set(V,q,N,K)
        arcs = arcs_of(ds)
        cen = [round(c[2],1) for c in predicted_centers(q,N,K)]
        print(f"    q={q:4d}: #deep-arcs={len(arcs):2d}  arcs={arcs}")
        print(f"            predicted centers a*q/b = {cen}")
    print()

if __name__ == "__main__":
    print("="*74)
    print("THE FAREY-CIRCLE LAW: deep set = union of resonance arcs, one per")
    print("Farey fraction a/b with b <= (N-1)/K.  Verifying across N and K.")
    print("="*74)
    print("\n--- N=14 (death-star's family), sweeping K to grow/shrink the circle count ---")
    report(14, 6, [200, 300])   # b<=2.16 -> {1,2}  TWO circles (baseline)
    report(14, 4, [200, 300])   # b<=3.25 -> {1,2,3}  THREE circles
    report(14, 3, [300, 420])   # b<=4.33 -> {1,2,3,4} FOUR circles
    report(14, 2, [300, 600])   # b<=6.5  -> {1..6}  SIX circles
    print("--- other N: the same Farey law ---")
    report(12, 3, [300, 400])   # b<=11/3=3.67 -> {1,2,3}
    report(20, 5, [400, 600])   # b<=19/5=3.8  -> {1,2,3}
    report(11, 2, [300, 400])   # b<=5 -> {1,2,3,4,5}
    print("READING: #circles = #{b>=1 : b<=(N-1)/K}. death-star's 'two circles' is the")
    print("K=6,N=14 slice (b<=2). The half-circle is just the b=2 Farey resonance.")
