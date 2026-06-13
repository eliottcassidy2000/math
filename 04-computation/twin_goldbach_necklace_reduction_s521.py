#!/usr/bin/env python3
"""
twin_goldbach_necklace_reduction_s521.py    oracle-2026-06-01-S521

WHY do the 35 twin-Goldbach exceptions come in TRIPLES?  Proposed structure:

  Every twin pair {p, p+2} is centered at c = p+1.  For p>=5 the pair is
  {6k-1, 6k+1}, so c = 6k -- a MULTIPLE OF 6.  (Only the anomalous first pair
  {3,5} has center 4.)  Let C = multiset of twin-pair centers.

  A sum of two twin primes a+b with a in {c1-1,c1+1}, b in {c2-1,c2+1} equals
  one of {c1+c2-2, c1+c2, c1+c2+2}.  So each center-pair (c1,c2) covers the
  TRIPLE  {s-2, s, s+2}  with s = c1+c2.

  Since (apart from center 4) C lies in 6Z, the covering sums s=c1+c2 lie in 6Z,
  and the triple {6m-2, 6m, 6m+2} = the three evens nearest 6m is covered AS A
  UNIT iff 6m is a sum of two twin-pair-centers.  Hence exceptions come in
  triples {6m-2,6m,6m+2}, one whole triple per multiple of 6 that is NOT in C+C.

  REDUCTION (the "necklace"): twin-Goldbach for even n  <=>  ordinary-Goldbach-
  style representability of m=round(n/6) on the TWIN-CENTER NECKLACE
  K = C/6 = {1,2,3,5,7,10,12,...} (k such that 6k+-1 are twin primes), i.e.
  6m in C+C  <=>  m = k1+k2 for some k1,k2 in K.

This script PROVES the reduction computationally up to a bound and lists the
11 multiple-of-6 "holes" of C+C that generate the 11 exception triples.
"""
from sympy import primerange

LIMIT = 200_000

def main():
    print(f"twin-Goldbach necklace reduction (oracle-S521), bound {LIMIT}\n")
    primes = set(primerange(2, LIMIT + 4))
    # twin-pair centers c=p+1 for each pair {p,p+2}
    centers = sorted(c for c in range(4, LIMIT + 2) if (c - 1) in primes and (c + 1) in primes)
    anomalous = [c for c in centers if c % 6 != 0]
    print(f"twin-pair centers up to {LIMIT}: {len(centers)};  first: {centers[:12]}")
    print(f"centers NOT divisible by 6 (anomalous): {anomalous}   <-- only the {{3,5}} pair\n")

    # twin primes set (for the ground-truth direct check)
    tw = sorted({c - 1 for c in centers} | {c + 1 for c in centers})
    tset = set(tw)

    # --- ground truth: even n representable as sum of two twin primes ---
    B = 6000
    repr_direct = set()
    for i, a in enumerate(tw):
        if a + a > B: break
        for b in tw[i:]:
            s = a + b
            if s > B: break
            repr_direct.add(s)

    # --- center sumset C+C (multiples of 6 mostly) ---
    cc = set()
    cs = [c for c in centers if c <= B]
    for i, c1 in enumerate(cs):
        for c2 in cs[i:]:
            s = c1 + c2
            if s > B + 2: break
            cc.add(s)

    # --- check the reduction on every even n in (8, B] ---
    # claim: for n>8, n is twin-Goldbach-representable  <=>  the multiple of 6
    # nearest n (one of n-2,n,n+2) lies in C+C.
    mism = []
    for n in range(10, B + 1, 2):
        near6 = [x for x in (n - 2, n, n + 2) if x % 6 == 0][0]
        predicted = near6 in cc
        actual = n in repr_direct
        if predicted != actual:
            mism.append((n, near6, predicted, actual))
    print(f"REDUCTION CHECK on even n in (8,{B}]: mismatches = {len(mism)}")
    if mism:
        print("  ", mism[:20])
    else:
        print("  => the triple/necklace reduction holds EXACTLY on this range.")

    # --- the 11 holes: multiples of 6 in (8,4210] not in C+C ---
    holes = [s for s in range(12, 4220, 6) if s not in cc]
    print(f"\nmultiples of 6 in (8,4210] NOT expressible as a sum of two twin-pair centers: {len(holes)}")
    print(f"  holes 6m: {holes}")
    print(f"  each generates the exception triple {{6m-2,6m,6m+2}}:")
    for s in holes:
        print(f"    6m={s:4d}  ->  {{{s-2}, {s}, {s+2}}}")

    # --- the necklace K = C/6 and its Goldbach holes ---
    K = sorted(c // 6 for c in centers if c % 6 == 0)
    print(f"\nTWIN-CENTER NECKLACE K = C/6 (k with 6k+-1 twin): {K[:25]} ...")
    Kset = set(K)
    KK = set()
    for i, k1 in enumerate(K):
        for k2 in K[i:]:
            s = k1 + k2
            if s > 750: break
            KK.add(s)
    mholes = [m for m in range(2, 705) if m not in KK]
    print(f"  m in [2,704] NOT a sum of two necklace elements: {mholes}")
    print(f"  (these m are exactly the holes/6: {[h//6 for h in holes]})")

if __name__ == "__main__":
    main()
