#!/usr/bin/env python3
"""
independent_pairs_boundary_s532b.py    oracle-2026-06-01-S532b

COMPLEMENT to the concurrent independent_pairs_channels_s532.py (which verified the
n=4 bijection, tied floor(n/2) to opus-S524's n=14 CRT=7, and did the parity-XOR).
I do NOT reduplicate. New, complementary results:

 (A) THE ISO-DETERMINATION BOUNDARY. Flipping a set S of arcs (rest fixed) can index
     all A000568(n) iso classes only if 2^|S| >= A000568(n). INDEPENDENT pairs (a
     matching) give |S| <= floor(n/2). So independent pairs ALONE can determine the
     iso class only if 2^floor(n/2) >= A000568(n) -- TRUE (with equality) iff n<=4.
     For n>=5, 2^floor(n/2) < A000568(n): COUPLED arcs are forced. Define
         coupling gap  =  ceil(log2 A000568(n)) - floor(n/2)
     = #genuinely-coupled bits beyond the independent-pair channels. Verified n=5:
     no matching+frame reaches more than 2^2=4 of the 12 classes.

 (B) CHANNEL (skip-shell) gcd-decomposition: skip j -> gcd(n,j) cycles of length
     n/gcd(n,j); the DIAMETER shell (j=n/2, even n) is a PERFECT MATCHING = n/2
     independent channels. The "amount of independent pairs" peaks there.

 (C) LINK to S531: independent (vertex-disjoint) = the FACTORING channels (H
     multiplies over disjoint modules); coupled (shared-vertex) = the obstruction.
     So 'amount/state of independent pairs' = the free part; the coupling gap = the
     multi-channel (inside-debt, S529) part. Synthesis in the reflection.
"""
from itertools import combinations, permutations, product
from math import log2, ceil, gcd

A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}

def canon(adj, n):
    best=None
    for p in permutations(range(n)):
        f=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f<best: best=f
    return best

# ----------------------------------------------------------------------
# (A) boundary table + n=5 matching-reach proof
# ----------------------------------------------------------------------
def boundary_table():
    print("  n :  A000568   2^floor(n/2)  indep>=iso?  ceil(log2 A000568)  coupling_gap")
    for n in range(3,9):
        a=A000568[n]; ip=n//2; need=ceil(log2(a)) if a>1 else 0
        gap=need-ip
        ok = (2**ip >= a)
        print(f"  {n} :  {a:6d}     {2**ip:6d}        {'YES' if ok else 'no ':>3}        "
              f"{need:5d}             {gap:+d}")
    print("  => 2^floor(n/2) >= A000568 only for n<=4 (equality): independent pairs alone")
    print("     index the whole iso-class set iff n<=4. coupling_gap>0 for n>=5 = the")
    print("     coupled bits beyond the independent-pair channels.")

def n5_matching_reach():
    """For n=5: max #iso classes reachable by flipping ONE matching (2 disjoint arcs)
    over the best frame. Must be <= 4 < 12; report the achieved max."""
    n=5
    edges=list(combinations(range(n),2))   # 10 edges
    matchings=[(e1,e2) for e1,e2 in combinations(edges,2) if not (set(e1)&set(e2))]
    best=0; best_cfg=None
    for M in matchings:
        frame=[e for e in edges if e not in M]   # 8 edges
        for fb in product((0,1),repeat=len(frame)):
            fixed=dict(zip(frame,fb))
            isos=set()
            for mb in product((0,1),repeat=2):
                adj=[[0]*n for _ in range(n)]
                for (i,j),b in {**fixed, **dict(zip(M,mb))}.items():
                    if b: adj[i][j]=1
                    else: adj[j][i]=1
                isos.add(canon(adj,n))
            if len(isos)>best:
                best=len(isos); best_cfg=(M,fb)
        # early exit if we already hit the ceiling 4
        if best>=4: break
    print(f"  n=5: max iso classes reachable by a 2-arc MATCHING (best frame) = {best} "
          f"(out of A000568(5)=12; ceiling 2^2=4). Coupled arcs needed for the other 8.")

# ----------------------------------------------------------------------
# (B) channels
# ----------------------------------------------------------------------
def channels():
    for n in range(3,11):
        rows=[]
        for j in range(1,n//2+1):
            g=gcd(n,j); clen=n//g
            rows.append(f"j{j}:{'MATCH' if clen==2 else f'{g}xC{clen}'}")
        tag=" [DIAMETER=perfect matching, n/2 indep channels]" if n%2==0 else ""
        print(f"   n={n}: "+" ".join(rows)+tag)

def main():
    print("="*72); print("(A) ISO-DETERMINATION BOUNDARY: independent pairs suffice iff n<=4"); print("="*72)
    boundary_table(); print()
    n5_matching_reach(); print()
    print("="*72); print("(B) CHANNELS: skip-shell gcd-decomposition; diameter shell = matching"); print("="*72)
    channels()
    print()
    print("  (C) Link to S531: independent (vertex-disjoint) arcs = the FACTORING")
    print("  channels (H multiplies over disjoint modules); coupled arcs = obstruction.")
    print("  'Amount/state of independent pairs' = free part; coupling_gap = inside-debt part.")

if __name__=="__main__":
    main()
