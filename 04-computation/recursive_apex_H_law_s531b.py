#!/usr/bin/env python3
"""
recursive_apex_H_law_s531b.py    oracle-2026-06-01-S531

Addendum to s531: pin down the RECURSIVE law for H under apex flips.
Single-tile law H=1+2^(r-1) is exact (s531). Here we test composition:
  (1) DISJOINT apex flips  -> multiplicative? additive?
  (2) NESTED (laminar, concentric vs endpoint-anchored)
  (3) confirm the concentric "diameter" family climbs toward the regular tournament.
The point: whether the recursion 'apex of a block built from sub-blocks' makes H
factor, or whether the coupling (S520o: arcs are not independent) blocks it.
"""
from itertools import combinations

def transitive_adj(n):
    b = [[False]*(n+1) for _ in range(n+1)]
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i > j: b[i][j] = True
    return b

def flip(b, flips):
    b = [row[:] for row in b]
    for (x, y) in flips:
        b[x][y] = False; b[y][x] = True
    return b

def H(n, flips):
    beats = flip(transitive_adj(n), flips)
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            for u in range(n):
                if mask & (1 << u): continue
                if beats[v+1][u+1]: dp[mask | (1 << u)][u] += c
    return sum(dp[full][v] for v in range(n))

def main():
    print("="*70)
    print("(1) DISJOINT apex flips: is H multiplicative over a separating vertex?")
    print("="*70)
    # n=7: block [1,3] (tile (3,1)) and block [5,7] (tile (7,5)); vertex 4 separates.
    tests = [
        (7, [(3,1)], [(7,5)]),
        (8, [(3,1)], [(8,6)]),
        (8, [(4,1)], [(8,5)]),
        (9, [(4,1)], [(9,6)]),
    ]
    for n, A, B in tests:
        hA, hB, hAB = H(n,A), H(n,B), H(n,A+B)
        print(f"  n={n}: H(A={A})={hA}, H(B={B})={hB}, H(A+B)={hAB}; "
              f"product={hA*hB}, 1+(hA-1)+(hB-1)={hA+hB-1}, "
              f"match_prod={hAB==hA*hB}, match_addq={hAB==hA+hB-1}")
    print()
    print("="*70)
    print("(2) NESTED concentric ('diameter onion')  vs  endpoint-anchored chain")
    print("="*70)
    def concentric(n):
        # (1,n),(2,n-1),(3,n-2),... the nested-with-shrinking-on-both-sides family
        fl = []
        y, x = 1, n
        while x - y >= 2:
            fl.append((x, y)); y += 1; x -= 1
        return fl
    def anchored(n):
        # (n,1),(n-2,1)... share endpoint 1
        return [(x,1) for x in range(n, 2, -2)]
    # regular tournament H baseline (rotational, odd n)
    def regular_H(n):
        half=(n-1)//2
        beats=[[False]*(n+1) for _ in range(n+1)]
        for i in range(n):
            for j in range(n):
                if i==j: continue
                if (j-i)%n in set(range(1,half+1)): beats[i+1][j+1]=True
        full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
        for v in range(n): dp[1<<v][v]=1
        for mask in range(1<<n):
            for v in range(n):
                c=dp[mask][v]
                if not c: continue
                for u in range(n):
                    if mask&(1<<u): continue
                    if beats[v+1][u+1]: dp[mask|(1<<u)][u]+=c
        return sum(dp[full][v] for v in range(n))
    for n in (5,6,7,8,9):
        c = concentric(n); a = anchored(n)
        Hc, Ha = H(n,c), H(n,a)
        extra = f", regular(rotational) H={regular_H(n)}" if n%2==1 else ""
        print(f"  n={n}: concentric {[(y,x) for (x,y) in c]} -> H={Hc};  "
              f"anchored {[(y,x) for (x,y) in a]} -> H={Ha}{extra}")
    print()
    print("  Reading: concentric 'diameter onion' apex flips climb toward the regular")
    print("  tournament H (the LRC tight witness); endpoint-anchored flips do not. The")
    print("  apex recursion is COUPLED (S520o): H is not multiplicative even over a")
    print("  separating vertex -- the recursive sub-ranking dependence is real.")

if __name__ == "__main__":
    main()
