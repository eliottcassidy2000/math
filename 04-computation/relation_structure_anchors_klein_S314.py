#!/usr/bin/env python3
"""
relation_structure_anchors_klein_S314.py — klein-2026-07-16-S314 (cont.8)
Anchors for the relation-as-object reflection:
 (1) COMPOSITION DEFECT: c3 = #{(u,w) pairs realizing R.R ∩ R^op} / counted per cyclic
     triple — precisely: the number of cyclic triples equals the number of ordered pairs
     (u,w) with w->u and exists v: u->v->w, divided by ... verified exactly: each cyclic
     triangle contributes its 3 "chords": |R∘R ∩ R^op| counted as ordered pairs with
     multiplicity = 3*c3 when counted with witnesses; as a SET of ordered pairs it needs
     multiplicity — we verify the exact identity  sum_{(u,w) in R^op} #{v: uRv, vRw} = 3*c3? 
     no: each cyclic triangle u->v->w->u gives the pairs (u,w),(v,u),(w,v) each with >=1
     witness; total WITNESS COUNT = # directed paths of length 2 closed by a reverse arc
     = 3*c3 exactly (each cyclic triple has exactly 3 such configurations).
 (2) FRAME DICTIONARY: under reversal T -> T^op: H invariant, c3 invariant, x invariant
     (acceleration-like); score sequence maps s -> n-1-s (velocity-like, covariant).
"""
import itertools, random
from math import comb
def paths2_closed(T, n):
    tot = 0
    for u in range(n):
        for w in range(n):
            if u != w and T[w][u]:
                tot += sum(1 for v in range(n) if v != u and v != w and T[u][v] and T[v][w])
    return tot
def c3(T, n):
    s = [sum(T[u]) for u in range(n)]
    return comb(n, 3) - sum(comb(x, 2) for x in s)
def H(T, n):
    # count Hamiltonian paths via DP
    full = 1 << n
    dp = [[0]*n for _ in range(full)]
    for v in range(n): dp[1 << v][v] = 1
    for m in range(full):
        for v in range(n):
            if dp[m][v]:
                for u in range(n):
                    if not (m >> u) & 1 and T[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
    return sum(dp[full-1][v] for v in range(n))
rng = random.Random(8)
ok1 = ok2 = True
for trial in range(60):
    n = rng.randint(4, 7)
    T = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            T[u][v] = rng.randint(0, 1); T[v][u] = 1 - T[u][v]
    if paths2_closed(T, n) != 3 * c3(T, n): ok1 = False
    Top = [[T[v][u] for v in range(n)] for u in range(n)]
    s = sorted(sum(T[u][v] for v in range(n) if v != u) for u in range(n))
    sop = sorted(sum(Top[u][v] for v in range(n) if v != u) for u in range(n))
    if not (H(T, n) == H(Top, n) and c3(T, n) == c3(Top, n)
            and sop == sorted(n - 1 - x for x in s)): ok2 = False
print("PASS (1) composition defect: #(2-paths closed by a reverse arc) = 3*c3 EXACTLY"
      if ok1 else "FAIL (1)")
print("PASS (2) frame dictionary: H, c3 reversal-INVARIANT; scores reversal-COVARIANT "
      "(s -> n-1-s)" if ok2 else "FAIL (2)")
