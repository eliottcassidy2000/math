#!/usr/bin/env python3
"""
lrc_p_adic_tree_s521.py   claudebox-2026-06-01-S521

The p-adic tree view of LRC (reflection: 07-reflections/lrc-p-adic-tree-s521.md).

TIME TREE: p-ary subdivision of [0,1); node [a/p^k,(a+1)/p^k] is "lonely-live" if it
contains a lonely time. Live subtree growth: FAT (->p, dim 1) = non-tight; the tight
extremizers have lonely point t=1/n = a non-dyadic infinite path (dim 0, captured by
no finite node).
SPEED TREE: doubling v_j=2v_i = 2-adic tree edge (benign); sum-relation v_i+v_j=v_k =
cross-tree additive (the equidistribution obstruction). First-even bridge n=2*odd =
the 2-adic valuation split.
"""
from fractions import Fraction as F
def dist(x):
    x = x % 1; return min(x, 1 - x)
def node_live(v, a, pk, n, depth=6):
    lo = F(a, pk); hi = F(a+1, pk)
    for j in range(0, 2**depth + 1):
        t = lo + (hi-lo)*F(j, 2**depth)
        if all(dist(F(vi)*t) >= F(1, n) for vi in v): return True
    return False
def tree_profile(v, p, n, K):
    return [sum(1 for a in range(p**k) if node_live(v, a, p**k, n)) for k in range(1, K+1)]
def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def main():
    print("TIME TREE (p=2): # lonely-live nodes per level. Fat(grows)=non-tight; 0=tight (path 1/n).")
    for v in [(1,2,4,7),(1,2,3,4),(1,3,4,5,9),(2,3,5,7,11)]:
        n = len(v)+1
        prof = tree_profile(list(v), 2, n, 6)
        tag = "non-tight (fat, dim 1)" if prof[-1] > prof[0] else "TIGHT (dim 0, path 1/n)"
        print(f"  v={v} (n={n}): live lvl1..6 = {prof}   -> {tag}")
    print("\nSPEED TREE: doublings (2-adic edges, benign) vs sum-relations (cross-tree, obstruction):")
    for v in [(1,2,3,4),(1,2,4,7),(2,3,5,7,11)]:
        doublings = [(a, b) for a in v for b in v if b == 2*a]
        sumrels = [(a, b, c) for a in v for b in v for c in v if a < b and a+b == c]
        vals = {x: v2(x) for x in v}
        print(f"  v={v}: 2-adic vals={vals}; doublings={doublings}; sum-relations={sumrels}"
              f"  => {'RESONANT' if sumrels else 'non-resonant'}")
    print("\n  Doublings are harmless 2-adic tree edges; the cross-tree additive triples v_i+v_j=v_k")
    print("  are the equidistribution obstruction (Galois-Weil). First-even bridge n=2*odd = 2-adic split.")

if __name__ == "__main__":
    main()
