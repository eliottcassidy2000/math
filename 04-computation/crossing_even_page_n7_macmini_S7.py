#!/usr/bin/env python3
"""
DECISIVE H2 TEST at n=7 (the canonical even-graph case; n odd => both pages can be
EVEN graphs simultaneously, and E_7 is the first non-perfect even-graph metagraph).
mac-mini-2026-06-20-S7.

Question: among 2-page book drawings of K_7 whose PAGE-0 is an EVEN GRAPH
(every vertex even degree on page 0; since n-1=6 even, page 1 is then also even),
is the MINIMUM crossing number still Guy's Z(7)=9?

If YES through n=7: the crossing optimum is attainable with an even-graph page at
the first nontrivial even-graph case -> a real (not surface) link to E_n's objects.
If NO: H2 fails at n=7; the link is only an n<=6 coincidence.

Method: the even-page constraint = page-0 indicator lies in the CYCLE SPACE of K_7
(GF(2) even subgraphs = cycle space, dimension C(n,2)-n+1 = 21-6 = 15). We enumerate
the cycle space by a basis and minimize crossings over all 2^15 even page-0 graphs.
Also: do the same for the FREE problem via local search to reconfirm Z(7)=9, and
report whether the even-constrained optimum's page-0 is one of the E_7 even-graph
classes (by degree sequence, a coarse class label).
"""

from itertools import combinations
from math import comb
import random

def guy(n):
    return (n//2)*((n-1)//2)*((n-2)//2)*((n-3)//2)//4

def edges_of(n):
    return [(i, j) for i in range(n) for j in range(i+1, n)]

def interleave(e, f):
    (a, b), (c, d) = e, f
    return (a < c < b < d) or (c < a < d < b)

def crossing_pairs(n):
    edges = edges_of(n)
    cps = []
    for x in range(len(edges)):
        for y in range(x+1, len(edges)):
            e, f = edges[x], edges[y]
            if len(set(e) | set(f)) == 4 and interleave(e, f):
                cps.append((x, y))
    return edges, cps

def cycle_space_basis(n, edges):
    """GF(2) cycle space basis of K_n. Spanning tree = star at vertex 0.
       Fundamental cycle for each non-tree edge {i,j} (i,j>=1) = {0-i, i-j, 0-j}."""
    eidx = {e: k for k, e in enumerate(edges)}
    m = len(edges)
    basis = []
    for i in range(1, n):
        for j in range(i+1, n):
            vec = 0
            for e in [(0, i), (i, j), (0, j)]:
                a, b = min(e), max(e)
                vec ^= (1 << eidx[(a, b)])
            basis.append(vec)
    return basis  # length C(n-1,2) = cycle-space dim

def crossings_of_mask(mask, cps):
    return sum(1 for (x, y) in cps if (((mask >> x) & 1) == ((mask >> y) & 1)))

def degseq(n, edges, mask):
    deg = [0]*n
    for k, (i, j) in enumerate(edges):
        if (mask >> k) & 1:
            deg[i] += 1; deg[j] += 1
    return tuple(sorted(deg))

def main():
    n = 7
    Z = guy(n)
    edges, cps = crossing_pairs(n)
    m = len(edges)
    basis = cycle_space_basis(n, edges)
    dim = len(basis)
    print(f"n={n}: edges={m}  Guy Z={Z}  cycle-space dim={dim} (expect {comb(n-1,2)})")
    print(f"Enumerating all 2^{dim}={1<<dim} even page-0 graphs...")

    best = None
    best_mask = None
    deg_of_best = set()
    # iterate all even subgraphs as GF(2) combos of basis
    for combo in range(1 << dim):
        mask = 0
        c = combo
        b = 0
        while c:
            if c & 1:
                mask ^= basis[b]
            c >>= 1; b += 1
        cr = crossings_of_mask(mask, cps)
        if best is None or cr < best:
            best = cr; best_mask = mask
            deg_of_best = {degseq(n, edges, mask)}
        elif cr == best:
            deg_of_best.add(degseq(n, edges, mask))

    print(f"\nMIN crossings over EVEN page-0 drawings: {best}")
    print(f"Guy Z(7) = {Z}   ->  H2 at n=7: {'HOLDS' if best == Z else 'FAILS'}")
    print(f"page-0 degree sequences achieving the even-optimum (a few): "
          f"{sorted(deg_of_best)[:6]}")

    # Reconfirm free optimum via local search
    involved = [[] for _ in range(m)]
    for pi, (x, y) in enumerate(cps):
        involved[x].append(pi); involved[y].append(pi)
    rng = random.Random(999)
    free_best = None
    for r in range(80):
        assign = [rng.randint(0, 1) for _ in range(m)]
        cur = sum(1 for (x, y) in cps if assign[x] == assign[y])
        T = 2.0
        for it in range(15000):
            k = rng.randrange(m)
            delta = 0
            for pi in involved[k]:
                x, y = cps[pi]
                same_now = (assign[x] == assign[y])
                newk = 1 - assign[k]
                same_new = (newk == assign[y]) if x == k else (assign[x] == newk)
                delta += (1 if same_new else 0) - (1 if same_now else 0)
            if delta <= 0 or rng.random() < pow(2.718281828, -delta/max(T,1e-9)):
                assign[k] = 1 - assign[k]; cur += delta
            T *= 0.9995
        if free_best is None or cur < free_best:
            free_best = cur
    print(f"\nFree 2-page optimum (local search): {free_best}  (Z={Z}, match={free_best==Z})")
    print(f"\nSUMMARY: even-constrained optimum {best} vs free optimum {free_best} vs Guy {Z}")
    if best == Z:
        print("=> The crossing optimum IS attainable with an EVEN-GRAPH page at n=7.")
    else:
        print(f"=> Even constraint COSTS {best-Z} crossings at n=7; link is n<=6 only.")

if __name__ == "__main__":
    main()
