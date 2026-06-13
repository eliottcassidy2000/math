"""
subdivision_lacunary_dyadic_kps2.py  (kind-pasteur-2026-06-09-S2, Branch II supplement)

THE LACUNARY HARD-CORE EXERCISE from the Erdos-Gyarfas reframe (HYP-2358/THM-456):
for a subdivision of a cubic multigraph H with edge lengths l_e, the cycle spectrum
is exactly { sum of l_e over E(C) : C a cycle of H }.  E-G for subdivisions of H
would say "every positive length assignment gives some cycle with dyadic total".
This script shows that is FALSE for EVERY topology, and finds the smallest witnesses.

(1) UNIVERSAL AVOIDANCE (all-lengths-3 rule): assigning every edge length 3 makes
    every cycle sum divisible by 3, hence never a power of 2.  So E-G fails on a
    subdivision of EVERY cubic multigraph -- but subdivision vertices have degree 2,
    so these graphs never have min degree 3: the delta>=3 hypothesis carries ALL the
    weight of Erdos-Gyarfas.  Verified by explicit spectrum DP on subdivided graphs.

(2) SMALLEST WITNESSES: exhaustive search over the loopless cubic multigraphs on
    n_H <= 4 vertices (theta3 = triple edge on 2 vertices; K4; the "ladder" L4 =
    4-cycle with two opposite edges doubled) and all length assignments, minimizing
    the subdivision order n = n_H + sum(l_e - 1), requiring the subdivision SIMPLE
    (no two parallel edges both of length 1) and all cycle sums avoiding
    {4, 8, 16, 32, ...}.

All claimed spectra are independently re-verified with the bitmask cycle-spectrum DP
(cycle anchored at min vertex -- correct for cycles; cf. MISTAKE-068) on the actual
subdivided simple graph.

Output -> 05-knowledge/results/subdivision_lacunary_dyadic_kps2.out
"""
import os, sys, time
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "05-knowledge", "results",
                   "subdivision_lacunary_dyadic_kps2.out")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
    def write(self, s):
        sys.__stdout__.write(s); self.f.write(s)
    def flush(self):
        sys.__stdout__.flush(); self.f.flush()

sys.stdout = Tee(OUT)
t0 = time.time()

# ---------------- the three loopless cubic multigraphs on <= 4 vertices --------
# (n_H = 3 impossible: odd degree sum.  Enumeration is classical; we hardcode.)
# Each: (name, n_H, edge list as (u,v), list of cycles as edge-index tuples)
THETA3 = ("theta3 (2 vertices, triple edge)", 2,
          [(0, 1), (0, 1), (0, 1)],
          [(0, 1), (0, 2), (1, 2)])
K4 = ("K4", 4,
      [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
      [(0, 1, 3), (0, 2, 4), (1, 2, 5), (3, 4, 5),       # 4 triangles
       (0, 3, 5, 2), (0, 4, 5, 1), (1, 3, 4, 2)])        # 3 quads
LADDER = ("L4 ladder (C4 with two opposite edges doubled)", 4,
          [(0, 1), (0, 1), (2, 3), (2, 3), (0, 2), (1, 3)],
          [(0, 1), (2, 3),                                # two 2-cycles
           (0, 4, 2, 5), (0, 4, 3, 5), (1, 4, 2, 5), (1, 4, 3, 5)])
MULTIGRAPHS = [THETA3, K4, LADDER]

def parallel_classes(edges):
    cls = {}
    for i, (u, v) in enumerate(edges):
        cls.setdefault((min(u, v), max(u, v)), []).append(i)
    return [v for v in cls.values() if len(v) > 1]

def is_simple_assignment(edges, lengths):
    for grp in parallel_classes(edges):
        ones = [i for i in grp if lengths[i] == 1]
        if len(ones) > 1:
            return False
    return True

def cycle_sums(cycles, lengths):
    return sorted(set(sum(lengths[i] for i in c) for c in cycles))

def dyadic_hit(sums):
    return sorted(s for s in sums if s >= 4 and (s & (s - 1)) == 0)

# ---------------- subdivision builder + spectrum DP ----------------------------
def subdivide(n_H, edges, lengths):
    """return adjacency bitmasks of the subdivision (edge i -> path of lengths[i])"""
    nxt = n_H
    adj_pairs = []
    for (u, v), L in zip(edges, lengths):
        prev = u
        for _ in range(L - 1):
            adj_pairs.append((prev, nxt)); prev = nxt; nxt += 1
        adj_pairs.append((prev, v))
    n = nxt
    adj = [0] * n
    for a, b in adj_pairs:
        if a == b or (adj[a] >> b) & 1:
            raise ValueError("subdivision not simple")
        adj[a] |= 1 << b; adj[b] |= 1 << a
    return adj

def cycle_spectrum(adj):
    """exact cycle-length set; subset DP anchored at the cycle's min vertex"""
    n = len(adj)
    ep = [0] * (1 << n)
    for v in range(n):
        ep[1 << v] = 1 << v
    spec = set()
    for S in range(1, 1 << n):
        e = ep[S]
        if not e:
            continue
        m = S & (-S)
        k = S.bit_count()
        if k >= 3 and (e & adj[m.bit_length() - 1]):
            spec.add(k)
        U = 0
        t = e
        while t:
            b = t & (-t); U |= adj[b.bit_length() - 1]; t ^= b
        c = U & ~S & ~((m << 1) - 1)
        while c:
            b = c & (-c); ep[S | b] |= b; c ^= b
    return sorted(spec)

# ---------------- (1) universal all-lengths-3 avoidance ------------------------
print("=" * 78)
print("(1) UNIVERSAL AVOIDANCE: every edge length 3  =>  all cycle sums in 3Z,")
print("    and 3k is never a power of 2.  DP verification on the actual graphs:")
print("=" * 78)
for name, n_H, edges, cycles in MULTIGRAPHS:
    L = [3] * len(edges)
    pred = cycle_sums(cycles, L)
    adj = subdivide(n_H, edges, L)
    spec = cycle_spectrum(adj)
    n = len(adj)
    ok = (spec == pred) and not dyadic_hit(spec)
    print(f"  {name}: subdivision n={n}, predicted sums {pred}, DP spectrum {spec},"
          f" dyadic hits {dyadic_hit(spec)}  -> {'OK' if ok else 'FAIL'}")
    assert ok

# ---------------- (2) smallest simple dyadic-avoiding subdivisions -------------
print()
print("=" * 78)
print("(2) SMALLEST SIMPLE SUBDIVISION of each cubic multigraph with ALL cycle")
print("    sums avoiding {4,8,16,32,...}   (exhaustive over length vectors)")
print("=" * 78)
MAXL = 6
global_best = None
for name, n_H, edges, cycles in MULTIGRAPHS:
    m = len(edges)
    best = None
    n_assign = 0
    for L in product(range(1, MAXL + 1), repeat=m):
        if not is_simple_assignment(edges, L):
            continue
        n_assign += 1
        sums = cycle_sums(cycles, L)
        if dyadic_hit(sums):
            continue
        n = n_H + sum(x - 1 for x in L)
        if best is None or n < best[0]:
            best = (n, L, sums)
    n, L, sums = best
    adj = subdivide(n_H, edges, L)
    spec = cycle_spectrum(adj)
    ok = (spec == sums) and not dyadic_hit(spec)
    deg2 = sum(1 for a in adj if a.bit_count() == 2)
    print(f"  {name}: searched {n_assign} simple assignments (l_e <= {MAXL})")
    print(f"     minimum subdivision order n = {n}, lengths {L}, "
          f"cycle spectrum {sums} (DP re-check {spec}: {'OK' if ok else 'FAIL'})")
    print(f"     degree-2 vertices: {deg2}/{n}  -> min degree 2, NOT an E-G graph")
    assert ok
    if global_best is None or n < global_best[0]:
        global_best = (n, name, L, sums)

n, name, L, sums = global_best
print()
print(f"GLOBAL MINIMUM: {name} with lengths {L}: a SIMPLE graph on {n} vertices")
print(f"  whose full cycle spectrum {sums} avoids every power of 2 -- this is the")
print(f"  generalized theta graph theta{tuple(sorted(L))}.")
print()
print("CONTRAST (THM-456 census): among genuine min-degree-3 graphs, NO graph on")
print("n <= 12 avoids {4,8,16} (exhaustive); Royle-Markstrom: none on n <= 16.")
print("Subdivisions dodge dyadic lengths with 6 vertices; delta>=3 graphs cannot")
print("even at 16. The entire difficulty of Erdos-Gyarfas problem 64 lives in the")
print("degree-2-free gap between these two worlds.")
print()
print(f"total time {time.time()-t0:.1f} s")
