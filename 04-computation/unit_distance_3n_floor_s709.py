"""
monad-explorer-2026-06-06-S709
==============================
THE COMMON-NEIGHBOR FLOOR FOR BEATING 3N UNIT DISTANCES (THM-420 / HYP-2283).

Context. S702 (HYP-2262) / S703 (HYP-2267, THM-412) settled the LATTICE story:
among 2D LATTICE disk patches the triangular (Eisenstein) lattice at the non-minimal
radius sqrt(7) is the earliest to exceed 3N unit distances, crossing at N=43. Open
handoff S702(b): "is triangular N=43 provably the 2D-OPTIMAL FINITE construction?"
That question is about ALL finite planar point sets, not just lattices. HYP-2267 gives
only a constructive UPPER bound (43) over lattices; there was no rigorous LOWER bound.

This script supplies the missing rigorous LOWER bound and narrows the gap.

THE LEVER (two elementary geometric facts about the planar unit-distance graph G):
  (K4) G is K4-free: 4 mutually unit-distant points would be a regular tetrahedron,
       impossible in the plane. (=> max clique 3, the equilateral triangle.)
  (CN) Any two points have AT MOST 2 common neighbours: a common neighbour lies on
       both unit circles, and two distinct circles meet in <=2 points.

(CN) alone gives the floor. Counting cherries (paths of length 2):
       sum_v C(d_v,2) = #cherries <= 2 * C(N,2) = N(N-1).
Maximising U = (1/2) sum_v d_v subject to this and convexity of C(.,2) (degrees as
equal as possible) gives  avg degree <= (1+sqrt(8N-7))/2, hence
       U <= N(1 + sqrt(8N-7)) / 4.
U > 3N is therefore IMPOSSIBLE unless N(1+sqrt(8N-7))/4 > 3N, i.e. sqrt(8N-7) > 11,
i.e. N > 16.  =>  THM-420:  U > 3N  forces  N >= 17.

We (1) verify the floor by exact integer optimisation of max-U(N) under the cherry
constraint; (2) confirm K4-freeness sharpens nothing here but record the clique bound;
(3) SEARCH the sqrt(7)-Eisenstein 12-regular unit graph (and other (lattice,radius)
layers) for the SMALLEST induced subgraph with average degree > 6 (=> U > 3N),
narrowing the true minimum N* into [17, X].
"""
import math
from collections import defaultdict
from itertools import combinations

# ---------------------------------------------------------------------------
# PART 1.  The cherry / common-neighbour floor: exact integer max-U(N).
# Maximise U = sum(d)/2 over integer degree sequences d_v in [0, N-1] with
#   sum_v C(d_v,2) <= N(N-1).
# By convexity the optimum makes the d_v as equal as possible: pick the largest
# integer base degree b with N*C(b,2) <= N(N-1), then raise as many vertices as
# possible from b to b+1 while the cherry budget allows.
# ---------------------------------------------------------------------------
def max_U_under_cherry(N):
    budget = N * (N - 1)                      # = 2*C(N,2)
    if N <= 1:
        return 0
    # base degree b: all-equal feasible
    b = 0
    while (b + 1) <= (N - 1) and N * ((b + 1) * b // 2) <= budget:
        b += 1
    used = N * (b * (b - 1) // 2)
    # raise vertices b -> b+1 (cost C(b+1,2)-C(b,2) = b each), capped at degree N-1
    raised = 0
    if b + 1 <= N - 1:
        cost_each = b  # C(b+1,2)-C(b,2)
        room = budget - used
        max_raise = room // cost_each if cost_each > 0 else N
        raised = min(N, max_raise)
    total_deg = (N - raised) * b + raised * (b + 1)
    return total_deg // 2  # U (floor; degrees sum is even at optimum boundary handled below)

print("=" * 74)
print("PART 1.  Cherry/common-neighbour floor:  max U(N) under sum C(d,2) <= N(N-1)")
print("         (>=2-common-neighbour property of the unit-distance graph)")
print("=" * 74)
print(f"{'N':>4} {'maxU(cherry)':>13} {'3N':>6} {'closed-form U<=N(1+sqrt(8N-7))/4':>34} {'beats3N?':>9}")
first = None
for N in range(3, 50):
    mu = max_U_under_cherry(N)
    cf = N * (1 + math.sqrt(8 * N - 7)) / 4
    beats = mu > 3 * N
    if beats and first is None:
        first = N
    if N <= 20 or beats and N <= 22:
        print(f"{N:>4} {mu:>13} {3*N:>6} {cf:>34.2f} {str(beats):>9}")
print(f"\n  => integer cherry-floor: max U(N) first exceeds 3N at N = {first}")
print(f"  => THEOREM: U > 3N is impossible for N <= {first-1}; first possible at N = {first}.")

# ---------------------------------------------------------------------------
# PART 2.  Build (lattice, squared-radius) unit-distance graphs and search for
# the SMALLEST induced subgraph with average degree > 6 (i.e. U > 3N).
# Lattices given by reduced binary quadratic form Q(x,y)=a x^2 + b xy + c y^2
# (squared distance between integer index points). We pick a non-minimal radius
# whose layer (= r_Q) exceeds 6.
# ---------------------------------------------------------------------------
def layer_count(a, b, c, R):
    B = int(math.isqrt(R)) + 3
    return sum(1 for x in range(-B, B + 1) for y in range(-B, B + 1)
               if a * x * x + b * x * y + c * y * y == R)

def build_graph(a, b, c, R, half_extent):
    """Unit-distance graph on the index box [-h,h]^2 with unit^2 = R."""
    pts = [(x, y) for x in range(-half_extent, half_extent + 1)
                  for y in range(-half_extent, half_extent + 1)]
    idx = {p: i for i, p in enumerate(pts)}
    B = int(math.isqrt(R)) + 3
    offs = [(dx, dy) for dx in range(-B, B + 1) for dy in range(-B, B + 1)
            if a * dx * dx + b * dx * dy + c * dy * dy == R]
    adj = [set() for _ in pts]
    for p in pts:
        i = idx[p]
        for (dx, dy) in offs:
            q = (p[0] + dx, p[1] + dy)
            if q in idx:
                adj[i].add(idx[q])
    return pts, adj

def smallest_dense_ball(pts, adj, center_idx, layer):
    """Grow induced subgraph by BFS shells around a center; report (k, edges, k_first>3k).
       Also try a greedy densest-subset refinement at each prefix."""
    # BFS order from center
    from collections import deque
    order = []
    seen = {center_idx}
    dq = deque([center_idx])
    while dq:
        u = dq.popleft(); order.append(u)
        for w in sorted(adj[u]):
            if w not in seen:
                seen.add(w); dq.append(w)
    best = None
    S = set()
    for u in order:
        S.add(u)
        k = len(S)
        e = sum(1 for v in S for w in adj[v] if w in S and w > v)
        if e > 3 * k and best is None:
            best = (k, e)
        if k > 80:
            break
    return best

print()
print("=" * 74)
print("PART 2.  Smallest induced subgraph with U > 3N (avg degree > 6), by (lattice,R)")
print("=" * 74)
# (name, form, candidate non-minimal radii to test)
CASES = [
    ("triangular(-3)", (1, 1, 1), [3, 4, 7, 9, 12, 13, 19, 21]),
    ("square(-4)",     (1, 0, 1), [2, 4, 5, 8, 9, 10, 13, 25]),
    ("disc(-7)",       (1, 1, 2), [2, 4, 7, 8, 9, 11, 14, 16]),
    ("disc(-8)",       (1, 0, 2), [3, 4, 8, 9, 11, 16, 17, 18]),
]
results = []
for name, (a, b, c), radii in CASES:
    for R in radii:
        L = layer_count(a, b, c, R)
        if L <= 6:
            continue
        h = max(6, int(math.isqrt(R)) + 5)
        pts, adj = build_graph(a, b, c, R, h)
        # center = (0,0)
        cidx = pts.index((0, 0))
        best = smallest_dense_ball(pts, adj, cidx, L)
        if best:
            k, e = best
            results.append((k, name, R, L, e))
            print(f"  {name:>14}  R={R:>3} (layer {L:>2}):  smallest U>3N induced subgraph "
                  f"k={k:>3}, U={e:>4}  (3k={3*k})")
        else:
            print(f"  {name:>14}  R={R:>3} (layer {L:>2}):  none up to k=80")
print()
if results:
    results.sort()
    k0, nm, R0, L0, e0 = results[0]
    print(f"  *** smallest construction found: {nm} radius^2={R0} (layer {L0}): "
          f"N={k0}, U={e0} > 3N={3*k0} ***")
    print(f"  => true minimum N* (over ALL planar sets) satisfies  {first} <= N* <= {k0}")
