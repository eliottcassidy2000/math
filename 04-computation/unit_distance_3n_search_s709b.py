"""
monad-explorer-2026-06-06-S709b
===============================
NARROWING [17, 43]: search for the smallest planar unit-distance config with U > 3N.

THM-420 (s709): U > 3N forces N >= 17 (>=2-common-neighbour cherry floor).
HYP-2267 (s703): triangular sqrt(7) DISK patch beats 3N at N=43 (lattice minimizer).

Here we search broadly within lattice unit-distance graphs (the only configs known to
beat 3N), over (lattice, radius=layer) AND patch SHAPE (Euclidean disks at several
centers, plus a greedy SHRINK refinement that removes low-surplus boundary vertices),
to test whether ANY config beats 3N below N=43 and to identify the best small shape.

Surplus identity: e(S) - 3|S| = (1/2) * sum_{v in S} (deg_S(v) - 6).  Maximise interior
(deg 12 in layer-12 graphs => +6 each) over boundary (negative). A round patch is
isoperimetrically near-optimal; we also test higher layers (18, 24) where interior
surplus per vertex is larger (+12, +18) against a thicker boundary shell.
"""
import math
from collections import defaultdict

def layer_count(a, b, c, R):
    B = int(math.isqrt(R)) + 3
    return sum(1 for x in range(-B, B + 1) for y in range(-B, B + 1)
               if a * x * x + b * x * y + c * y * y == R)

def offsets(a, b, c, R):
    B = int(math.isqrt(R)) + 3
    return [(dx, dy) for dx in range(-B, B + 1) for dy in range(-B, B + 1)
            if a * dx * dx + b * dx * dy + c * dy * dy == R]

def disk_smallest(a, b, c, R, centers, half):
    """Over candidate (fractional) centers, order box points by squared Euclidean
       distance to the center, grow the prefix, return smallest prefix with U>3N.
       Squared distance uses the true metric of the form:
         ||(x,y)||^2 = a x^2 + b xy + c y^2  (the form IS the Gram quadratic form)."""
    offs = offsets(a, b, c, R)
    box = [(x, y) for x in range(-half, half + 1) for y in range(-half, half + 1)]
    idx = {p: i for i, p in enumerate(box)}
    adj = [[] for _ in box]
    for p in box:
        i = idx[p]
        for (dx, dy) in offs:
            q = (p[0] + dx, p[1] + dy)
            j = idx.get(q)
            if j is not None:
                adj[i].append(j)
    best = None
    for (cx, cy) in centers:
        def d2(p):
            x = p[0] - cx; y = p[1] - cy
            return a * x * x + b * x * y + c * y * y
        order = sorted(box, key=d2)
        S = set(); e = 0
        for p in order:
            i = idx[p]
            S.add(i)
            e += sum(1 for j in adj[i] if j in S)  # j already in S
            k = len(S)
            if e > 3 * k:
                if best is None or k < best[0]:
                    best = (k, e, (cx, cy))
                break
    return best, idx, adj

def shrink(idx, adj, seed_pts_indices):
    """Greedy shrink: from a vertex set with U>3N, repeatedly delete the vertex whose
       removal keeps U>3(N-1) and most increases surplus density, until no deletion
       preserves U>3N. Returns the smallest set reached (k, e)."""
    S = set(seed_pts_indices)
    deg = {v: sum(1 for w in adj[v] if w in S) for v in S}
    e = sum(deg.values()) // 2
    improved = True
    while improved and len(S) > 3:
        improved = False
        # candidate: delete v if e - deg(v) > 3*(|S|-1)
        best_v = None; best_newsurplus = None
        k = len(S)
        for v in S:
            ne = e - deg[v]; nk = k - 1
            if ne > 3 * nk:
                surplus = ne - 3 * nk
                if best_newsurplus is None or surplus > best_newsurplus:
                    best_newsurplus = surplus; best_v = v
        if best_v is not None:
            v = best_v
            for w in adj[v]:
                if w in S:
                    deg[w] -= 1
            e -= deg[v]
            S.discard(v); deg.pop(v)
            improved = True
    return len(S), e

print("=" * 78)
print("Smallest unit-distance config with U > 3N: disk-ordering + shrink refinement")
print("=" * 78)
# Triangular form x^2+xy+y^2. Natural centers: lattice point (0,0); deep hole
# (1/3,1/3) & (2/3,2/3); edge midpoint (1/2,0).
TRI_CENTERS = [(0, 0), (1/3, 1/3), (2/3, 2/3), (1/2, 0), (1/2, 1/2)]
SQ_CENTERS = [(0, 0), (1/2, 1/2), (1/2, 0)]
CASES = [
    ("triangular(-3)", (1, 1, 1), [7, 13, 19, 21, 28, 31, 37, 39, 43, 49, 61, 63, 67, 73, 79, 84, 91, 93, 97], TRI_CENTERS),
    ("square(-4)",     (1, 0, 1), [5, 10, 13, 25, 50, 65, 85], SQ_CENTERS),
]
overall = None
for name, (a, b, c), radii, centers in CASES:
    print(f"\n--- {name} ---")
    for R in radii:
        L = layer_count(a, b, c, R)
        if L <= 6:
            continue
        half = int(math.isqrt(R)) + 6
        best, idx, adj = disk_smallest(a, b, c, R, centers, half)
        if best is None:
            print(f"  R={R:>3} layer {L:>2}: no disk beats 3N within box")
            continue
        k0, e0, ctr = best
        # rebuild the seed disk set at the winning center and shrink
        offs = offsets(a, b, c, R)
        box = list(idx.keys())
        cx, cy = ctr
        def d2(p):
            x = p[0] - cx; y = p[1] - cy
            return a * x * x + b * x * y + c * y * y
        order = sorted(box, key=d2)
        seed = [idx[p] for p in order[:k0]]
        ks, es = shrink(idx, adj, seed)
        tag = f"  R={R:>3} layer {L:>2}: disk N={k0:>3} U={e0:>4} (3N={3*k0}) @c={ctr}"
        tag += f" | shrink-> N={ks:>3} U={es:>4} (3N={3*ks})"
        print(tag)
        cand = min((k0, name, R, L, "disk"), (ks, name, R, L, "shrink"))
        if overall is None or cand[0] < overall[0]:
            overall = cand
print()
print("=" * 78)
if overall:
    print(f"BEST small config: {overall[1]} radius^2={overall[2]} (layer {overall[3]}) "
          f"via {overall[4]}: N = {overall[0]}")
    print(f"=> minimum N* over all planar sets:  17 <= N* <= {overall[0]}   (THM-420 floor 17)")
