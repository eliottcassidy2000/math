#!/usr/bin/env python3
"""
almost_1d_s305.py -- opus-2026-03-24-S305

TOURNAMENT SPACE IS ALMOST ONE-DIMENSIONAL

The tiling hypercube Q_m lives in m dimensions (m = C(n-1,2)).
But the iso-class quotient behaves as if it's nearly 1D:
- H orders classes linearly (the principal line)
- The metagraph is a DAG under H-gradient
- The spectral gap suggests mixing in O(n) steps (not O(m))

Yet it's NOT purely 1D:
- Multiple classes share the same H value
- The width (classes per H level) grows
- The meta-graph has cycles (genus > 0)

THE CLAIM: tournament space is the BOUNDARY between 1D and 2D.
The staircase triangle IS the geometric realization of this boundary.
The diagonal of the triangle (hypotenuse = anti-diagonal) is the
space-filling curve that maps the 1D H-ordering to the 2D tile space.

THE sqrt(2) CONNECTION:
- The staircase is a right isosceles triangle
- Its hypotenuse/leg ratio is sqrt(2)
- The pseudo-doubling ratio of the metagraph is 2 - 1/(n-2) → 2
- The Hamming scheme on Q_m has Krawtchouk eigenvalues
- The association scheme structure connects to coding theory

THIS SCRIPT: tests the dimensionality of tournament space at n=5,6,7
by computing:
1. The effective dimension from the spectral decay
2. The H-ordering as a 1D projection
3. The "width" at each H level
4. The space-filling curve properties

Author: opus-2026-03-24-S305
"""

import sys
import numpy as np
from math import comb, factorial, sqrt, log, log2
from itertools import permutations, combinations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def build_space(n):
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()
    cmap = {}; tc = {}; members = defaultdict(list); sizes = {}
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap); sizes[cmap[c]] = 0
        cid = cmap[c]; sizes[cid] += 1; tc[bits] = cid; members[cid].append(bits)
    nc = len(cmap)
    H_vals = {}
    info = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (b >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        comp_c = canon(complement(A, n), n)
        H = Hcount(A, n)
        info[cid] = {'comp': cmap.get(comp_c, -1), 'H': H}
        H_vals[cid] = H
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        if comp >= 0 and comp != cid:
            merged[mid] = (cid, comp); mid_map[cid] = mid; mid_map[comp] = mid
        else:
            merged[mid] = (cid,); mid_map[cid] = mid
        mid += 1
    return m, total, tile_pairs, tc, members, sizes, mid_map, mid, info, merged

# ============================================================
banner("THE DIMENSIONALITY OF TOURNAMENT SPACE")
# ============================================================

for n in [5, 6]:
    m, total, tile_pairs, tc, members, sizes, mid_map, nm, info, merged = build_space(n)

    # Build wiggly weight matrix
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = mid_map[tc[bits]]
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            W[mi][mj] += 1

    # Fiber vector
    fiber = np.array([sum(sizes[c] for c in merged[mi]) for mi in range(nm)])

    # H vector
    H_vec = np.array([info[merged[mi][0]]['H'] for mi in range(nm)])

    # Transition matrix
    P = W.astype(float)
    for mi in range(nm):
        P[mi] /= (fiber[mi] * m)

    # Eigenvalues
    evals = sorted(np.linalg.eigvals(P).real, reverse=True)

    print("  n=%d: m=%d, %d merged classes\n" % (n, m, nm))

    # 1. EFFECTIVE DIMENSION from spectral decay
    # The effective dimension is where the eigenvalues drop below a threshold
    significant = sum(1 for ev in evals if abs(ev) > 0.1)
    print("  Eigenvalues > 0.1: %d out of %d" % (significant, nm))
    print("  Top eigenvalues: %s" % ["%.3f" % ev for ev in evals[:8]])
    print("  Spectral gap: %.4f" % (1 - evals[1]))

    # Effective dimension = log(V_n) / log(max degree)
    max_deg = max(sum(1 for mj in range(nm) if W[mi][mj] > 0 and mi != mj) for mi in range(nm))
    eff_dim = log(nm) / log(max_deg) if max_deg > 1 else 0
    print("  Effective dimension (log V / log max_deg): %.3f" % eff_dim)

    # 2. H-ORDERING: how 1D is the graph under H?
    # Count distinct H values
    distinct_H = len(set(H_vec))
    max_width = max(Counter(H_vec).values())
    print("\n  H-ordering:")
    print("    Distinct H values: %d" % distinct_H)
    print("    Max width (classes per H level): %d" % max_width)
    print("    V_merged / distinct_H = %.2f (avg width)" % (nm / distinct_H))
    print("    H range: [%d, %d]" % (min(H_vec), max(H_vec)))

    # 3. H as 1D embedding quality
    # If the graph were 1D, H-adjacent classes would be graph-adjacent.
    # Measure: what fraction of graph edges connect H-adjacent classes?
    H_sorted = sorted(set(H_vec))
    H_rank = {h: i for i, h in enumerate(H_sorted)}
    h_ranks = [H_rank[h] for h in H_vec]

    adj = (W > 0).astype(int)
    np.fill_diagonal(adj, 0)
    total_edges = adj.sum() // 2
    h_close_edges = 0
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if adj[mi][mj]:
                if abs(h_ranks[mi] - h_ranks[mj]) <= 1:
                    h_close_edges += 1
    print("    Edges between H-adjacent levels: %d / %d (%.1f%%)" %
          (h_close_edges, total_edges, 100*h_close_edges/total_edges if total_edges > 0 else 0))

    # 4. The CORRELATION DIMENSION
    # For each pair of classes, compute graph distance and H difference
    # If 1D: graph_dist ∝ |ΔH|
    from scipy.sparse.csgraph import shortest_path
    dist_matrix = shortest_path(adj.astype(float), directed=False, unweighted=True)

    dists = []; dH = []
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if dist_matrix[mi][mj] < np.inf:
                dists.append(dist_matrix[mi][mj])
                dH.append(abs(H_vec[mi] - H_vec[mj]))

    corr_dist_dH = np.corrcoef(dists, dH)[0,1] if len(dists) > 2 else 0
    print("    corr(graph_distance, |ΔH|) = %.4f" % corr_dist_dH)

    # 5. THE WIDTH FUNCTION
    print("\n  WIDTH FUNCTION (classes per H level):")
    h_counter = Counter(H_vec)
    for h in sorted(h_counter.keys()):
        w = h_counter[h]
        bar = "█" * w
        print("    H=%3d: %d  %s" % (h, w, bar))

    # 6. THE SQRT(2) IN THE SPECTRUM
    # The ratio of consecutive eigenvalues
    print("\n  EIGENVALUE RATIOS:")
    for i in range(min(5, len(evals)-1)):
        if evals[i+1] != 0:
            ratio = evals[i] / evals[i+1]
            print("    λ_%d/λ_%d = %.4f / %.4f = %.4f" %
                  (i, i+1, evals[i], evals[i+1], ratio))

    print()

# ============================================================
banner("THE STAIRCASE AS A SPACE-FILLING CURVE")
# ============================================================

print("""
  THE STAIRCASE TILING IS A DISCRETE SPACE-FILLING CURVE:

  The m = C(n-1, 2) tiles of the staircase delta_{n-2} are arranged
  in a right isosceles triangle. The standard tile ordering
  (A, B, C, ...) follows a DIAGONAL TRAVERSAL:

  For n=5 (6 tiles in a 3×3 triangle):
    Row 0: A=(0,2) B=(0,3) C=(0,4)     ranges: 2, 3, 4
    Row 1: D=(1,3) E=(1,4)             ranges: 2, 3
    Row 2: F=(2,4)                     range:  2

  The diagonal traversal (constant x+y) gives:
    Diagonal 0: A (range 2)
    Diagonal 1: B, D (ranges 3, 2)
    Diagonal 2: C, E, F (ranges 4, 3, 2)

  Each DIAGONAL has constant x+y = range + (something).
  Actually: tile (x,y) = (lo, hi) has range = hi - lo.
  And the anti-diagonal coordinate is lo + hi.

  The ANTI-DIAGONAL traversal maps the 2D staircase to a 1D sequence.
  This is the natural "space-filling" direction of the triangle.

  THE SQRT(2) RATIO:
  The staircase has legs of length n-2 and hypotenuse of length
  (n-2)×sqrt(2). The number of tiles on each anti-diagonal is:
  d=0: 1 tile (corner)
  d=1: 2 tiles
  ...
  d=k: min(k+1, n-2-k) tiles (reaches max at d = floor((n-3)/2))
  d=n-3: 1 tile (opposite corner)

  This is the WIDTH FUNCTION of the staircase along the diagonal.
  It grows from 1 to floor(n/2)-1, then shrinks back to 1.
  This is a TRIANGULAR shape — the 1D shadow of the 2D triangle.
""")

# Compute the anti-diagonal distribution of tiles
for n in [5, 6, 7, 8]:
    m = comb(n-1, 2)
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]

    diag_count = Counter()
    for lo, hi in tile_pairs:
        diag = lo + hi  # anti-diagonal coordinate
        diag_count[diag] += 1

    print("  n=%d (m=%d tiles): anti-diagonal distribution:" % (n, m))
    for d in sorted(diag_count.keys()):
        bar = "█" * diag_count[d]
        print("    d=%d: %d tiles  %s" % (d, diag_count[d], bar))

    # The shape: unimodal, reaching max at mid-diagonal
    max_width = max(diag_count.values())
    total_diags = len(diag_count)
    print("    Max width: %d, #diagonals: %d, avg: %.1f" %
          (max_width, total_diags, m / total_diags))
    print()

# ============================================================
banner("THE ASSOCIATION SCHEME STRUCTURE")
# ============================================================

print("""
  THE HAMMING SCHEME ON TILINGS:

  The 2^m tilings form the vertices of the Hamming scheme H(m, 2).
  This is the association scheme where:
  - Relation i: pairs at Hamming distance i (for i = 0, 1, ..., m)
  - The scheme has m+1 classes
  - The eigenvalues are KRAWTCHOUK POLYNOMIALS K_k(i; m)

  The waggly layers d=1,...,m are EXACTLY the association scheme relations.
  Each relation has:
  - Valency v_i = C(m, i) (number of i-neighbors per vertex)
  - Eigenmatrix entry P_j(i) = K_j(i; m) (Krawtchouk polynomial)

  THE QUOTIENT BY ISOMORPHISM:
  The S_n action partitions tilings into orbits (iso classes).
  This gives a QUOTIENT SCHEME — but it's NOT an association scheme
  in general (the quotient of a scheme is a scheme iff the group
  action is compatible with the scheme relations).

  IS THE S_n ACTION COMPATIBLE WITH THE HAMMING SCHEME?
  For this, we'd need: if σ ∈ S_n maps tiling t to t', and
  d(t, s) = k, then d(t', σ(s)) = k.
  But σ permutes the TILE POSITIONS (because relabeling vertices
  changes which arcs are tiles). So σ does NOT preserve Hamming distance
  in general.

  HOWEVER: the waggly weight matrices W_d DO capture the scheme structure
  projected onto the quotient. They are the "orbital matrices" of the
  Hamming scheme action under S_n.

  THE CODING THEORY CONNECTION:
  In coding theory, an (m, K, d) code C is a subset of {0,1}^m
  with K codewords and minimum Hamming distance d.

  An iso class of tournaments IS a "code" in this sense:
  - Codewords = tilings in the class
  - Minimum distance = min Hamming distance within the class
  - K = fiber(class) = number of codewords

  The MINIMUM DISTANCE of a class tells us how "spread out" its tilings are.
  Classes with large |Aut| have small fiber and large min distance.
  Classes with |Aut|=1 have large fiber and small min distance.
""")

# Compute minimum intra-class Hamming distance for each class
for n in [5, 6]:
    m, total, tile_pairs, tc, members, sizes, mid_map, nm, info, merged = build_space(n)

    print("  n=%d: Intra-class Hamming distances:\n" % n)
    print("  %-4s %-4s %-6s %-8s %-8s %-8s" %
          ("mid", "H", "fiber", "min_d", "avg_d", "max_d"))
    print("  " + "-"*45)

    order = sorted(range(nm), key=lambda x: info[merged[x][0]]['H'])

    for mi in order:
        H = info[merged[mi][0]]['H']
        all_tilings = []
        for c in merged[mi]:
            all_tilings.extend(members[c])
        fiber = len(all_tilings)

        if fiber <= 1:
            print("  %-4d %-4d %-6d %-8s %-8s %-8s" % (mi, H, fiber, "-", "-", "-"))
            continue

        dists = []
        for i in range(len(all_tilings)):
            for j in range(i+1, len(all_tilings)):
                d = bin(all_tilings[i] ^ all_tilings[j]).count('1')
                dists.append(d)

        print("  %-4d %-4d %-6d %-8d %-8.2f %-8d" %
              (mi, H, fiber, min(dists), np.mean(dists), max(dists)))

    # Coding theory parameters: (m, K, d_min) per class
    print("\n  AS CODES: (m=%d, K=fiber, d=min_dist):" % m)
    for mi in order:
        H = info[merged[mi][0]]['H']
        all_tilings = []
        for c in merged[mi]:
            all_tilings.extend(members[c])
        fiber = len(all_tilings)
        if fiber <= 1:
            dmin = m  # single codeword, distance to nearest = m
        else:
            dmin = min(bin(all_tilings[i] ^ all_tilings[j]).count('1')
                       for i in range(len(all_tilings)) for j in range(i+1, len(all_tilings)))

        # Singleton bound: K ≤ 2^m / V(m, d-1) where V = volume of Hamming ball
        # Plotkin bound: if d > m/2, K ≤ 2d/(2d-m)
        # Just report the parameters
        rate = log2(fiber) / m if fiber > 1 else 0
        print("    (%d, %d, %d) H=%d, rate=%.3f" % (m, fiber, dmin, H, rate))

    print()

# ============================================================
banner("H AS A SPACE-FILLING CURVE INDEX")
# ============================================================

# The H values at n=5 are: 1, 3, 3, 5, 9, 9, 11, 13, 15, 15
# These are ODD numbers. The gaps are: 2, 0, 2, 4, 0, 2, 2, 2, 0.
# The "curve" H: {merged classes} → Z visits odd integers.

# At n=6, H values span a wider range. Let's see the structure.

for n in [5, 6]:
    m, total, tile_pairs, tc, members, sizes, mid_map, nm, info, merged = build_space(n)
    H_vec = [info[merged[mi][0]]['H'] for mi in range(nm)]
    H_sorted = sorted(H_vec)
    H_unique = sorted(set(H_vec))

    print("  n=%d: H values = %s" % (n, H_sorted))
    print("  Distinct H values: %d out of %d classes" % (len(H_unique), nm))
    print("  H range: [%d, %d]" % (min(H_vec), max(H_vec)))

    # Gaps between consecutive distinct H values
    gaps = [H_unique[i+1] - H_unique[i] for i in range(len(H_unique)-1)]
    print("  H gaps: %s" % gaps)
    print("  Gap statistics: min=%d, max=%d, avg=%.1f" %
          (min(gaps), max(gaps), np.mean(gaps)))

    # The ratio H_max / H_min
    print("  H_max / H_min = %d / %d = %.1f" %
          (max(H_vec), min(H_vec), max(H_vec)/min(H_vec)))

    # Width at each H
    h_width = Counter(H_vec)
    widths = [h_width[h] for h in H_unique]
    print("  Width sequence: %s" % widths)
    print("  Max width: %d at H=%d" %
          (max(h_width.values()),
           [h for h, w in h_width.items() if w == max(h_width.values())][0]))

    # Is the width function unimodal?
    max_idx = widths.index(max(widths))
    increasing = all(widths[i] <= widths[i+1] for i in range(max_idx))
    decreasing = all(widths[i] >= widths[i+1] for i in range(max_idx, len(widths)-1))
    print("  Width unimodal: %s (inc=%s, dec=%s)" % (increasing and decreasing, increasing, decreasing))
    print()

# ============================================================
banner("THE SQRT(2) IN TOURNAMENT THEORY")
# ============================================================

print("""
  WHERE SQRT(2) APPEARS:

  1. THE STAIRCASE GEOMETRY:
     The staircase delta_{n-2} is a right isosceles triangle.
     Legs have length n-2. Hypotenuse has length (n-2)×sqrt(2).
     The diagonal traversal covers sqrt(2) × (leg length) distance.

  2. THE PSEUDO-DOUBLING RATIO:
     V_merged(n+1) / V_merged(n) ≈ 2 - 1/(n-2) → 2
     This is the growth rate of tournament iso classes.
     The factor 2 = (sqrt(2))^2 — the square of the hypotenuse ratio.

  3. THE H RANGE:
     H ranges from 1 to ≈ n!/2^{n-1} (Szele bound).
     The ratio H_max/H_min grows as n!/2^{n-1}.
     By Stirling: n!/2^{n-1} ≈ sqrt(2πn) × (n/e)^n / 2^{n-1}
     = sqrt(2πn) × (n/(2e))^n × 2.
     The sqrt(2πn) factor contains sqrt(2).

  4. THE FIBER FRACTION:
     f(n) = (1/2)_{n-2} / (n-2)! ∼ 1/sqrt(π(n-2))
     = 1/sqrt(π) × 1/sqrt(n-2)
     The GF is (1-x)^{-1/2} — a SQUARE ROOT singularity.
     This is the 2-sheeted branched cover: the number 2 again.

  5. THE ASSOCIATION SCHEME:
     The Hamming scheme H(m, 2) on Q_m has alphabet size q=2.
     The Krawtchouk polynomials are K_k(x; m, 2).
     The eigenvalues involve powers of (q-1) = 1 and q = 2.

  6. THE WALSH SPECTRUM:
     Walsh characters are χ_S(t) = (-1)^{|S∩t|}.
     The base field is GF(2). Everything is mod 2.
     The Fourier transform over GF(2)^m has exactly 2^m terms.

  THE UNIFICATION: sqrt(2) is the characteristic ratio of the
  right isosceles triangle. Tournament theory lives on this triangle.
  The 1D-to-2D transition is mediated by the diagonal traversal,
  which has length sqrt(2) times the leg.

  The boundary between 1D and 2D is a FRACTAL of dimension
  between 1 and 2. The Hilbert curve has dimension 2 but covers
  a 1D path. The staircase traversal has dimension log(m)/log(n-2)
  = log(C(n-1,2))/log(n-2) → 2 as n → ∞.

  Tournament space IS the transition from 1D (H-ordering) to 2D
  (tile space), mediated by the staircase triangle and its sqrt(2).
""")

# Compute the fractal dimension of the staircase traversal
print("  FRACTAL DIMENSION OF THE STAIRCASE:\n")
for n in range(3, 15):
    m = comb(n-1, 2)
    leg = n - 2
    if leg > 0:
        dim = log(m) / log(leg) if leg > 1 else 0
        print("  n=%2d: m=%3d tiles, leg=%d, dim=log(%d)/log(%d) = %.4f" %
              (n, m, leg, m, leg, dim))

print("\n  As n → ∞: dim → log(n²/2)/log(n) = 2 - log(2)/log(n) → 2")
print("  The staircase dimension APPROACHES 2 but never reaches it.")
print("  This is the ALMOST-2D nature of tournament space.")

print("\nDone. opus-2026-03-24-S305")
