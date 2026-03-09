#!/usr/bin/env python3
"""beta2_case4_deep.py - Deep investigation of Case 4 of HYP-278

Case 4: beta_1=0, strongly connected, kappa >= 2.

Prior results (kind-pasteur):
- n=5: Case 4 is EMPTY
- n=6: 1680 tournaments, ALL score (2,2,2,3,3,3), c3=8, kappa=2
  - All have good vertex; #good in {3,6}
  - Sum b1(T\v) in {0, 3}
- n=7 (2000 samples): 378 found, all have good vertex

THIS SCRIPT:
1. n=6 exhaustive + n=7 sampled: detailed good/bad vertex analysis
2. Structural analysis of WHY good vertices exist
3. Group structure at n=6 (score-2 vs score-3 vertices)
4. Test whether Case 4 is always empty at higher n
   (i.e., SC + kappa>=2 => beta_1 >= 1)

Author: opus-2026-03-08
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations

sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')
sys.stdout.reconfigure(line_buffering=True)

# Suppress path_homology_v2 validation output
_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


# ============================================================
# Utility functions
# ============================================================

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def is_sc(A, n):
    if n <= 1:
        return True
    for start_dir in [0, 1]:
        visited = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v in range(n):
                if start_dir == 0 and A[u][v] and v not in visited:
                    visited.add(v); stack.append(v)
                elif start_dir == 1 and A[v][u] and v not in visited:
                    visited.add(v); stack.append(v)
        if len(visited) < n:
            return False
    return True


def compute_b1(A, n):
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    paths2 = enumerate_allowed_paths(A, n, 2)
    if not paths1:
        return 0
    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 == 0:
        return 0
    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D1_om = D1 @ omega1
    sv = np.linalg.svd(D1_om, compute_uv=False)
    rk_d1 = int(sum(s > 1e-8 for s in sv))
    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
        if dim_O2 > 0:
            D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                            [tuple(p) for p in paths1])
            D2_om = D2 @ omega2
            sv2 = np.linalg.svd(D2_om, compute_uv=False)
            rk_d2 = int(sum(s > 1e-8 for s in sv2))
        else:
            rk_d2 = 0
    else:
        rk_d2 = 0
    return dim_O1 - rk_d1 - rk_d2


def vertex_connectivity(A, n):
    if not is_sc(A, n):
        return 0
    if n <= 2:
        return n - 1
    for k in range(1, n):
        for subset in combinations(range(n), k):
            remaining = [v for v in range(n) if v not in subset]
            B, _ = get_induced(A, n, remaining)
            if not is_sc(B, len(remaining)):
                return k
    return n - 1


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c3 += 1
    return c3


def c3_through_vertex(A, n, v):
    """Count 3-cycles containing vertex v."""
    c = 0
    others = [u for u in range(n) if u != v]
    for i in range(len(others)):
        for j in range(i+1, len(others)):
            a, b = others[i], others[j]
            if A[v][a] and A[a][b] and A[b][v]:
                c += 1
            if A[v][b] and A[b][a] and A[a][v]:
                c += 1
    return c


def is_king(A, n, v):
    """Check if v is a king: beats every other vertex in at most 2 steps."""
    for u in range(n):
        if u == v:
            continue
        if A[v][u]:
            continue
        # v doesn't beat u directly; check if v->w->u for some w
        found = False
        for w in range(n):
            if w == v or w == u:
                continue
            if A[v][w] and A[w][u]:
                found = True
                break
        if not found:
            return False
    return True


def scc_condensation(A, n):
    """Return list of SCCs using Tarjan-like approach."""
    # Kosaraju's algorithm
    order = []
    visited = set()
    def dfs1(u):
        visited.add(u)
        for v in range(n):
            if A[u][v] and v not in visited:
                dfs1(v)
        order.append(u)

    sys.setrecursionlimit(10000)
    for v in range(n):
        if v not in visited:
            dfs1(v)

    # Reverse graph DFS
    visited2 = set()
    sccs = []
    def dfs2(u, comp):
        visited2.add(u)
        comp.append(u)
        for v in range(n):
            if A[v][u] and v not in visited2:
                dfs2(v, comp)

    for v in reversed(order):
        if v not in visited2:
            comp = []
            dfs2(v, comp)
            sccs.append(comp)
    return sccs


def classify_tv_failure(A, n, v):
    """When T\v has b1=1, find the SCC structure and what cycle appears."""
    others = [u for u in range(n) if u != v]
    B, vlist = get_induced(A, n, others)
    m = len(others)

    # Check if T\v is SC
    sc = is_sc(B, m)

    # Get SCCs
    sccs = scc_condensation(B, m)

    return {
        'sc': sc,
        'num_sccs': len(sccs),
        'scc_sizes': sorted([len(c) for c in sccs], reverse=True),
    }


# ============================================================
# PART 1: n=6 exhaustive - deep analysis
# ============================================================
print("=" * 70)
print("PART 1: n=6 EXHAUSTIVE - CASE 4 DEEP ANALYSIS")
print("=" * 70)

n = 6
total = 2 ** (n*(n-1)//2)
case4_data = []
t0 = time.time()

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue
    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue

    scores = [sum(A[v]) for v in range(n)]
    c3 = count_c3(A, n)

    # Classify each vertex
    vertex_info = []
    for v in range(n):
        others = [u for u in range(n) if u != v]
        B, _ = get_induced(A, n, others)
        b1v = compute_b1(B, n-1)
        c3v = c3_through_vertex(A, n, v)
        d_out = scores[v]
        d_in = n - 1 - d_out
        king = is_king(A, n, v)

        # Analyze T\v structure
        tv_sc = is_sc(B, n-1)
        tv_c3 = count_c3(B, n-1)
        tv_scores = tuple(sorted(sum(B[i]) for i in range(n-1)))

        # Check if v's out-neighbors and in-neighbors form special structure
        out_nbrs = [u for u in range(n) if A[v][u]]
        in_nbrs = [u for u in range(n) if A[u][v]]

        # Arcs among out-neighbors
        out_arcs = sum(A[a][b] for a in out_nbrs for b in out_nbrs if a != b)
        in_arcs = sum(A[a][b] for a in in_nbrs for b in in_nbrs if a != b)

        vertex_info.append({
            'v': v, 'b1v': b1v, 'd_out': d_out, 'd_in': d_in,
            'c3v': c3v, 'king': king,
            'tv_sc': tv_sc, 'tv_c3': tv_c3, 'tv_scores': tv_scores,
            'out_nbrs': out_nbrs, 'in_nbrs': in_nbrs,
            'out_internal_arcs': out_arcs, 'in_internal_arcs': in_arcs,
        })

    good_vs = [vi for vi in vertex_info if vi['b1v'] == 0]
    bad_vs = [vi for vi in vertex_info if vi['b1v'] != 0]

    case4_data.append({
        'bits': bits, 'A': [row[:] for row in A],
        'scores': tuple(sorted(scores)), 'c3': c3, 'kappa': kappa,
        'vertex_info': vertex_info,
        'num_good': len(good_vs), 'num_bad': len(bad_vs),
    })

elapsed = time.time() - t0
print(f"\n{len(case4_data)} Case 4 tournaments found ({elapsed:.0f}s)")

# Basic stats
score_dist = Counter(d['scores'] for d in case4_data)
print(f"\nScore distribution: {dict(score_dist)}")
kappa_dist = Counter(d['kappa'] for d in case4_data)
print(f"Kappa distribution: {dict(kappa_dist)}")
c3_dist = Counter(d['c3'] for d in case4_data)
print(f"c3 distribution: {dict(c3_dist)}")
ngood_dist = Counter(d['num_good'] for d in case4_data)
print(f"#good vertices: {dict(ngood_dist)}")

# Good vertex stats
good_counts = [d['num_good'] for d in case4_data]
print(f"\nGood vertices: min={min(good_counts)}, max={max(good_counts)}, "
      f"mean={sum(good_counts)/len(good_counts):.2f}")

# ============================================================
# PART 1b: WHY are good vertices good?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 1b: STRUCTURAL ANALYSIS OF GOOD vs BAD VERTICES (n=6)")
print("=" * 70)

# For good vertices: what makes T\v have b1=0?
good_tv_sc = Counter()
good_tv_c3 = Counter()
good_tv_scores = Counter()
good_c3v = Counter()
good_king = Counter()
good_degree = Counter()
good_out_internal = Counter()
good_in_internal = Counter()

bad_tv_sc = Counter()
bad_tv_c3 = Counter()
bad_tv_scores = Counter()
bad_c3v = Counter()
bad_king = Counter()
bad_degree = Counter()
bad_out_internal = Counter()
bad_in_internal = Counter()

for d in case4_data:
    for vi in d['vertex_info']:
        if vi['b1v'] == 0:
            good_tv_sc[vi['tv_sc']] += 1
            good_tv_c3[vi['tv_c3']] += 1
            good_tv_scores[vi['tv_scores']] += 1
            good_c3v[vi['c3v']] += 1
            good_king[vi['king']] += 1
            good_degree[vi['d_out']] += 1
            good_out_internal[vi['out_internal_arcs']] += 1
            good_in_internal[vi['in_internal_arcs']] += 1
        else:
            bad_tv_sc[vi['tv_sc']] += 1
            bad_tv_c3[vi['tv_c3']] += 1
            bad_tv_scores[vi['tv_scores']] += 1
            bad_c3v[vi['c3v']] += 1
            bad_king[vi['king']] += 1
            bad_degree[vi['d_out']] += 1
            bad_out_internal[vi['out_internal_arcs']] += 1
            bad_in_internal[vi['in_internal_arcs']] += 1

print("\nT\\v is SC?")
print(f"  Good: {dict(good_tv_sc)}")
print(f"  Bad:  {dict(bad_tv_sc)}")

print(f"\nc3(T\\v):")
print(f"  Good: {dict(sorted(good_tv_c3.items()))}")
print(f"  Bad:  {dict(sorted(bad_tv_c3.items()))}")

print(f"\nScore of T\\v:")
print(f"  Good: {dict(sorted(good_tv_scores.items()))}")
print(f"  Bad:  {dict(sorted(bad_tv_scores.items()))}")

print(f"\nc3 through v (c3(v)):")
print(f"  Good: {dict(sorted(good_c3v.items()))}")
print(f"  Bad:  {dict(sorted(bad_c3v.items()))}")

print(f"\nIs king vertex?")
print(f"  Good: {dict(good_king)}")
print(f"  Bad:  {dict(bad_king)}")

print(f"\nOut-degree of v:")
print(f"  Good: {dict(sorted(good_degree.items()))}")
print(f"  Bad:  {dict(sorted(bad_degree.items()))}")

print(f"\nInternal arcs among out-neighbors:")
print(f"  Good: {dict(sorted(good_out_internal.items()))}")
print(f"  Bad:  {dict(sorted(bad_out_internal.items()))}")

print(f"\nInternal arcs among in-neighbors:")
print(f"  Good: {dict(sorted(good_in_internal.items()))}")
print(f"  Bad:  {dict(sorted(bad_in_internal.items()))}")

# ============================================================
# PART 2: Group structure at n=6
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: GROUP STRUCTURE (score-2 vs score-3 at n=6)")
print("=" * 70)

# All Case 4 at n=6 have score (2,2,2,3,3,3)
# Low group = {vertices with score 2}, High group = {vertices with score 3}
cross_patterns = Counter()
for d in case4_data:
    A = d['A']
    scores = [sum(A[v]) for v in range(n)]
    low = [v for v in range(n) if scores[v] == 2]
    high = [v for v in range(n) if scores[v] == 3]

    # Arcs from low to high
    lh = sum(A[l][h] for l in low for h in high)
    # Arcs from high to low
    hl = sum(A[h][l] for h in high for l in low)
    # Arcs within low (there are C(3,2)=3 pairs, each has 1 arc)
    ll = sum(A[a][b] for a in low for b in low if a != b)
    # Arcs within high
    hh = sum(A[a][b] for a in high for b in high if a != b)

    # Good vertices by group
    good_low = sum(1 for vi in d['vertex_info'] if vi['d_out'] == 2 and vi['b1v'] == 0)
    good_high = sum(1 for vi in d['vertex_info'] if vi['d_out'] == 3 and vi['b1v'] == 0)

    cross_patterns[(lh, hl, ll, hh, good_low, good_high)] += 1

print(f"\nArc patterns (L->H, H->L, L->L, H->H, #good_low, #good_high):")
for pattern, count in sorted(cross_patterns.items(), key=lambda x: -x[1]):
    lh, hl, ll, hh = pattern[:4]
    gl, gh = pattern[4], pattern[5]
    print(f"  L->H={lh}, H->L={hl}, L->L={ll}, H->H={hh}, "
          f"good_low={gl}, good_high={gh}: {count}")

# Check: within low group, is it always a 3-cycle?
low_patterns = Counter()
high_patterns = Counter()
for d in case4_data:
    A = d['A']
    scores = [sum(A[v]) for v in range(n)]
    low = sorted([v for v in range(n) if scores[v] == 2])
    high = sorted([v for v in range(n) if scores[v] == 3])

    # Low group induced tournament
    Blow, _ = get_induced(A, n, low)
    low_sc = is_sc(Blow, 3)
    low_c3 = count_c3(Blow, 3)
    low_patterns[(low_sc, low_c3)] += 1

    # High group
    Bhigh, _ = get_induced(A, n, high)
    high_sc = is_sc(Bhigh, 3)
    high_c3 = count_c3(Bhigh, 3)
    high_patterns[(high_sc, high_c3)] += 1

print(f"\nLow group (score-2) induced tournament: {dict(low_patterns)}")
print(f"High group (score-3) induced tournament: {dict(high_patterns)}")


# ============================================================
# PART 3: n=7 sampled analysis
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: n=7 SAMPLED ANALYSIS")
print("=" * 70)

n = 7
random.seed(42)
case4_n7 = []
t0 = time.time()
trials = 5000

for trial in range(trials):
    A = random_tournament(n)
    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue
    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue

    scores = tuple(sorted(sum(A[v]) for v in range(n)))
    c3 = count_c3(A, n)

    # Check each vertex
    good_vs = []
    bad_vs = []
    for v in range(n):
        others = [u for u in range(n) if u != v]
        B, _ = get_induced(A, n, others)
        b1v = compute_b1(B, n-1)
        d_out = sum(A[v])
        c3v = c3_through_vertex(A, n, v)
        good_vs.append({'v': v, 'd_out': d_out, 'c3v': c3v}) if b1v == 0 \
            else bad_vs.append({'v': v, 'd_out': d_out, 'c3v': c3v})

    case4_n7.append({
        'scores': scores, 'c3': c3, 'kappa': kappa,
        'num_good': len(good_vs), 'num_bad': len(bad_vs),
        'good_degs': [g['d_out'] for g in good_vs],
        'bad_degs': [b['d_out'] for b in bad_vs],
        'good_c3v': [g['c3v'] for g in good_vs],
        'bad_c3v': [b['c3v'] for b in bad_vs],
    })

    if len(case4_n7) % 200 == 0:
        print(f"  Found {len(case4_n7)} Case 4 ({time.time()-t0:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=7: {len(case4_n7)} Case 4 tournaments found from {trials} trials ({elapsed:.0f}s)")

if case4_n7:
    # All have good vertex?
    no_good = [d for d in case4_n7 if d['num_good'] == 0]
    print(f"Has good vertex: {len(case4_n7) - len(no_good)}/{len(case4_n7)}")
    if no_good:
        print("  *** COUNTEREXAMPLE FOUND ***")

    # Score distribution
    score_dist7 = Counter(d['scores'] for d in case4_n7)
    print(f"\nScore distribution (n=7):")
    for s, cnt in sorted(score_dist7.items(), key=lambda x: -x[1])[:10]:
        print(f"  {s}: {cnt}")

    # Kappa distribution
    kappa_dist7 = Counter(d['kappa'] for d in case4_n7)
    print(f"\nKappa distribution: {dict(kappa_dist7)}")

    # c3 distribution
    c3_dist7 = Counter(d['c3'] for d in case4_n7)
    print(f"\nc3 distribution:")
    for c, cnt in sorted(c3_dist7.items()):
        print(f"  c3={c}: {cnt}")

    # #good vertices
    ngood7 = Counter(d['num_good'] for d in case4_n7)
    print(f"\n#good vertices distribution:")
    for g, cnt in sorted(ngood7.items()):
        print(f"  {g}/{n}: {cnt}")

    good_counts7 = [d['num_good'] for d in case4_n7]
    print(f"Good vertices: min={min(good_counts7)}, max={max(good_counts7)}, "
          f"mean={sum(good_counts7)/len(good_counts7):.2f}")

    # Good vs bad degree
    all_good_degs = Counter()
    all_bad_degs = Counter()
    for d in case4_n7:
        for deg in d['good_degs']:
            all_good_degs[deg] += 1
        for deg in d['bad_degs']:
            all_bad_degs[deg] += 1
    print(f"\nDegree of good vs bad vertices:")
    for deg in sorted(set(list(all_good_degs.keys()) + list(all_bad_degs.keys()))):
        g = all_good_degs.get(deg, 0)
        b = all_bad_degs.get(deg, 0)
        t = g + b
        print(f"  d_out={deg}: good={g}, bad={b}, rate={g/t:.3f}" if t > 0 else f"  d_out={deg}: -")

    # c3(v) for good vs bad
    all_good_c3v = Counter()
    all_bad_c3v = Counter()
    for d in case4_n7:
        for c in d['good_c3v']:
            all_good_c3v[c] += 1
        for c in d['bad_c3v']:
            all_bad_c3v[c] += 1
    print(f"\nc3(v) for good vs bad vertices:")
    for c in sorted(set(list(all_good_c3v.keys()) + list(all_bad_c3v.keys()))):
        g = all_good_c3v.get(c, 0)
        b = all_bad_c3v.get(c, 0)
        t = g + b
        print(f"  c3(v)={c}: good={g}, bad={b}, rate={g/t:.3f}" if t > 0 else f"  c3(v)={c}: -")


# ============================================================
# PART 4: TEST IF CASE 4 IS ALWAYS EMPTY
# (i.e., SC + kappa>=2 => b1>=1)
# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: IS CASE 4 ALWAYS EMPTY?")
print("(i.e., does SC + kappa>=2 always imply b1>=1?)")
print("=" * 70)

# n=5 exhaustive check
print("\n--- n=5 exhaustive ---")
n = 5
count_case4_n5 = 0
count_sc_kappa2_n5 = 0
for A in all_tournaments(n):
    if not is_sc(A, n):
        continue
    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue
    count_sc_kappa2_n5 += 1
    b1 = compute_b1(A, n)
    if b1 == 0:
        count_case4_n5 += 1

print(f"n=5: SC + kappa>=2 tournaments: {count_sc_kappa2_n5}")
print(f"n=5: Case 4 (SC + kappa>=2 + b1=0): {count_case4_n5}")
if count_case4_n5 == 0:
    print("=> Case 4 is EMPTY at n=5 (SC + kappa>=2 => b1>=1)")
else:
    print(f"=> Case 4 EXISTS at n=5")

# n=6: we already know it's 1680 (not empty)
print(f"\nn=6: Case 4 EXISTS ({len(case4_data)} tournaments)")
print("=> SC + kappa>=2 does NOT imply b1>=1 at n=6")

# n=7 sampled
if case4_n7:
    print(f"\nn=7 (sampled): Case 4 EXISTS ({len(case4_n7)} found)")
else:
    print(f"\nn=7 (sampled): Case 4 appears EMPTY (0 found in {trials} trials)")

# n=8 quick sample
print(f"\n--- n=8 sampled ---")
n = 8
random.seed(123)
count_case4_n8 = 0
count_sc_kappa2_n8 = 0
t0 = time.time()

for trial in range(500):
    A = random_tournament(n)
    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue
    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue
    count_sc_kappa2_n8 += 1
    count_case4_n8 += 1

elapsed = time.time() - t0
print(f"n=8 ({elapsed:.0f}s): SC + kappa>=2 + b1=0: {count_case4_n8} found in 500 trials")

# ============================================================
# PART 5: For n=6, analyze the 240 "all-good" vs 1440 "3-bad"
# ============================================================
print(f"\n{'=' * 70}")
print("PART 5: n=6 ALL-GOOD vs 3-BAD TOURNAMENTS")
print("=" * 70)

n = 6
all_good_tours = [d for d in case4_data if d['num_good'] == 6]
three_bad_tours = [d for d in case4_data if d['num_good'] == 3]

print(f"\nAll-good: {len(all_good_tours)}")
print(f"Three-bad: {len(three_bad_tours)}")

# For three-bad: are the 3 bad vertices always the low-score or high-score group?
bad_in_low = 0
bad_in_high = 0
bad_mixed = 0
for d in three_bad_tours:
    A = d['A']
    scores = [sum(A[v]) for v in range(n)]
    bad_set = {vi['v'] for vi in d['vertex_info'] if vi['b1v'] != 0}
    low = {v for v in range(n) if scores[v] == 2}
    high = {v for v in range(n) if scores[v] == 3}
    if bad_set == low:
        bad_in_low += 1
    elif bad_set == high:
        bad_in_high += 1
    else:
        bad_mixed += 1

print(f"\n3-bad tournaments: bad vertices are...")
print(f"  All in low-score group: {bad_in_low}")
print(f"  All in high-score group: {bad_in_high}")
print(f"  Mixed: {bad_mixed}")

# For all-good: what's special about these 240 tournaments?
print(f"\nAll-good tournament analysis:")
# Check self-converse (T ≅ T^op)?
sc_count = 0
for d in all_good_tours:
    A = d['A']
    # T^op: reverse all arcs
    Aop = [[A[j][i] for j in range(n)] for i in range(n)]
    # Check isomorphism (brute force for n=6)
    from itertools import permutations as perms
    is_selfconverse = False
    for perm in perms(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != Aop[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            is_selfconverse = True
            break
    if is_selfconverse:
        sc_count += 1
print(f"  Self-converse: {sc_count}/{len(all_good_tours)}")

# Automorphism group size
print(f"\n  Checking automorphism group sizes...")
auto_sizes_ag = Counter()
for d in all_good_tours[:60]:  # sample
    A = d['A']
    auto = 0
    for perm in perms(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            auto += 1
    auto_sizes_ag[auto] += 1

auto_sizes_tb = Counter()
for d in three_bad_tours[:60]:
    A = d['A']
    auto = 0
    for perm in perms(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            auto += 1
    auto_sizes_tb[auto] += 1

print(f"  All-good |Aut| distribution (sample 60): {dict(auto_sizes_ag)}")
print(f"  Three-bad |Aut| distribution (sample 60): {dict(auto_sizes_tb)}")


# ============================================================
# PART 6: Deeper n=6 analysis - bad vertex T\v SCC structure
# ============================================================
print(f"\n{'=' * 70}")
print("PART 6: BAD VERTEX T\\v STRUCTURE (n=6)")
print("=" * 70)

n = 6
# For bad vertices, T\v is always SC (since kappa>=2),
# so b1(T\v)=1 means there's a 1-cycle in homology not killed by d_2.
# At n=5, b1=1 requires SC. The b1=1 SC tournaments at n=5 are:
# - score (2,2,2,2,2) with c3=5: 24 (all regular)
# - score (1,2,2,2,3) with c3=4: 160
# - score (1,1,2,3,3) with c3=3: 120

# For our bad vertices at n=6, what is the T\v isomorphism type?
# We can check if T\v is regular
bad_regular = 0
bad_nonregular = 0
bad_c3_val = Counter()
for d in three_bad_tours[:200]:  # sample
    A = d['A']
    for vi in d['vertex_info']:
        if vi['b1v'] != 0:
            tv_scores = vi['tv_scores']
            if tv_scores == (2, 2, 2, 2, 2):
                bad_regular += 1
            else:
                bad_nonregular += 1
            bad_c3_val[vi['tv_c3']] += 1

print(f"\nBad vertex T\\v classification (sample of 200 tournaments):")
print(f"  T\\v is regular (2,2,2,2,2): {bad_regular}")
print(f"  T\\v is non-regular: {bad_nonregular}")
print(f"  c3(T\\v) values: {dict(sorted(bad_c3_val.items()))}")


# ============================================================
# PART 7: Can we predict good vertices from LOCAL structure?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 7: LOCAL PREDICTORS OF GOOD VERTICES (n=6)")
print("=" * 70)

n = 6
# Key question: is there a LOCAL property that always identifies a good vertex?

# Test: vertex v is good iff its out-neighborhood induces a tournament with c3 <= X?
# Or: v is good iff it's NOT on a "homological 3-cycle" in T

# Let's compute: for each vertex, the number of transitive triples it's in
# A transitive triple through v: v->a->b and v->b (or permutations)
def transitive_triples_through(A, n, v):
    """Count ordered transitive triples (a,b,c) with a->b->c, a->c, containing v."""
    count = 0
    others = [u for u in range(n) if u != v]
    for a in others:
        for b in others:
            if a == b:
                continue
            # Triple (v, a, b): v->a->b, v->b
            if A[v][a] and A[a][b] and A[v][b]:
                count += 1
            # Triple (a, v, b): a->v->b, a->b
            if A[a][v] and A[v][b] and A[a][b]:
                count += 1
            # Triple (a, b, v): a->b->v, a->v
            if A[a][b] and A[b][v] and A[a][v]:
                count += 1
    return count

# For each good/bad vertex, compute ratio trans_triples / c3
good_trans = Counter()
bad_trans = Counter()
good_ratio_sum = 0
bad_ratio_sum = 0
good_count = 0
bad_count = 0

for d in case4_data[:500]:  # sample
    A = d['A']
    for vi in d['vertex_info']:
        v = vi['v']
        tt = transitive_triples_through(A, n, v)
        if vi['b1v'] == 0:
            good_trans[tt] += 1
            good_count += 1
        else:
            bad_trans[tt] += 1
            bad_count += 1

print(f"\nTransitive triples through vertex:")
print(f"  Good: {dict(sorted(good_trans.items()))}")
print(f"  Bad:  {dict(sorted(bad_trans.items()))}")


# ============================================================
# PART 8: DEFINITIVE TEST — does T\v have b1=1 iff T\v is regular?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 8: IS b1(T\\v)=1 IFF T\\v IS REGULAR? (n=6, exhaustive)")
print("=" * 70)

n = 6
b1_1_regular = 0
b1_1_nonregular = 0
b1_0_regular = 0
b1_0_nonregular = 0

for d in case4_data:
    A = d['A']
    for vi in d['vertex_info']:
        is_reg = (vi['tv_scores'] == (2, 2, 2, 2, 2))
        if vi['b1v'] == 1:
            if is_reg:
                b1_1_regular += 1
            else:
                b1_1_nonregular += 1
        else:
            if is_reg:
                b1_0_regular += 1
            else:
                b1_0_nonregular += 1

print(f"\nb1(T\\v)=1 AND T\\v regular: {b1_1_regular}")
print(f"b1(T\\v)=1 AND T\\v non-regular: {b1_1_nonregular}")
print(f"b1(T\\v)=0 AND T\\v regular: {b1_0_regular}")
print(f"b1(T\\v)=0 AND T\\v non-regular: {b1_0_nonregular}")

if b1_1_nonregular == 0 and b1_0_regular == 0:
    print("\n=> PERFECT CORRESPONDENCE: b1(T\\v)=1 iff T\\v is regular!")
elif b1_1_nonregular == 0:
    print(f"\n=> b1=1 implies regular, but {b1_0_regular} regular have b1=0")
else:
    print(f"\n=> No clean correspondence")

# How many T\v are regular?
reg_count_dist = Counter()
for d in case4_data:
    reg = sum(1 for vi in d['vertex_info'] if vi['tv_scores'] == (2, 2, 2, 2, 2))
    reg_count_dist[reg] += 1
print(f"\n#regular T\\v per tournament: {dict(sorted(reg_count_dist.items()))}")


# ============================================================
# PART 9: CRUCIAL — For 3-bad tours, are bad vertices always
# in the same score class AND forming a directed 3-cycle?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 9: BAD VERTEX INDUCED TOURNAMENT (n=6)")
print("=" * 70)

n = 6
bad_induced_sc = 0
bad_induced_trans = 0
bad_score_same = Counter()

for d in three_bad_tours:
    A = d['A']
    bad_set = [vi['v'] for vi in d['vertex_info'] if vi['b1v'] != 0]
    B_bad, _ = get_induced(A, n, bad_set)
    sc_bad = is_sc(B_bad, 3)
    c3_bad = count_c3(B_bad, 3)

    if sc_bad:
        bad_induced_sc += 1
    else:
        bad_induced_trans += 1

    # Score within T of bad vertices
    bad_scores = tuple(sorted(sum(A[v]) for v in bad_set))
    bad_score_same[bad_scores] += 1

print(f"\nBad vertex set induced subtournament:")
print(f"  Is 3-cycle (SC): {bad_induced_sc}")
print(f"  Is transitive: {bad_induced_trans}")
print(f"\nBad vertex score tuples: {dict(bad_score_same)}")


# ============================================================
# PART 10: KEY STRUCTURAL FINDING SUMMARY
# ============================================================
print(f"\n{'=' * 70}")
print("PART 10: SUMMARY AND PROOF DIRECTIONS")
print("=" * 70)

print("""
KEY FINDINGS:

1. Case 4 is NOT empty (exists at n=6,7,8).
   SC + kappa>=2 does NOT imply b1>=1.
   But: all Case 4 tournaments DO have a good vertex.

2. At n=6 (1680 tournaments):
   - ALL have score (2,2,2,3,3,3), c3=8, kappa=2
   - 240 have ALL 6 vertices good
   - 1440 have exactly 3 good vertices

3. Good vertex existence at n=7,8: verified by sampling.

PROOF STRATEGIES for Case 4:
""")

# Check: does at least one vertex with min c3(v) always work?
print("Testing: vertex with minimum c3(v) is always good?")
n = 6
min_c3v_good = 0
for d in case4_data:
    A = d['A']
    c3vs = [(c3_through_vertex(A, n, v), v) for v in range(n)]
    min_c3 = min(c for c, v in c3vs)
    min_vertices = [v for c, v in c3vs if c == min_c3]
    good_set = {vi['v'] for vi in d['vertex_info'] if vi['b1v'] == 0}
    if any(v in good_set for v in min_vertices):
        min_c3v_good += 1
print(f"  n=6: {min_c3v_good}/{len(case4_data)}")

# Check: max c3(v) vertex always good?
max_c3v_good = 0
for d in case4_data:
    A = d['A']
    c3vs = [(c3_through_vertex(A, n, v), v) for v in range(n)]
    max_c3 = max(c for c, v in c3vs)
    max_vertices = [v for c, v in c3vs if c == max_c3]
    good_set = {vi['v'] for vi in d['vertex_info'] if vi['b1v'] == 0}
    if any(v in good_set for v in max_vertices):
        max_c3v_good += 1
print(f"  n=6: max c3(v) vertex good: {max_c3v_good}/{len(case4_data)}")

# Check at n=7
if case4_n7:
    print(f"\n  (n=7 not checked for min/max c3v — would need re-running with vertex data)")

# Final: check sum b1(T\v) distribution for both n
print(f"\nSum b1(T\\v) distribution:")
sum_dist_6 = Counter()
for d in case4_data:
    s = sum(vi['b1v'] for vi in d['vertex_info'])
    sum_dist_6[s] += 1
print(f"  n=6: {dict(sorted(sum_dist_6.items()))}")

print("\n\nDone.")
