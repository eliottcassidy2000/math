#!/usr/bin/env python3
"""beta2_case4_attack.py - Attack Case 4 of HYP-278

Case 4: kappa>=2 + SC + b1=0 -> exists good vertex.

NEW APPROACHES:
1. For kappa>=2, ALL T\v are SC. So b1(T\v) in {0,1}.
   b1(T\v)=1 iff T\v is "special" SC tournament.
   We need: NOT ALL T\v have b1=1.

2. Key insight from data: Sum_v b1(T\v) in {0, 3} at n=6.
   At n=7 (sampled): Sum_v b1(T\v) <= 3.

3. NEW: Try to understand b1(T\v)=1 in terms of T's structure.
   b1=1 iff the cycle space Z_1 properly contains image(d_2).
   For SC tournaments: Z_1 has dimension (n-2)(n-3)/2 + 1 (or similar).

4. NEW: Look at the COMPLEMENT graph Omega(T\v) and count independent sets.
   b1 connects to independence polynomial evaluation.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def is_sc(A, n):
    if n <= 1:
        return True
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[u][v] and v not in visited:
                visited.add(v)
                stack.append(v)
    if len(visited) < n:
        return False
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[v][u] and v not in visited:
                visited.add(v)
                stack.append(v)
    return len(visited) == n


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


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c3 += 1
    return c3


# ============================================================
# Part 1: Deep analysis of WHICH T\v have b1=1 at n=6
# ============================================================
print("=" * 70)
print("WHICH T\\v HAVE b1=1? (n=6, kappa>=2 case)")
print("=" * 70)

n = 6
total = 2 ** (n * (n - 1) // 2)
case4_data = []
t0 = time.time()

for bits in range(total):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
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

    # Check kappa >= 2: all single deletions SC
    all_sc = True
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        if not is_sc(B, n - 1):
            all_sc = False
            break
    if not all_sc:
        continue

    # Now check b1 for each deletion
    degs = [sum(A[v]) for v in range(n)]
    b1_vals = []
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        b1v = compute_b1(B, n - 1)
        b1_vals.append(b1v)

    bad_vs = [v for v in range(n) if b1_vals[v] == 1]
    good_vs = [v for v in range(n) if b1_vals[v] == 0]

    # For bad vertices, analyze the T\v structure
    bad_info = []
    for v in bad_vs:
        others = [x for x in range(n) if x != v]
        B, vlist = get_induced(A, n, others)
        m = n - 1
        scores_sub = tuple(sorted(sum(B[i]) for i in range(m)))
        c3_sub = count_c3(B, m)

        # What are the in-neighbors and out-neighbors of v?
        in_nbrs = [x for x in range(n) if A[x][v] and x != v]
        out_nbrs = [x for x in range(n) if A[v][x] and x != v]

        bad_info.append({
            'v': v, 'deg': degs[v],
            'scores_sub': scores_sub, 'c3_sub': c3_sub,
            'in_nbrs': len(in_nbrs), 'out_nbrs': len(out_nbrs)
        })

    case4_data.append({
        'bits': bits, 'degs': degs, 'b1_vals': b1_vals,
        'num_bad': len(bad_vs), 'bad_info': bad_info,
        'bad_vs': bad_vs, 'good_vs': good_vs
    })

elapsed = time.time() - t0
print(f"\n{len(case4_data)} kappa>=2 tournaments ({elapsed:.0f}s)")
print(f"  Sum=0 (all good): {sum(1 for d in case4_data if d['num_bad']==0)}")
print(f"  Sum=3 (3 bad): {sum(1 for d in case4_data if d['num_bad']==3)}")

# Score sequences of bad deletions
bad_sub_scores = Counter()
for d in case4_data:
    for info in d['bad_info']:
        bad_sub_scores[info['scores_sub']] += 1

print(f"\nScore sequences of T\\v when v is BAD:")
for s, cnt in sorted(bad_sub_scores.items(), key=lambda x: -x[1]):
    print(f"  {s}: {cnt}")

# Score sequences of good deletions
good_sub_scores = Counter()
for d in case4_data:
    for v in d['good_vs']:
        others = [x for x in range(n) if x != v]
        bits = d['bits']
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        B, _ = get_induced(A, n, others)
        scores_sub = tuple(sorted(sum(B[i]) for i in range(n-1)))
        good_sub_scores[scores_sub] += 1

print(f"\nScore sequences of T\\v when v is GOOD:")
for s, cnt in sorted(good_sub_scores.items(), key=lambda x: -x[1]):
    print(f"  {s}: {cnt}")


# ============================================================
# Part 2: Analyze the 3-cycle structure at bad vertices
# ============================================================
print(f"\n{'=' * 70}")
print("3-CYCLE STRUCTURE OF BAD VERTICES")
print("=" * 70)

# For each bad vertex v, count c3 in T\v AND how many 3-cycles v participates in
c3_participation = Counter()
c3_sub_dist = Counter()
for d in case4_data:
    bits = d['bits']
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    for v in d['bad_vs']:
        # Count 3-cycles containing v
        c3v = 0
        for i in range(n):
            for j in range(i + 1, n):
                if i == v or j == v:
                    continue
                # Check if {v, i, j} form a 3-cycle
                triple = [v, i, j]
                for perm in [(0,1,2), (0,2,1)]:
                    a, b, c = triple[perm[0]], triple[perm[1]], triple[perm[2]]
                    if A[a][b] and A[b][c] and A[c][a]:
                        c3v += 1
                        break  # count each triple once
        c3_participation[c3v] += 1

        # c3 of T\v
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        c3_sub_dist[count_c3(B, n-1)] += 1

    for v in d['good_vs']:
        c3v = 0
        for i in range(n):
            for j in range(i + 1, n):
                if i == v or j == v:
                    continue
                triple = [v, i, j]
                for perm in [(0,1,2), (0,2,1)]:
                    a, b, c = triple[perm[0]], triple[perm[1]], triple[perm[2]]
                    if A[a][b] and A[b][c] and A[c][a]:
                        c3v += 1
                        break
        c3_participation[('good', c3v)] = c3_participation.get(('good', c3v), 0) + 1

print(f"\n3-cycle participation of BAD vertices:")
for k, cnt in sorted((k2, v2) for k2, v2 in c3_participation.items() if isinstance(k2, int)):
    print(f"  c3(v)={k}: {cnt}")
print(f"\n3-cycle participation of GOOD vertices:")
for k, cnt in sorted((k2, v2) for k2, v2 in c3_participation.items() if isinstance(k2, tuple)):
    print(f"  c3(v)={k[1]}: {cnt}")

print(f"\nc3 of T\\v at BAD vertices:")
for c, cnt in sorted(c3_sub_dist.items()):
    print(f"  c3={c}: {cnt}")


# ============================================================
# Part 3: KEY TEST - Try to prove Case 4 via counting argument
# ============================================================
print(f"\n{'=' * 70}")
print("COUNTING ARGUMENT FOR CASE 4")
print("=" * 70)

# At n=6: c3(T) = 8 always. Each 3-cycle uses 3 vertices.
# b1(T)=0 means "enough" of these cycles are boundary.
# Key: Can we show that c3(T\v)=4 for some v?
# (b1=0 at n=5 seems to require c3 <= 4)

# Let's check: what is c3(T\v) for good vs bad at n=5?
print("\nb1 vs c3 at n=5 (exhaustive):")
n5_b1_c3 = Counter()
n5 = 5
for bits in range(2 ** (n5 * (n5 - 1) // 2)):
    A = [[0] * n5 for _ in range(n5)]
    idx = 0
    for i in range(n5):
        for j in range(i + 1, n5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    b1 = compute_b1(A, n5)
    c3 = count_c3(A, n5)
    n5_b1_c3[(b1, c3)] += 1

for (b, c), cnt in sorted(n5_b1_c3.items()):
    print(f"  b1={b}, c3={c}: {cnt}")


# ============================================================
# Part 4: For kappa>=2 case at n=6, check c3(T\v) distribution
# ============================================================
print(f"\n{'=' * 70}")
print("c3(T\\v) FOR KAPPA>=2 CASE AT n=6")
print("=" * 70)

c3_tv_good = Counter()
c3_tv_bad = Counter()

for d in case4_data:
    bits = d['bits']
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        c3v = count_c3(B, n - 1)
        if d['b1_vals'][v] == 0:
            c3_tv_good[c3v] += 1
        else:
            c3_tv_bad[c3v] += 1

print(f"\nc3(T\\v) at GOOD vertices (kappa>=2, n=6):")
for c, cnt in sorted(c3_tv_good.items()):
    print(f"  c3={c}: {cnt}")

print(f"\nc3(T\\v) at BAD vertices (kappa>=2, n=6):")
for c, cnt in sorted(c3_tv_bad.items()):
    print(f"  c3={c}: {cnt}")

# c3 total is 8. If v participates in k of 8 3-cycles, then c3(T\v) = 8 - k.
# Check this:
print(f"\nVerify: c3(T\\v) = c3(T) - c3_containing_v?")
verified = 0
for d in case4_data[:10]:  # Just check first 10
    bits = d['bits']
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        c3_sub = count_c3(B, n-1)

        # Count 3-cycles containing v
        c3_v = 0
        for i in range(n):
            for j in range(i+1, n):
                if i == v or j == v:
                    continue
                triple = sorted([v, i, j])
                a, b, c = triple
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3_v += 1
        if c3_sub == 8 - c3_v:
            verified += 1
        else:
            print(f"  MISMATCH: c3(T\\v)={c3_sub}, 8-c3_v={8-c3_v}")
print(f"  Verified {verified}/60")


# ============================================================
# Part 5: Try b1(T\v) = 1 iff c3(T\v) >= threshold
# ============================================================
print(f"\n{'=' * 70}")
print("b1(T\\v) vs c3(T\\v) RELATIONSHIP")
print("=" * 70)

# At n=5, if b1=1 iff c3 >= some threshold...
print("\nn=5 exhaustive: b1 vs c3 (combining with SC)")
n5_data = Counter()
for bits in range(2 ** (n5 * (n5 - 1) // 2)):
    A = [[0] * n5 for _ in range(n5)]
    idx = 0
    for i in range(n5):
        for j in range(i + 1, n5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    b1 = compute_b1(A, n5)
    c3 = count_c3(A, n5)
    sc = is_sc(A, n5)
    n5_data[(sc, b1, c3)] += 1

for (sc, b, c), cnt in sorted(n5_data.items()):
    print(f"  SC={sc}, b1={b}, c3={c}: {cnt}")


# ============================================================
# Part 6: KING VERTEX APPROACH
# ============================================================
print(f"\n{'=' * 70}")
print("KING VERTEX APPROACH (kappa>=2 case)")
print("=" * 70)

# A king vertex is one that reaches every other in at most 2 steps.
# In tournaments, every vertex of max outdegree is a king.
# For near-regular (score (2,2,2,3,3,3) at n=6), the degree-3 vertices are kings.

# Question: Is a king vertex always good?
king_good = 0
king_bad = 0
non_king_good = 0
non_king_bad = 0

for d in case4_data:
    bits = d['bits']
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    for v in range(n):
        # Check if v is a king
        reachable = set()
        for u in range(n):
            if u != v and A[v][u]:
                reachable.add(u)
        for u in range(n):
            if u != v and u not in reachable:
                for w in reachable:
                    if A[w][u]:
                        reachable.add(u)
                        break

        is_king = (len(reachable) == n - 1)

        if is_king:
            if d['b1_vals'][v] == 0:
                king_good += 1
            else:
                king_bad += 1
        else:
            if d['b1_vals'][v] == 0:
                non_king_good += 1
            else:
                non_king_bad += 1

print(f"  King vertex good: {king_good}, bad: {king_bad}")
print(f"  Non-king vertex good: {non_king_good}, bad: {non_king_bad}")
if king_good + king_bad > 0:
    print(f"  King good rate: {king_good/(king_good+king_bad):.3f}")
if non_king_good + non_king_bad > 0:
    print(f"  Non-king good rate: {non_king_good/(non_king_good+non_king_bad):.3f}")


# ============================================================
# Part 7: DOMINATING SET approach — does every tournament with
# kappa>=2 and b1=0 have a vertex whose OUT-neighborhood
# contains a vertex with b1(T\v)=0?
# ============================================================
print(f"\n{'=' * 70}")
print("INDEPENDENCE NUMBER OF 'BAD' SET")
print("=" * 70)

# Key question: for the 3-bad-vertex case, do bad vertices form
# an independent set in the tournament? (No edges between them?)
# Or do they dominate each other?

bad_patterns = Counter()
for d in case4_data:
    if d['num_bad'] == 0:
        continue
    bits = d['bits']
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    bvs = d['bad_vs']
    # Tournament on bad vertices
    if len(bvs) == 3:
        # Is the sub-tournament on bad vertices a 3-cycle or transitive?
        a, b, c = bvs
        if A[a][b] and A[b][c] and A[c][a]:
            bad_patterns['3-cycle'] += 1
        elif A[a][c] and A[c][b] and A[b][a]:
            bad_patterns['3-cycle'] += 1
        else:
            bad_patterns['transitive'] += 1

        # Degrees of bad vertices in T
        bad_degs = tuple(sorted(sum(A[v]) for v in bvs))
        bad_patterns[('bad_degs', bad_degs)] += 1

print(f"\nBad vertex sub-tournament patterns (kappa>=2, n=6):")
for k, cnt in sorted(bad_patterns.items(), key=lambda x: str(x)):
    print(f"  {k}: {cnt}")


print("\n\nDone.")
