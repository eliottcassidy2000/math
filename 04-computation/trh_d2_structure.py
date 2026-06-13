#!/usr/bin/env python3
"""
trh_d2_structure.py — Investigate structural property making TRH d²=0 hold.

TRH uses:
- Regular paths: (v_0,...,v_m) where v_i→v_{i+1} AND v_{i-1}→v_{i+1} (skip-one)
- Interior-only boundary: ∂(v_0,...,v_m) = sum_{i=1}^{m-1} (-1)^i (v_0,...,v̂_i,...,v_m)

We investigate:
1. Full characterization of n=5 failures (score, 3-cycles, subtournaments, etc.)
2. Whether d²=0 relates to local regularity / doubly regular subtournaments
3. Explicit d² matrix entries for failing tournaments
4. Attempt to find equivalent combinatorial condition
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

# ============================================================
# Core functions
# ============================================================

def get_regular_paths(A, m):
    """Find all regular paths of length m in tournament A."""
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path:
                continue
            if not A[last][v]:
                continue
            # skip-one condition: prev → v (prev beats v)
            if depth >= 1 and not A[prev][v]:
                continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def interior_boundary_matrix(paths_m, paths_m1):
    """Build the interior boundary matrix ∂: paths_m → paths_{m-1}."""
    if not paths_m or not paths_m1:
        return np.zeros((len(paths_m1) if paths_m1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    idx = {p: i for i, p in enumerate(paths_m1)}
    m = len(paths_m[0]) - 1  # path length
    B = np.zeros((len(paths_m1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):  # interior vertices only
            face = path[:i] + path[i+1:]
            if face in idx:
                B[idx[face], j] += (-1)**i
    return B

def check_d2_detailed(A):
    """Check d²=0 for all degrees. Return dict with details."""
    n = A.shape[0]
    all_paths = {}
    for m in range(n):
        all_paths[m] = get_regular_paths(A, m)

    results = {}
    for m in range(2, n):
        p_m = all_paths[m]
        p_m1 = all_paths[m-1]
        p_m2 = all_paths[m-2]
        if p_m and p_m1 and p_m2:
            B_m = interior_boundary_matrix(p_m, p_m1)
            B_m1 = interior_boundary_matrix(p_m1, p_m2)
            d2 = B_m1 @ B_m
            max_val = np.max(np.abs(d2))
            results[m] = {
                'ok': max_val == 0,
                'd2': d2,
                'max_abs': int(max_val),
                'n_paths_m': len(p_m),
                'n_paths_m1': len(p_m1),
                'n_paths_m2': len(p_m2),
                'nonzero_count': int(np.count_nonzero(d2)),
            }
        else:
            results[m] = {
                'ok': True,
                'd2': None,
                'max_abs': 0,
                'n_paths_m': len(p_m),
                'n_paths_m1': len(p_m1) if m-1 in all_paths else 0,
                'n_paths_m2': len(p_m2) if m-2 in all_paths else 0,
                'nonzero_count': 0,
            }
    return results, all_paths

def score_sequence(A):
    n = A.shape[0]
    scores = sorted([int(sum(A[v])) for v in range(n)])
    return tuple(scores)

def count_3cycles(A):
    n = A.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if i,j,k form a 3-cycle
                edges = A[i][j] + A[j][k] + A[k][i]
                rev_edges = A[j][i] + A[k][j] + A[i][k]
                if edges == 3 or rev_edges == 3:
                    count += 1
    return count

def is_vertex_transitive(A):
    """Check if tournament is vertex-transitive (via automorphism group)."""
    n = A.shape[0]
    # Check if all vertices have the same score first
    scores = [int(sum(A[v])) for v in range(n)]
    if len(set(scores)) > 1:
        return False
    # For small n, try all permutations
    if n > 8:
        return None  # too expensive
    for target in range(1, n):
        found = False
        for perm in permutations(range(n)):
            if perm[0] != target:
                continue
            # Check if perm is an automorphism
            is_auto = True
            for i in range(n):
                for j in range(n):
                    if A[i][j] != A[perm[i]][perm[j]]:
                        is_auto = False
                        break
                if not is_auto:
                    break
            if is_auto:
                found = True
                break
        if not found:
            return False
    return True

def get_subtournament(A, vertices):
    """Extract subtournament on given vertices."""
    k = len(vertices)
    B = np.zeros((k, k), dtype=int)
    for i in range(k):
        for j in range(k):
            B[i][j] = A[vertices[i]][vertices[j]]
    return B

def subtournament_iso_type(A):
    """Return canonical form of tournament (minimum over relabelings)."""
    n = A.shape[0]
    min_repr = None
    for perm in permutations(range(n)):
        repr_bits = []
        for i in range(n):
            for j in range(i+1, n):
                repr_bits.append(int(A[perm[i]][perm[j]]))
        repr_tuple = tuple(repr_bits)
        if min_repr is None or repr_tuple < min_repr:
            min_repr = repr_tuple
    return min_repr

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(pairs)):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def tournament_to_bits(A):
    """Convert tournament to bit string for identification."""
    n = A.shape[0]
    bits = []
    for i in range(n):
        for j in range(i+1, n):
            bits.append(int(A[i][j]))
    return tuple(bits)

def has_property_W(A):
    """Check if tournament has 'property W': for every triple (a,b,c),
    if a→b→c and a→c, then there exists d beaten by all of a,b,c
    or d that beats all of a,b,c."""
    n = A.shape[0]
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            for c in range(n):
                if c == a or c == b:
                    continue
                # Check a→b→c and a→c (transitive triple)
                if A[a][b] and A[b][c] and A[a][c]:
                    # Check if there's a vertex beaten by all three
                    found = False
                    for d in range(n):
                        if d == a or d == b or d == c:
                            continue
                        if A[a][d] and A[b][d] and A[c][d]:
                            found = True
                            break
                        if A[d][a] and A[d][b] and A[d][c]:
                            found = True
                            break
                    if not found:
                        return False
    return True

def common_out_neighborhood(A, u, v):
    """Vertices beaten by both u and v."""
    n = A.shape[0]
    return {w for w in range(n) if w != u and w != v and A[u][w] and A[v][w]}

def common_in_neighborhood(A, u, v):
    """Vertices that beat both u and v."""
    n = A.shape[0]
    return {w for w in range(n) if w != u and w != v and A[w][u] and A[w][v]}

# ============================================================
# ANALYSIS 1: Full n=5 characterization
# ============================================================

print("=" * 70)
print("ANALYSIS 1: FULL n=5 CHARACTERIZATION")
print("=" * 70)

n = 5
pass_list = []
fail_list = []
fail_details = []

for A in all_tournaments(n):
    results, all_paths = check_d2_detailed(A)
    d2_ok = all(results[m]['ok'] for m in results)
    ss = score_sequence(A)
    t3 = count_3cycles(A)
    bits = tournament_to_bits(A)

    if d2_ok:
        pass_list.append((bits, ss, t3))
    else:
        # Find first failing degree
        first_fail_m = None
        for m in sorted(results.keys()):
            if not results[m]['ok']:
                first_fail_m = m
                break
        fail_list.append((bits, ss, t3))
        fail_details.append({
            'bits': bits,
            'ss': ss,
            't3': t3,
            'first_fail_m': first_fail_m,
            'results': results,
            'A': A.copy(),
        })

print(f"\nn=5: {len(pass_list)} pass, {len(fail_list)} fail out of {len(pass_list)+len(fail_list)}")

# Score sequence distribution
print("\n--- Score sequence distribution ---")
pass_scores = Counter(x[1] for x in pass_list)
fail_scores = Counter(x[1] for x in fail_list)
all_scores = sorted(set(list(pass_scores.keys()) + list(fail_scores.keys())))
for ss in all_scores:
    p = pass_scores.get(ss, 0)
    f = fail_scores.get(ss, 0)
    pct = f / (p + f) * 100 if (p + f) > 0 else 0
    print(f"  {ss}: {p} pass, {f} fail ({pct:.1f}% fail)")

# 3-cycle distribution
print("\n--- 3-cycle count distribution ---")
pass_t3 = Counter(x[2] for x in pass_list)
fail_t3 = Counter(x[2] for x in fail_list)
all_t3 = sorted(set(list(pass_t3.keys()) + list(fail_t3.keys())))
for t3 in all_t3:
    p = pass_t3.get(t3, 0)
    f = fail_t3.get(t3, 0)
    pct = f / (p + f) * 100 if (p + f) > 0 else 0
    print(f"  t3={t3}: {p} pass, {f} fail ({pct:.1f}% fail)")

# First failing degree
print("\n--- First failing degree m ---")
fail_m_counts = Counter(d['first_fail_m'] for d in fail_details)
for m in sorted(fail_m_counts.keys()):
    print(f"  m={m}: {fail_m_counts[m]} tournaments")

# Vertex transitivity (only for regular tournaments)
print("\n--- Vertex transitivity check (regular tournaments only) ---")
regular_pass = [d for d in pass_list if d[1] == (2, 2, 2, 2, 2)]
regular_fail = [d for d in fail_details if d['ss'] == (2, 2, 2, 2, 2)]
print(f"  Regular (score (2,2,2,2,2)): {len(regular_pass)} pass, {len(regular_fail)} fail")
for d in regular_fail[:3]:
    vt = is_vertex_transitive(d['A'])
    print(f"    Failing regular tournament: VT={vt}, t3={d['t3']}")
for bits, ss, t3 in regular_pass[:3]:
    # Reconstruct A
    A_tmp = np.zeros((n, n), dtype=int)
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for k, (i, j) in enumerate(pairs):
        if bits[k]:
            A_tmp[i][j] = 1
        else:
            A_tmp[j][i] = 1
    vt = is_vertex_transitive(A_tmp)
    print(f"    Passing regular tournament: VT={vt}, t3={t3}")

# ============================================================
# ANALYSIS 2: Subtournament patterns in failing tournaments
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 2: SUBTOURNAMENT PATTERNS")
print("=" * 70)

# Classify all 4-vertex subtournaments
print("\n--- 4-vertex subtournament types in pass vs fail ---")
# Canonical forms of T4 tournaments
t4_types = {}
for A4 in all_tournaments(4):
    canon = subtournament_iso_type(A4)
    if canon not in t4_types:
        ss4 = score_sequence(A4)
        t3_4 = count_3cycles(A4)
        t4_types[canon] = {'ss': ss4, 't3': t3_4, 'count': 0}

# For each n=5 tournament, count subtournament types
def count_sub4_types(A):
    n = A.shape[0]
    type_counts = Counter()
    for verts in combinations(range(n), 4):
        B = get_subtournament(A, verts)
        canon = subtournament_iso_type(B)
        type_counts[canon] += 1
    return type_counts

# Sample: compare pass vs fail
print("\nSubtournament-4 type counts (averaged):")
pass_type_totals = Counter()
fail_type_totals = Counter()
n_pass_sample = min(50, len(pass_list))
n_fail_sample = min(50, len(fail_details))

pairs5 = [(i, j) for i in range(5) for j in range(i+1, 5)]
for idx in range(n_pass_sample):
    bits = pass_list[idx][0]
    A_tmp = np.zeros((5, 5), dtype=int)
    for k, (i, j) in enumerate(pairs5):
        if bits[k]:
            A_tmp[i][j] = 1
        else:
            A_tmp[j][i] = 1
    tc = count_sub4_types(A_tmp)
    for t, c in tc.items():
        pass_type_totals[t] += c

for idx in range(n_fail_sample):
    A_tmp = fail_details[idx]['A']
    tc = count_sub4_types(A_tmp)
    for t, c in tc.items():
        fail_type_totals[t] += c

for canon in sorted(t4_types.keys()):
    info = t4_types[canon]
    p_avg = pass_type_totals.get(canon, 0) / max(n_pass_sample, 1)
    f_avg = fail_type_totals.get(canon, 0) / max(n_fail_sample, 1)
    print(f"  T4 type ss={info['ss']} t3={info['t3']}: pass_avg={p_avg:.2f}, fail_avg={f_avg:.2f}")

# ============================================================
# ANALYSIS 3: Explicit d² matrix for a failing tournament
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 3: EXPLICIT d² MATRIX (first failing tournament)")
print("=" * 70)

d = fail_details[0]
A = d['A']
print(f"\nTournament adjacency matrix:")
print(A)
print(f"Score sequence: {d['ss']}")
print(f"3-cycles: {d['t3']}")
print(f"First failing degree: m={d['first_fail_m']}")

m_fail = d['first_fail_m']
results, all_paths = check_d2_detailed(A)

paths_m = all_paths[m_fail]
paths_m1 = all_paths[m_fail - 1]
paths_m2 = all_paths[m_fail - 2]

print(f"\nPaths at degree {m_fail}: {len(paths_m)}")
print(f"Paths at degree {m_fail-1}: {len(paths_m1)}")
print(f"Paths at degree {m_fail-2}: {len(paths_m2)}")

B_m = interior_boundary_matrix(paths_m, paths_m1)
B_m1 = interior_boundary_matrix(paths_m1, paths_m2)
d2 = B_m1 @ B_m

print(f"\nd² = ∂_{m_fail-1} ∘ ∂_{m_fail} matrix ({d2.shape[0]}×{d2.shape[1]}):")
print(f"Number of nonzero entries: {np.count_nonzero(d2)}")
print(f"Max absolute value: {np.max(np.abs(d2))}")
print(f"Nonzero entries (row_path, col_path, value):")

nonzero_rows, nonzero_cols = np.nonzero(d2)
for r, c in zip(nonzero_rows, nonzero_cols):
    print(f"  d²[{paths_m2[r]}, {paths_m[c]}] = {d2[r,c]}")
    # Show the chain: path_m → faces → faces of faces
    path = paths_m[c]
    print(f"    Source path: {path}")
    print(f"    ∂({path}) contributions:")
    m_len = len(path) - 1
    for i in range(1, m_len):
        face = path[:i] + path[i+1:]
        if face in {p for p in paths_m1}:
            sign = (-1)**i
            print(f"      ({sign:+d}) * {face}")
            # Now boundary of this face
            face_len = len(face) - 1
            for j in range(1, face_len):
                fface = face[:j] + face[j+1:]
                if fface in {p for p in paths_m2}:
                    sign2 = (-1)**j
                    print(f"          → ({sign*sign2:+d}) * {fface}")

# ============================================================
# ANALYSIS 4: What makes faces of faces fail to cancel?
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 4: WHY FACES FAIL TO CANCEL")
print("=" * 70)

# For standard simplicial d²=0:
# deleting i then j (j<i) gives (-1)^i * (-1)^j
# deleting j then i-1 (since i shifted) gives (-1)^j * (-1)^{i-1}
# These cancel: (-1)^{i+j} + (-1)^{i+j-1} = 0
#
# For TRH, the issue is that deleting vertex i might make the resulting
# path NOT regular, so the face is missing from the chain group.

print("\nFor the first failing path, trace which faces are regular vs not:")
for c_idx in range(min(3, len(nonzero_cols))):
    c = nonzero_cols[c_idx]
    path = paths_m[c]
    print(f"\n  Path: {path}")
    m_len = len(path) - 1

    # Check all interior deletions
    for i in range(1, m_len):
        face = path[:i] + path[i+1:]
        is_regular = face in set(paths_m1)
        print(f"    Delete pos {i} (vertex {path[i]}): {face} — regular? {is_regular}")

        if is_regular:
            # Check which faces of this face are regular
            face_len = len(face) - 1
            for j in range(1, face_len):
                fface = face[:j] + face[j+1:]
                is_reg2 = fface in set(paths_m2)
                print(f"      Delete pos {j} (vertex {face[j]}): {fface} — regular? {is_reg2}")

# ============================================================
# ANALYSIS 5: Common out-neighborhood condition
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 5: NEIGHBORHOOD CONDITIONS")
print("=" * 70)

print("\nFor each edge a→b, compute |N+(a) ∩ N+(b)| (common out-neighbors):")

# Check if d²=0 correlates with common neighborhood regularity
def edge_common_stats(A):
    n = A.shape[0]
    stats = []
    for a in range(n):
        for b in range(n):
            if a != b and A[a][b]:
                co = len(common_out_neighborhood(A, a, b))
                ci = len(common_in_neighborhood(A, a, b))
                stats.append((a, b, co, ci))
    return stats

# Compare a passing and failing tournament
print("\nFailing tournament:")
stats_fail = edge_common_stats(fail_details[0]['A'])
co_vals = [s[2] for s in stats_fail]
ci_vals = [s[3] for s in stats_fail]
print(f"  Common out-neighbor counts: {Counter(co_vals)}")
print(f"  Common in-neighbor counts: {Counter(ci_vals)}")

# Find a passing tournament with same score sequence
target_ss = fail_details[0]['ss']
for bits, ss, t3 in pass_list:
    if ss == target_ss:
        A_tmp = np.zeros((5, 5), dtype=int)
        for k, (i, j) in enumerate(pairs5):
            if bits[k]:
                A_tmp[i][j] = 1
            else:
                A_tmp[j][i] = 1
        print(f"\nPassing tournament with same score seq {target_ss}:")
        stats_pass = edge_common_stats(A_tmp)
        co_vals = [s[2] for s in stats_pass]
        ci_vals = [s[3] for s in stats_pass]
        print(f"  Common out-neighbor counts: {Counter(co_vals)}")
        print(f"  Common in-neighbor counts: {Counter(ci_vals)}")
        break

# ============================================================
# ANALYSIS 6: Test specific combinatorial conditions
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 6: COMBINATORIAL CONDITIONS")
print("=" * 70)

# Condition 1: "skip-closure" — if a→b→c with a→c, and b→d→e with b→e,
# and c→d, is a→d forced?
# I.e., does the regularity condition propagate transitively?

# Condition 2: Check if d²=0 ⟺ every regular 3-path has ALL its faces regular
def all_faces_regular(A, m):
    """Check if every regular m-path has all interior faces being regular (m-1)-paths."""
    paths_m_list = get_regular_paths(A, m)
    paths_m1_set = set(get_regular_paths(A, m - 1))

    bad_paths = []
    for path in paths_m_list:
        path_len = len(path) - 1  # = m
        for i in range(1, path_len):
            face = path[:i] + path[i+1:]
            if face not in paths_m1_set:
                bad_paths.append((path, i, face))
    return bad_paths

print("\n--- Condition: All interior faces of regular paths are regular ---")
n = 5
face_pass = 0
face_fail = 0
face_equiv_d2 = 0  # count where this condition ⟺ d²=0

d2_results_all = {}

for idx, A in enumerate(all_tournaments(n)):
    bits = tournament_to_bits(A)
    results, all_paths = check_d2_detailed(A)
    d2_ok = all(results[m]['ok'] for m in results)

    # Check if all faces are regular
    all_faces_ok = True
    for m in range(2, n):
        bad = all_faces_regular(A, m)
        if bad:
            all_faces_ok = False
            break

    if all_faces_ok:
        face_pass += 1
    else:
        face_fail += 1

    if all_faces_ok == d2_ok:
        face_equiv_d2 += 1

    d2_results_all[bits] = d2_ok

total = face_pass + face_fail
print(f"  All-faces-regular: {face_pass} pass, {face_fail} fail")
print(f"  d²=0 ⟺ all-faces-regular? Agreement: {face_equiv_d2}/{total} ({face_equiv_d2/total*100:.1f}%)")

# Condition 3: "Regularity closure" - weaker: check only at the degree where d² fails
print("\n--- Condition: Regularity at degree m where d² is tested ---")
# For d²=0 at degree m, we need ∂_{m-1} ∘ ∂_m = 0
# This means: for each m-path p and (m-2)-path q, the sum
#   sum_{i interior of p} (-1)^i * [face_i(p) is regular] *
#     sum_{j interior of face_i(p)} (-1)^j * [face_j(face_i(p)) = q]
# must equal zero.
# If all faces of faces exist, this is standard simplicial cancellation.
# The issue is when some face_i(p) is NOT regular, breaking the pairing.

# Let's check: is "all faces regular at all degrees" EQUIVALENT to d²=0?
# Already computed above. Let's also check if "all faces regular" is STRONGER.
only_d2 = sum(1 for bits in d2_results_all if d2_results_all[bits] and
              bits not in {tournament_to_bits(A2) for A2 in [] })

# Actually let's directly compare
print(f"\n  all-faces-regular → d²=0? ", end="")
af_implies_d2 = True
d2_implies_af = True
for A in all_tournaments(n):
    bits = tournament_to_bits(A)
    d2_ok = d2_results_all[bits]

    af_ok = True
    for m in range(2, n):
        if all_faces_regular(A, m):
            af_ok = False
            break
    af_ok = not af_ok  # fix: all_faces_regular returns BAD paths

    # Recompute properly
    af_ok2 = True
    for m in range(2, n):
        bad = all_faces_regular(A, m)
        if bad:
            af_ok2 = False
            break

    if af_ok2 and not d2_ok:
        af_implies_d2 = False
    if d2_ok and not af_ok2:
        d2_implies_af = False

print(f"{'YES' if af_implies_d2 else 'NO'}")
print(f"  d²=0 → all-faces-regular? {'YES' if d2_implies_af else 'NO'}")
if af_implies_d2 and d2_implies_af:
    print("  *** EQUIVALENT! d²=0 ⟺ all interior faces are regular ***")

# ============================================================
# ANALYSIS 7: Deeper look — when does a face become non-regular?
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 7: REGULARITY FAILURES IN FACES")
print("=" * 70)

print("\nFor a failing tournament, show WHY faces are non-regular:")
A = fail_details[0]['A']
print(f"Tournament:\n{A}")

for m in range(2, 5):
    bad = all_faces_regular(A, m)
    if bad:
        print(f"\n  At m={m}, {len(bad)} faces are non-regular:")
        for path, i, face in bad[:5]:
            print(f"    Path {path}, delete pos {i} → {face}")
            # Check which regularity condition fails
            for k in range(len(face) - 1):
                v_k = face[k]
                v_k1 = face[k + 1]
                if not A[v_k][v_k1]:
                    print(f"      Edge {v_k}→{v_k1} missing (reversed: {v_k1}→{v_k})")
                if k >= 1:
                    v_km1 = face[k - 1]
                    if not A[v_km1][v_k1]:
                        print(f"      Skip-one {v_km1}→{v_k1} FAILS at pos {k}")

# ============================================================
# ANALYSIS 8: Is d²=0 related to the 3-cycle graph?
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 8: 3-CYCLE GRAPH STRUCTURE")
print("=" * 70)

def three_cycle_graph(A):
    """Build the graph where 3-cycles are vertices, edges connect sharing 2 vertices."""
    n = A.shape[0]
    cycles = []
    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][j] and A[j][k] and A[k][i]:
                    cycle = tuple(sorted([i, j, k]))
                    if cycle not in cycles:
                        cycles.append(cycle)

    # Build adjacency
    edges = []
    for a in range(len(cycles)):
        for b in range(a+1, len(cycles)):
            shared = len(set(cycles[a]) & set(cycles[b]))
            if shared == 2:
                edges.append((a, b))
    return cycles, edges

print("\nFailing tournament 3-cycle graph:")
cycles_f, edges_f = three_cycle_graph(fail_details[0]['A'])
print(f"  {len(cycles_f)} cycles, {len(edges_f)} edges")

# Find a passing tournament with same t3
target_t3 = fail_details[0]['t3']
for bits, ss, t3 in pass_list:
    if t3 == target_t3:
        A_tmp = np.zeros((5, 5), dtype=int)
        for k, (i, j) in enumerate(pairs5):
            if bits[k]:
                A_tmp[i][j] = 1
            else:
                A_tmp[j][i] = 1
        cycles_p, edges_p = three_cycle_graph(A_tmp)
        print(f"Passing tournament with same t3={target_t3}:")
        print(f"  {len(cycles_p)} cycles, {len(edges_p)} edges")
        break

# ============================================================
# ANALYSIS 9: Isomorphism classes
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 9: ISOMORPHISM CLASSES AT n=5")
print("=" * 70)

# Group by isomorphism type
iso_pass = Counter()
iso_fail = Counter()

for bits, ss, t3 in pass_list:
    A_tmp = np.zeros((5, 5), dtype=int)
    for k, (i, j) in enumerate(pairs5):
        if bits[k]:
            A_tmp[i][j] = 1
        else:
            A_tmp[j][i] = 1
    canon = subtournament_iso_type(A_tmp)
    iso_pass[canon] += 1

for d in fail_details:
    canon = subtournament_iso_type(d['A'])
    iso_fail[canon] += 1

all_isos = sorted(set(list(iso_pass.keys()) + list(iso_fail.keys())))
print(f"\n{len(all_isos)} isomorphism classes at n=5")
print(f"  Pass-only classes: {sum(1 for c in all_isos if c not in iso_fail)}")
print(f"  Fail-only classes: {sum(1 for c in all_isos if c not in iso_pass)}")
print(f"  Mixed classes: {sum(1 for c in all_isos if c in iso_pass and c in iso_fail)}")

# Show each class
print("\nDetails by iso class:")
for canon in all_isos:
    p = iso_pass.get(canon, 0)
    f = iso_fail.get(canon, 0)
    # Get score sequence and t3 for this class
    A_tmp = np.zeros((5, 5), dtype=int)
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            A_tmp[i][j] = canon[idx]
            A_tmp[j][i] = 1 - canon[idx]
            idx += 1
    ss = score_sequence(A_tmp)
    t3 = count_3cycles(A_tmp)
    status = "PASS" if f == 0 else ("FAIL" if p == 0 else "MIXED")
    print(f"  ss={ss} t3={t3}: {p} pass + {f} fail = {p+f} total [{status}]")

# ============================================================
# ANALYSIS 10: Explicit face regularity check — the KEY condition
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 10: THE ALGEBRAIC OBSTACLE — DETAILED")
print("=" * 70)

# Pick the simplest failing case and trace d² ≠ 0 completely
A = fail_details[0]['A']
m_fail = fail_details[0]['first_fail_m']

paths_m_list = get_regular_paths(A, m_fail)
paths_m1_list = get_regular_paths(A, m_fail - 1)
paths_m1_set = set(paths_m1_list)
paths_m2_list = get_regular_paths(A, m_fail - 2)
paths_m2_set = set(paths_m2_list)

print(f"\nFailing at m={m_fail}")
print(f"Regular {m_fail}-paths: {paths_m_list}")
print(f"Regular {m_fail-1}-paths: {paths_m1_list}")
print(f"Regular {m_fail-2}-paths: {paths_m2_list}")

# For each m-path, compute ∂² contribution manually
print(f"\nManual ∂² computation:")
for path in paths_m_list:
    print(f"\n  ∂({path}):")
    m_len = len(path) - 1
    boundary_terms = []
    for i in range(1, m_len):
        face = path[:i] + path[i+1:]
        if face in paths_m1_set:
            sign = (-1)**i
            boundary_terms.append((sign, face))
            print(f"    ({sign:+d}) * {face}")
        else:
            print(f"    pos {i}: {face} NOT REGULAR — missing from chain group")

    print(f"  ∂²({path}):")
    d2_contribs = defaultdict(int)
    for sign1, face in boundary_terms:
        face_len = len(face) - 1
        for j in range(1, face_len):
            fface = face[:j] + face[j+1:]
            if fface in paths_m2_set:
                sign2 = (-1)**j
                d2_contribs[fface] += sign1 * sign2
                print(f"    ∂({face}): ({sign1*sign2:+d}) * {fface}")

    # Show net contributions
    net = {k: v for k, v in d2_contribs.items() if v != 0}
    if net:
        print(f"  *** NET NONZERO: {dict(net)} ***")
    else:
        print(f"  Net: all cancel to 0")

# ============================================================
# ANALYSIS 11: "Flag" condition — for each pair (i,j) in {1,...,m-1},
# does deleting i then j always give a regular path?
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 11: DOUBLE-DELETION REGULARITY")
print("=" * 70)

print("\nIn standard simplicial homology, d²=0 because deleting pos i then j")
print("pairs with deleting pos j then i-1 (or i then j+1), and signs cancel.")
print("In TRH, this fails when one deletion produces a non-regular path.")
print("\nChecking: for PASSING tournaments, are ALL double-deletions regular?")

n = 5
pass_all_double = 0
pass_some_missing = 0
fail_all_double = 0
fail_some_missing = 0

sample_pass = []
sample_fail = []

for A in all_tournaments(n):
    bits = tournament_to_bits(A)
    d2_ok = d2_results_all[bits]

    has_missing = False
    for m in range(2, n):
        paths_m_list2 = get_regular_paths(A, m)
        paths_m1_set2 = set(get_regular_paths(A, m - 1))
        for path in paths_m_list2:
            path_len = len(path) - 1
            for i in range(1, path_len):
                face = path[:i] + path[i+1:]
                if face not in paths_m1_set2:
                    has_missing = True
                    break
            if has_missing:
                break
        if has_missing:
            break

    if d2_ok:
        if has_missing:
            pass_some_missing += 1
            if len(sample_pass) < 3:
                sample_pass.append(A.copy())
        else:
            pass_all_double += 1
    else:
        if has_missing:
            fail_some_missing += 1
        else:
            fail_all_double += 1
            if len(sample_fail) < 3:
                sample_fail.append(A.copy())

print(f"\n  d²=0 AND all faces regular: {pass_all_double}")
print(f"  d²=0 BUT some face non-regular: {pass_some_missing}")
print(f"  d²≠0 AND some face non-regular: {fail_some_missing}")
print(f"  d²≠0 BUT all faces regular: {fail_all_double}")
print()
if pass_some_missing > 0:
    print(f"  → d²=0 does NOT require all faces regular ({pass_some_missing} counterexamples)")
    print(f"  → d²=0 can hold even when ∂ maps to paths outside the chain group")
    print(f"     (meaning some terms are 'projected out' but the remaining ones still cancel)")
if fail_all_double > 0:
    print(f"  → all-faces-regular does NOT imply d²=0 ({fail_all_double} counterexamples)")
    print(f"     This would be very surprising!")

# ============================================================
# ANALYSIS 12: What if d²=0 ⟺ the tournament avoids a specific
# set of forbidden subtournaments?
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 12: FORBIDDEN SUBTOURNAMENT CHARACTERIZATION")
print("=" * 70)

# At n=5, there are 12 isomorphism classes of tournaments
# Check which iso classes of T4 subtournaments appear ONLY in failing tournaments

# Get T4 iso classes
t4_canon_set = set()
for A4 in all_tournaments(4):
    t4_canon_set.add(subtournament_iso_type(A4))
t4_canon_list = sorted(t4_canon_set)
print(f"\n{len(t4_canon_list)} isomorphism classes of 4-vertex tournaments")

# For each n=5 tournament, check which T4 types it contains
t4_in_pass = set()
t4_in_fail = set()

for bits, ss, t3 in pass_list:
    A_tmp = np.zeros((5, 5), dtype=int)
    for k, (i, j) in enumerate(pairs5):
        if bits[k]:
            A_tmp[i][j] = 1
        else:
            A_tmp[j][i] = 1
    for verts in combinations(range(5), 4):
        B = get_subtournament(A_tmp, verts)
        t4_in_pass.add(subtournament_iso_type(B))

for d in fail_details:
    for verts in combinations(range(5), 4):
        B = get_subtournament(d['A'], verts)
        t4_in_fail.add(subtournament_iso_type(B))

only_in_fail = t4_in_fail - t4_in_pass
only_in_pass = t4_in_pass - t4_in_fail
in_both = t4_in_pass & t4_in_fail

print(f"  T4 types in pass-only: {len(only_in_pass)}")
print(f"  T4 types in fail-only: {len(only_in_fail)}")
print(f"  T4 types in both: {len(in_both)}")

if only_in_fail:
    print("\n  T4 types appearing ONLY in failing tournaments:")
    for canon in sorted(only_in_fail):
        A_tmp = np.zeros((4, 4), dtype=int)
        idx = 0
        for i in range(4):
            for j in range(i+1, 4):
                A_tmp[i][j] = canon[idx]
                A_tmp[j][i] = 1 - canon[idx]
                idx += 1
        ss = score_sequence(A_tmp)
        t3 = count_3cycles(A_tmp)
        print(f"    ss={ss} t3={t3}: {canon}")

# Check: does AVOIDING those T4 types ⟺ d²=0?
if only_in_fail:
    print("\n  Testing: d²=0 ⟺ avoids these T4 types?")
    equiv_count = 0
    for A in all_tournaments(5):
        bits = tournament_to_bits(A)
        d2_ok = d2_results_all[bits]
        has_forbidden = False
        for verts in combinations(range(5), 4):
            B = get_subtournament(A, verts)
            if subtournament_iso_type(B) in only_in_fail:
                has_forbidden = True
                break
        avoids = not has_forbidden
        if avoids == d2_ok:
            equiv_count += 1
    total = 2**10
    print(f"    Agreement: {equiv_count}/{total} ({equiv_count/total*100:.1f}%)")

# ============================================================
# ANALYSIS 13: Check "Property W" and other named properties
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 13: NAMED TOURNAMENT PROPERTIES")
print("=" * 70)

print("\n--- Property W (every transitive triple has a common friend/enemy) ---")
w_pass = 0
w_fail = 0
d2_and_w = 0
d2_and_not_w = 0
not_d2_and_w = 0
not_d2_and_not_w = 0

for A in all_tournaments(5):
    bits = tournament_to_bits(A)
    d2_ok = d2_results_all[bits]
    w = has_property_W(A)
    if d2_ok and w:
        d2_and_w += 1
    elif d2_ok and not w:
        d2_and_not_w += 1
    elif not d2_ok and w:
        not_d2_and_w += 1
    else:
        not_d2_and_not_w += 1

print(f"  d²=0 AND W: {d2_and_w}")
print(f"  d²=0 AND ¬W: {d2_and_not_w}")
print(f"  d²≠0 AND W: {not_d2_and_w}")
print(f"  d²≠0 AND ¬W: {not_d2_and_not_w}")

# ============================================================
# ANALYSIS 14: "Strong path closure" — the regularity condition
# propagates through the skip-one condition
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 14: SKIP-ONE CLOSURE PROPERTY")
print("=" * 70)

def skip_closure_check(A):
    """Check: for every triple a→b→c with a→c (transitive),
    and every d with b→d and c→d, is a→d?
    This would mean: extending a regular path preserves skip-one."""
    n = A.shape[0]
    violations = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b:
                    continue
                if not A[b][c] or not A[a][c]:
                    continue
                # a→b→c with a→c (transitive triple, regular path a,b,c)
                for d in range(n):
                    if d in (a, b, c):
                        continue
                    if A[c][d] and A[b][d]:
                        # b→d and c→d, so (b,c,d) continues with skip-one
                        # For (a,b,c,d) to be regular: need a→c (yes) and b→d (yes)
                        # Actually regularity of (a,b,c,d) needs:
                        # a→b, b→c, c→d (edges)
                        # a→c, b→d (skip-one)
                        # So this is automatically regular! No violation here.
                        pass
    return violations

# More interesting: check if the FACE regularity condition has a nice form
# A face of (v0,v1,...,vm) deleting vi gives (v0,...,v_{i-1},v_{i+1},...,vm)
# This is regular iff:
#   v_{j}→v_{j+1} for all j (edges) — these are inherited from original except
#     the new edge v_{i-1}→v_{i+1} which may not be in the tournament
#   v_{j-1}→v_{j+1} for all j≥1 (skip-one) — new skip-ones:
#     v_{i-2}→v_{i+1} (if i≥2) and v_{i-1}→v_{i+2} (if i≤m-2)

print("\nFace regularity requires THREE conditions when deleting interior vertex v_i:")
print("  1. v_{i-1} → v_{i+1} (new edge in face)")
print("  2. v_{i-2} → v_{i+1} (new skip-one, if i≥2)")
print("  3. v_{i-1} → v_{i+2} (new skip-one, if i≤m-2)")
print("\nCondition 1 (edge): v_{i-1} beats v_{i+1}")
print("Condition 2 (left skip): v_{i-2} beats v_{i+1}")
print("Condition 3 (right skip): v_{i-1} beats v_{i+2}")

print("\nFor d²=0, we need: whenever conditions 1,2,3 create a 'missing' face,")
print("the missing contributions must cancel in the alternating sum.")
print("This is an algebraic cancellation, NOT just a combinatorial property.")

# ============================================================
# ANALYSIS 15: Direct combinatorial characterization attempt
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 15: TESTING CANDIDATE EQUIVALENT CONDITIONS")
print("=" * 70)

# Condition A: "2-regular" — every vertex has the same out-degree
# Condition B: doubly regular — every pair has same # common out-neighbors
# Condition C: every 3-element subset spans a 3-cycle
# Condition D: the complement tournament also has d²=0

def is_doubly_regular(A):
    """Every pair of vertices has same # common out-neighbors."""
    n = A.shape[0]
    vals = set()
    for i in range(n):
        for j in range(i+1, n):
            co = len(common_out_neighborhood(A, i, j))
            vals.add(co)
    return len(vals) <= 2  # at most 2 values (for i→j vs j→i)

def complement_tournament(A):
    """Reverse all edges."""
    return 1 - A - np.eye(A.shape[0], dtype=int)

print("\n--- Condition: doubly regular ---")
dr_and_d2 = 0
dr_and_not_d2 = 0
not_dr_and_d2 = 0
not_dr_and_not_d2 = 0

for A in all_tournaments(5):
    bits = tournament_to_bits(A)
    d2_ok = d2_results_all[bits]
    dr = is_doubly_regular(A)
    if dr and d2_ok:
        dr_and_d2 += 1
    elif dr and not d2_ok:
        dr_and_not_d2 += 1
    elif not dr and d2_ok:
        not_dr_and_d2 += 1
    else:
        not_dr_and_not_d2 += 1

print(f"  DR ∧ d²=0: {dr_and_d2}")
print(f"  DR ∧ d²≠0: {dr_and_not_d2}")
print(f"  ¬DR ∧ d²=0: {not_dr_and_d2}")
print(f"  ¬DR ∧ d²≠0: {not_dr_and_not_d2}")

print("\n--- Condition: complement also has d²=0 ---")
comp_same = 0
for A in all_tournaments(5):
    bits = tournament_to_bits(A)
    d2_ok = d2_results_all[bits]
    A_comp = complement_tournament(A)
    bits_comp = tournament_to_bits(A_comp)
    d2_comp = d2_results_all.get(bits_comp, None)
    if d2_comp is not None and d2_ok == d2_comp:
        comp_same += 1

print(f"  d²=0(T) ⟺ d²=0(T^op): {comp_same}/1024 ({comp_same/1024*100:.1f}%)")

# ============================================================
# ANALYSIS 16: Check n=6 statistics for confirmation
# ============================================================

print(f"\n{'='*70}")
print("ANALYSIS 16: n=6 CONFIRMATION (sampling)")
print("=" * 70)

import random
random.seed(42)

n = 6
pairs6 = [(i, j) for i in range(n) for j in range(i+1, n)]
n_sample = 500
n6_pass = 0
n6_fail = 0
n6_fail_m = Counter()
n6_pass_scores = Counter()
n6_fail_scores = Counter()

# Also track face-regularity vs d²=0 correlation
n6_face_d2_agree = 0

for _ in range(n_sample):
    bits = random.randint(0, 2**15 - 1)
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs6):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    results, all_paths = check_d2_detailed(A)
    d2_ok = all(results[m_]['ok'] for m_ in results)
    ss = score_sequence(A)

    # Face regularity
    af_ok = True
    for m_ in range(2, n):
        bad = all_faces_regular(A, m_)
        if bad:
            af_ok = False
            break

    if af_ok == d2_ok:
        n6_face_d2_agree += 1

    if d2_ok:
        n6_pass += 1
        n6_pass_scores[ss] += 1
    else:
        n6_fail += 1
        for m_ in sorted(results.keys()):
            if not results[m_]['ok']:
                n6_fail_m[m_] += 1
                break
        n6_fail_scores[ss] += 1

print(f"\nn=6 sample ({n_sample}): {n6_pass} pass, {n6_fail} fail ({n6_fail/n_sample*100:.1f}% fail)")
print(f"  First failing degree: {dict(n6_fail_m)}")
print(f"  face-regular ⟺ d²=0 agreement: {n6_face_d2_agree}/{n_sample} ({n6_face_d2_agree/n_sample*100:.1f}%)")

# ============================================================
# SUMMARY
# ============================================================

print(f"\n{'='*70}")
print("SUMMARY OF FINDINGS")
print("=" * 70)

print("""
Key results:
1. At n=5: d²=0 depends on the isomorphism class of the tournament.
2. The "all faces regular" condition is related but may not be equivalent.
3. d²=0 is preserved under complementation (T → T^op).
4. The algebraic obstacle: when deleting interior vertex v_i from a regular path,
   the face may fail to be regular because:
   (a) v_{i-1} → v_{i+1} (new edge) might go the wrong way
   (b) v_{i-2} → v_{i+1} or v_{i-1} → v_{i+2} (new skip-ones) might fail
5. For d²=0, missing faces must produce cancellation in the alternating sum.
   This is an algebraic condition, not purely combinatorial.
""")

print("\nDONE.")
