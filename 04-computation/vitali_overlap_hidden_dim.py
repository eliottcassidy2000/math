#!/usr/bin/env python3
"""
vitali_overlap_hidden_dim.py -- kind-pasteur-2026-03-13-S61

Deep exploration: How does the {2,1,0} overlap weight structure between
3-cycles connect to the Vitali non-measurable dimension (c7_dir)?

KEY QUESTION: The labeled lambda graph has entries lambda_{uv} in {0,1,2,...,5}.
These are the overlap weights between VERTEX PAIRS. But the conflict graph
Omega(T) has overlap weights between CYCLES. The lambda graph is a "projection"
of the full cycle overlap structure onto vertex pairs.

What information is lost in this projection?

ANALYSIS PLAN:
1. For the two deepest-ambiguity tournaments (bits 4728, 4658):
   - Full {2,1,0} overlap matrix between ALL cycles (not just 3-cycles)
   - Spectrum of the overlap matrix
   - Which 5-cycles and 7-cycles differ?

2. The "hidden parity" interpretation:
   - c7_dir mod 2 = parity of Hamiltonian cycle count
   - Does this parity have a topological/combinatorial interpretation?

3. The overlap weight as a metric:
   - Do the overlap weights define a metric on cycles?
   - Is there a natural "distance" that separates c7=8 from c7=9?

4. Higher-order overlap invariants:
   - Products and sums of overlap weights around cliques in Omega
   - "Overlap energy" = sum of W[i,j]^2

5. The {2,1,0} structure as a graded ring:
   - Overlap 2 = "strongly conflicting" (2 shared vertices)
   - Overlap 1 = "weakly conflicting" (1 shared vertex)
   - Overlap 0 = "independent" (disjoint)
   - Do these form a filtration that captures the hidden dim?

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
from collections import defaultdict
import math


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 1:
        return 0
    if k == 2:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def get_all_directed_cycles(A, n, max_k=None):
    """Get ALL directed cycles as (vertex_set, length, orientation) tuples.
    orientation = +1 for 'canonical' direction, -1 for reverse."""
    if max_k is None:
        max_k = n
    cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            n_cyc = count_directed_ham_cycles_on_subset(A, verts)
            if n_cyc > 0:
                # For each vertex set with cycles, record it with the count
                cycles.append((frozenset(subset), k, n_cyc))
    return cycles


def get_labeled_lambda(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val
    return lam, c3_sets


def get_3cycle_orientation(A, c3):
    """Return +1 for canonical direction (a->b->c->a where a<b<c), -1 for reverse."""
    a, b, c = sorted(c3)
    if A[a][b] and A[b][c] and A[c][a]:
        return +1
    else:
        return -1


def overlap_weight(c1, c2):
    """Overlap weight = number of shared vertices."""
    return len(c1 & c2)


n = 7

print("=" * 70)
print("VITALI OVERLAP HIDDEN DIMENSION ANALYSIS")
print("=" * 70)

# The two deepest-ambiguity tournaments
targets = [(4728, "H=109, c7=8"), (4658, "H=111, c7=9")]

all_data = {}

for bits, label in targets:
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    lam, c3_sets = get_labeled_lambda(A, n)

    # Get ALL directed cycles with their counts
    cycles = get_all_directed_cycles(A, n)

    c3_count = sum(n_cyc for fs, k, n_cyc in cycles if k == 3)
    c5_count = sum(n_cyc for fs, k, n_cyc in cycles if k == 5)
    c7_count = sum(n_cyc for fs, k, n_cyc in cycles if k == 7)

    c3_vsets = [fs for fs, k, _ in cycles if k == 3]
    c5_vsets = [fs for fs, k, _ in cycles if k == 5]
    c7_vsets = [fs for fs, k, _ in cycles if k == 7]

    print(f"\n{'='*60}")
    print(f"  {label} (bits={bits})")
    print(f"{'='*60}")
    print(f"  H = {H}")
    print(f"  3-cycle vertex sets: {len(c3_vsets)}, directed: {c3_count}")
    print(f"  5-cycle vertex sets: {len(c5_vsets)}, directed: {c5_count}")
    print(f"  7-cycle vertex sets: {len(c7_vsets)}, directed: {c7_count}")

    # 3-cycle orientations
    orientations = {}
    for c3 in c3_sets:
        orientations[c3] = get_3cycle_orientation(A, c3)

    pos_count = sum(1 for v in orientations.values() if v == +1)
    neg_count = sum(1 for v in orientations.values() if v == -1)
    print(f"\n  3-cycle orientations: +1={pos_count}, -1={neg_count}")

    # ======================================================================
    # ANALYSIS 1: Full overlap weight matrix between 3-cycle vertex sets
    # ======================================================================
    print(f"\n  --- 3-CYCLE OVERLAP WEIGHT MATRIX ---")

    n3 = len(c3_sets)
    ov_counts = defaultdict(int)
    for i in range(n3):
        for j in range(i+1, n3):
            ov = overlap_weight(c3_sets[i], c3_sets[j])
            ov_counts[ov] += 1

    for ov in sorted(ov_counts):
        print(f"    overlap={ov}: {ov_counts[ov]} pairs")

    # ======================================================================
    # ANALYSIS 2: Overlap between 3-cycles and 5-cycles
    # ======================================================================
    print(f"\n  --- 3-CYCLE x 5-CYCLE OVERLAP ---")

    c35_ov = defaultdict(int)
    for c3 in c3_sets:
        for c5 in c5_vsets:
            ov = overlap_weight(c3, c5)
            c35_ov[ov] += 1

    for ov in sorted(c35_ov):
        print(f"    overlap={ov}: {c35_ov[ov]} pairs")

    # ======================================================================
    # ANALYSIS 3: Overlap between 5-cycles (this is where c7 hides)
    # ======================================================================
    print(f"\n  --- 5-CYCLE x 5-CYCLE OVERLAP ---")

    n5 = len(c5_vsets)
    c55_ov = defaultdict(int)
    for i in range(n5):
        for j in range(i+1, n5):
            ov = overlap_weight(c5_vsets[i], c5_vsets[j])
            c55_ov[ov] += 1

    for ov in sorted(c55_ov):
        print(f"    overlap={ov}: {c55_ov[ov]} pairs")

    # ======================================================================
    # ANALYSIS 4: Per-5-cycle-vertex-set directed cycle counts
    # ======================================================================
    print(f"\n  --- PER-5-CYCLE-VERTEX-SET DETAILS ---")

    c5_details = []
    for fs, k, n_cyc in cycles:
        if k == 5:
            c5_details.append((sorted(fs), n_cyc))

    for verts, nc in sorted(c5_details):
        print(f"    {verts}: {nc} directed 5-cycles")

    # ======================================================================
    # ANALYSIS 5: Overlap energy = sum of W[i,j]^2
    # ======================================================================
    print(f"\n  --- OVERLAP ENERGY ---")

    # 3-cycle overlap energy
    ov_energy_3 = sum(ov**2 * cnt for ov, cnt in ov_counts.items())
    ov_linear_3 = sum(ov * cnt for ov, cnt in ov_counts.items())
    print(f"    3-cycle overlap energy (sum W^2): {ov_energy_3}")
    print(f"    3-cycle overlap linear (sum W): {ov_linear_3}")

    # ======================================================================
    # ANALYSIS 6: The "signed overlap" (orientation-weighted)
    # ======================================================================
    print(f"\n  --- SIGNED OVERLAP (orientation-weighted) ---")

    # For each pair of 3-cycles with overlap w, compute:
    # sigma(i,j) = w * orient(i) * orient(j)
    signed_ov_sum = {0: 0, 1: 0, 2: 0, 3: 0}
    signed_ov_count = {0: 0, 1: 0, 2: 0, 3: 0}

    for i in range(n3):
        for j in range(i+1, n3):
            ov = overlap_weight(c3_sets[i], c3_sets[j])
            sign_prod = orientations[c3_sets[i]] * orientations[c3_sets[j]]
            signed_ov_sum[ov] += sign_prod
            signed_ov_count[ov] += 1

    for ov in sorted(signed_ov_sum):
        total = signed_ov_count[ov]
        if total > 0:
            print(f"    overlap={ov}: sign_sum={signed_ov_sum[ov]:+d} (of {total} pairs)")

    total_signed = sum(signed_ov_sum.values())
    print(f"    TOTAL signed overlap: {total_signed:+d}")

    # ======================================================================
    # ANALYSIS 7: Vertex participation in cycles of each length
    # ======================================================================
    print(f"\n  --- VERTEX PARTICIPATION ---")

    for length_label, cycle_list in [("3-cycles", c3_sets), ("5-cycles", c5_vsets)]:
        vertex_count = [0] * n
        for c in cycle_list:
            for v in c:
                vertex_count[v] += 1
        print(f"    {length_label}: {vertex_count}")

    # c7 participation (trivially all vertices)
    print(f"    7-cycles: [trivially all {n} vertices in each of {c7_count} directed cycles]")

    # ======================================================================
    # ANALYSIS 8: The "overlap graph" at each weight level
    # ======================================================================
    print(f"\n  --- OVERLAP GRAPHS AT EACH WEIGHT LEVEL ---")

    for w_level in [0, 1, 2]:
        # Build graph on 3-cycle vertex sets where edge iff overlap == w_level
        edge_count = 0
        degree_at_w = [0] * n3
        for i in range(n3):
            for j in range(i+1, n3):
                if overlap_weight(c3_sets[i], c3_sets[j]) == w_level:
                    edge_count += 1
                    degree_at_w[i] += 1
                    degree_at_w[j] += 1
        deg_seq = tuple(sorted(degree_at_w))
        print(f"    W={w_level}: {edge_count} edges, deg_seq={deg_seq}")

    # ======================================================================
    # ANALYSIS 9: The lambda-0 edges (the crucial ones)
    # ======================================================================
    print(f"\n  --- LAMBDA-0 EDGES (non-3-cycle edges) ---")

    lambda0_edges = []
    for u in range(n):
        for v in range(u+1, n):
            if lam[u][v] == 0:
                lambda0_edges.append((u, v))

    print(f"    Lambda-0 edges: {lambda0_edges}")

    # For each lambda-0 edge, how many 5-cycles contain both vertices?
    for u, v in lambda0_edges:
        c5_containing = sum(1 for c5 in c5_vsets if u in c5 and v in c5)
        print(f"      ({u},{v}): in {c5_containing} 5-cycle vertex sets")

    # ======================================================================
    # ANALYSIS 10: The "hidden parity" -- c7 mod 2
    # ======================================================================
    print(f"\n  --- HIDDEN PARITY ---")
    print(f"    c7_dir = {c7_count}")
    print(f"    c7_dir mod 2 = {c7_count % 2}")

    # c3 + c5 + c7 = alpha_1 (from OCF)
    alpha_1 = c3_count + c5_count + c7_count
    print(f"    alpha_1 = c3+c5+c7 = {c3_count}+{c5_count}+{c7_count} = {alpha_1}")
    print(f"    alpha_1 mod 2 = {alpha_1 % 2}")

    # alpha_2 (disjoint pairs)
    alpha_2 = sum(1 for i in range(n3) for j in range(i+1, n3) if not (c3_sets[i] & c3_sets[j]))
    print(f"    alpha_2 (disjoint 3-cycle pairs) = {alpha_2}")

    print(f"    H = 1 + 2*{alpha_1} + 4*{alpha_2} = {1 + 2*alpha_1 + 4*alpha_2}")

    all_data[bits] = {
        'H': H, 'c3_count': c3_count, 'c5_count': c5_count,
        'c7_count': c7_count, 'c3_sets': c3_sets, 'c5_vsets': c5_vsets,
        'orientations': orientations, 'lam': lam, 'ov_counts': ov_counts,
        'signed_ov_sum': signed_ov_sum, 'alpha_2': alpha_2, 'alpha_1': alpha_1
    }


# ========================================================================
# COMPARATIVE ANALYSIS
# ========================================================================
print(f"\n\n{'='*70}")
print("COMPARATIVE ANALYSIS: WHAT DIFFERS BETWEEN c7=8 and c7=9?")
print("=" * 70)

d1 = all_data[4728]  # H=109, c7=8
d2 = all_data[4658]  # H=111, c7=9

print(f"\n  c3_count: {d1['c3_count']} vs {d2['c3_count']} -> delta = {d2['c3_count'] - d1['c3_count']}")
print(f"  c5_count: {d1['c5_count']} vs {d2['c5_count']} -> delta = {d2['c5_count'] - d1['c5_count']}")
print(f"  c7_count: {d1['c7_count']} vs {d2['c7_count']} -> delta = {d2['c7_count'] - d1['c7_count']}")
print(f"  alpha_1:  {d1['alpha_1']} vs {d2['alpha_1']} -> delta = {d2['alpha_1'] - d1['alpha_1']}")
print(f"  alpha_2:  {d1['alpha_2']} vs {d2['alpha_2']} -> delta = {d2['alpha_2'] - d1['alpha_2']}")

print(f"\n  3-cycle overlap distribution:")
for ov in sorted(set(list(d1['ov_counts'].keys()) + list(d2['ov_counts'].keys()))):
    v1 = d1['ov_counts'].get(ov, 0)
    v2 = d2['ov_counts'].get(ov, 0)
    delta = v2 - v1
    marker = " <-- DIFFERENT" if delta != 0 else ""
    print(f"    overlap={ov}: {v1} vs {v2} (delta={delta:+d}){marker}")

print(f"\n  Signed overlap by weight level:")
for ov in sorted(set(list(d1['signed_ov_sum'].keys()) + list(d2['signed_ov_sum'].keys()))):
    v1 = d1['signed_ov_sum'].get(ov, 0)
    v2 = d2['signed_ov_sum'].get(ov, 0)
    delta = v2 - v1
    marker = " <-- DIFFERENT" if delta != 0 else ""
    print(f"    overlap={ov}: {v1:+d} vs {v2:+d} (delta={delta:+d}){marker}")

# Check: which 5-cycle vertex sets differ in directed cycle count?
print(f"\n  5-cycle vertex set comparison:")
c5_data_1 = {}
c5_data_2 = {}
A1 = binary_to_tournament(4728, n)
A2 = binary_to_tournament(4658, n)

for subset in combinations(range(n), 5):
    verts = list(subset)
    nc1 = count_directed_ham_cycles_on_subset(A1, verts)
    nc2 = count_directed_ham_cycles_on_subset(A2, verts)
    c5_data_1[frozenset(subset)] = nc1
    c5_data_2[frozenset(subset)] = nc2
    if nc1 != nc2:
        print(f"    {sorted(subset)}: {nc1} vs {nc2} (delta={nc2-nc1:+d})")

same_5 = sum(1 for s in c5_data_1 if c5_data_1[s] == c5_data_2[s])
diff_5 = sum(1 for s in c5_data_1 if c5_data_1[s] != c5_data_2[s])
print(f"\n  5-cycle vertex sets: {same_5} same, {diff_5} different")

# ========================================================================
# THE DEEP STRUCTURE: Overlap weights as a filtration
# ========================================================================
print(f"\n\n{'='*70}")
print("THE {2,1,0} FILTRATION AND THE HIDDEN DIMENSION")
print("=" * 70)

print("""
The overlap weight W(C_i, C_j) = |V(C_i) intersection V(C_j)| defines
a graded structure on pairs of 3-cycles:

  W = 2: "strongly conflicting" (2 shared vertices out of 3)
         These pairs share an EDGE of the tournament.
         Lambda_{uv} counts exactly how many 3-cycles contain edge (u,v).

  W = 1: "weakly conflicting" (1 shared vertex)
         These pairs share a VERTEX but not an edge.

  W = 0: "independent" (disjoint)
         These pairs are non-adjacent in Omega(T).
         Alpha_2 = #{W=0 pairs}.

The lambda graph captures the W=2 structure completely (lambda_{uv} =
number of 3-cycles containing edge (u,v), and two 3-cycles with overlap 2
share exactly one edge).

THEOREM: Lambda captures W=2 and W=0 counts exactly. So:
  - |{pairs with W=2}| = sum_{u<v} C(lambda_{uv}, 2) [depends only on lambda]
  - |{pairs with W=0}| = alpha_2 [depends only on lambda at n<=8]
  - |{pairs with W=1}| = C(c3,2) - |W=2| - |W=0| [determined by above]

So the TOTAL overlap weight distribution {2,1,0} is lambda-determined!
This means the hidden dimension (c7_dir) is NOT in the overlap weights
of 3-cycles at all. It's in the ORIENTATIONS of 5-cycles.
""")

# Verify this theorem
for bits, label in targets:
    A = binary_to_tournament(bits, n)
    lam, c3_sets = get_labeled_lambda(A, n)
    n3 = len(c3_sets)

    # Count W=2 pairs from lambda
    w2_from_lambda = sum(lam[u][v] * (lam[u][v] - 1) // 2 for u in range(n) for v in range(u+1, n))

    # Count W=2 pairs directly
    w2_direct = 0
    for i in range(n3):
        for j in range(i+1, n3):
            if overlap_weight(c3_sets[i], c3_sets[j]) == 2:
                w2_direct += 1

    print(f"  {label}: W=2 from lambda = {w2_from_lambda}, direct = {w2_direct}, match = {w2_from_lambda == w2_direct}")


# ========================================================================
# THE 5-CYCLE ORIENTATION STRUCTURE
# ========================================================================
print(f"\n\n{'='*70}")
print("5-CYCLE ORIENTATION STRUCTURE: WHERE c7 HIDES")
print("=" * 70)

print("""
Since the {2,1,0} overlap weights on 3-cycles are lambda-determined,
the hidden dimension must live in the 5-cycle or 7-cycle orientations.

For each 5-vertex subset S, the sub-tournament T|_S has some number
of directed 5-cycles. This number can be 0, 1, 2, or 3 (max for 5-cycle
types in a 5-vertex tournament).

The 7-cycle count (c7_dir) is the number of directed Hamiltonian cycles
on the full 7 vertices. It equals 2 * (# undirected HC) if all HCs go
in pairs, but for tournaments, each undirected HC has exactly one
direction as a directed HC... wait, no: in a tournament, each undirected
HC has EXACTLY ONE direction that is a valid directed path.

So c7_dir = number of orientations of the complete graph's Hamiltonian
cycles that are consistent with the tournament.
""")

# Compute the 5-cycle orientation structure more carefully
for bits, label in targets:
    A = binary_to_tournament(bits, n)

    print(f"\n  {label} (bits={bits}):")

    # For each 5-vertex subset, compute: c5 count and c3 count
    for subset in combinations(range(n), 5):
        verts = list(subset)
        c5 = count_directed_ham_cycles_on_subset(A, verts)
        c3 = sum(count_directed_ham_cycles_on_subset(A, list(sub3))
                 for sub3 in combinations(verts, 3))
        # Score sequence of sub-tournament
        scores = sorted(sum(A[verts[i]][verts[j]] for j in range(5) if i != j) for i in range(5))
        print(f"    {sorted(subset)}: c3={c3}, c5={c5}, scores={scores}")

# ========================================================================
# DEEPER: The CHIRALITY of 5-cycle orientations
# ========================================================================
print(f"\n\n{'='*70}")
print("5-CYCLE CHIRALITY: orientation x multiplicity")
print("=" * 70)

for bits, label in targets:
    A = binary_to_tournament(bits, n)

    print(f"\n  {label}:")

    chirality_sum = 0
    for subset in combinations(range(n), 5):
        verts = list(subset)
        c5 = count_directed_ham_cycles_on_subset(A, verts)

        if c5 == 0:
            continue

        # Compute "canonical orientation" = direction following vertex order
        # A 5-cycle on {a,b,c,d,e} (sorted) has canonical direction a->b->c->d->e->a
        # Check if it goes this way or the reverse
        a, b, c, d, e = sorted(subset)

        # Count forward vs backward 5-cycles
        # This is subtle: there are multiple Hamiltonian cycles on 5 vertices
        # (up to 12 directed = 6 undirected per graph, but tournament limits this)

        # Actually for a tournament on 5 vertices, we just count directed HC directly
        # The "chirality" question is: how many go in each direction?
        # But "direction" is not well-defined for k>3...

        # Better approach: use the PARITY of the permutation
        # A directed HC = a permutation cycle. Its sign is (-1)^{k-1}.
        # For k=5: sign = (-1)^4 = +1 (always even). So sign doesn't help.

        # Instead: use the DETERMINANT approach.
        # For each directed HC visiting v0->v1->v2->v3->v4->v0,
        # compute the sign of the permutation (v0,v1,v2,v3,v4) relative to sorted order.

        # Find all directed HCs
        start = 0
        dp = {}
        dp[(1 << start, start)] = [[start]]
        for mask in range(1, 1 << 5):
            if not (mask & (1 << start)):
                continue
            for v in range(5):
                if not (mask & (1 << v)):
                    continue
                key = (mask, v)
                if key not in dp:
                    continue
                for w in range(5):
                    if mask & (1 << w):
                        continue
                    if A[verts[v]][verts[w]]:
                        nk = (mask | (1 << w), w)
                        if nk not in dp:
                            dp[nk] = []
                        for path in dp[key]:
                            dp[nk].append(path + [w])

        full = (1 << 5) - 1
        hc_perms = []
        for v in range(1, 5):
            key = (full, v)
            if key in dp and A[verts[v]][verts[start]]:
                for path in dp[key]:
                    hc_perms.append(tuple(path))

        # Compute sign of each HC permutation
        signs = []
        for perm in hc_perms:
            # perm is a sequence of indices into verts
            # Sign = sign of the permutation
            inversions = sum(1 for a2 in range(5) for b2 in range(a2+1, 5) if perm[a2] > perm[b2])
            sign = (-1) ** inversions
            signs.append(sign)

        total_sign = sum(signs)
        chirality_sum += total_sign

    print(f"    5-cycle chirality sum (perm sign): {chirality_sum}")


# ========================================================================
# THE KEY INSIGHT: Signed adjacency determinant
# ========================================================================
print(f"\n\n{'='*70}")
print("SIGNED ADJACENCY STRUCTURE: det and permanent")
print("=" * 70)

for bits, label in targets:
    A = binary_to_tournament(bits, n)

    # Skew-adjacency matrix: S[i][j] = A[i][j] - A[j][i]
    S = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]

    # Note: full det computed below after det_small is defined


def det_small(M, n):
    """Determinant via Leibniz formula for n <= 7."""
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]

    result = 0
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [[M[r][c] for c in range(n) if c != j] for r in range(1, n)]
        cofactor = ((-1) ** j) * M[0][j] * det_small(minor, n-1)
        result += cofactor
    return result


# Redo with proper determinant
print(f"\n\n{'='*70}")
print("SIGNED ADJACENCY DETERMINANTS")
print("=" * 70)

for bits, label in targets:
    A = binary_to_tournament(bits, n)

    # Skew-adjacency matrix
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

    print(f"\n  {label}:")

    # 6x6 principal minors
    for del_v in range(n):
        remaining = [i for i in range(n) if i != del_v]
        S6 = [[S[remaining[i]][remaining[j]] for j in range(6)] for i in range(6)]
        det6 = det_small(S6, 6)
        # Pfaffian^2 = det for skew-symmetric => Pfaffian = sqrt(|det|)
        pf_sq = det6
        print(f"    det(S del {del_v}) = {det6}, sqrt = {abs(det6)**0.5:.1f}")

    # det(I+A)
    IA = [[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        IA[i][i] = 1
    det_IA = det_small(IA, n)
    print(f"    det(I+A) = {det_IA}")

    # det(I+2A) -- relates to I(Omega, 2)?
    I2A = [[2*A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        I2A[i][i] = 1
    det_I2A = det_small(I2A, n)
    print(f"    det(I+2A) = {det_I2A}")


print(f"\n\n{'='*70}")
print("DONE.")
print("=" * 70)
