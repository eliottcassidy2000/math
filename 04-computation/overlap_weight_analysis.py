#!/usr/bin/env python3
"""
overlap_weight_analysis.py -- Deep analysis of cycle overlap structure

The conflict graph Omega(T) determines H(T) = I(Omega(T), 2).
This script analyzes WHY Interval has more vertex-disjoint cycle pairs
despite having fewer total cycles.

Key analyses:
1. Full overlap weight matrix W[i,j] = |V(C_i) ∩ V(C_j)| and its spectrum
2. 3-cycle anatomy: which triples form cycles, disjointness by arc gaps
3. Additive combinatorics: relation to sum-free sets, additive energy
4. Higher-order independence: full alpha_j decomposition
5. Fractional independence / weighted chromatic structure
6. Neighborhood concentration: how do Omega(T) degrees differ?

Author: kind-pasteur-2026-03-12-S59
"""

import time
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    """Build adjacency matrix for circulant tournament on Z_p with connection set S."""
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def enumerate_directed_cycles(A, p, max_k=None):
    """Enumerate all directed odd cycles, returning list of (frozenset(vertices), length).

    For each vertex subset of size k, count Hamiltonian directed cycles
    within the induced sub-tournament. Each directed cycle on k vertices
    appears once per vertex set (we count directed cycles, not vertex sets).
    """
    if max_k is None:
        max_k = p
    cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            # Count directed Hamiltonian cycles on this k-vertex subset
            verts = list(subset)
            n_cyc = count_directed_ham_cycles(A, verts)
            for _ in range(n_cyc):
                cycles.append((frozenset(subset), k))
    return cycles


def count_directed_ham_cycles(A, verts):
    """Count directed Hamiltonian cycles on vertex subset using DP."""
    k = len(verts)
    if k == 3:
        a, b, c = verts
        fwd = A[a][b] * A[b][c] * A[c][a]
        bwd = A[a][c] * A[c][b] * A[b][a]
        return fwd + bwd

    # General Held-Karp for directed Ham cycles starting and returning to verts[0]
    start = 0  # index within verts
    # dp[mask][v] = # directed paths from start through mask ending at v
    dp = {}
    dp[(1 << start, start)] = 1

    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt

    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def compute_overlap_matrix(cycles):
    """Compute the full overlap weight matrix W[i,j] = |V(C_i) ∩ V(C_j)|."""
    n = len(cycles)
    W = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            ov = len(cycles[i][0] & cycles[j][0])
            W[i][j] = ov
            W[j][i] = ov
    return W


def omega_adjacency(cycles):
    """Build the adjacency matrix of the conflict graph Omega(T).
    Two cycles are adjacent iff they share a vertex."""
    n = len(cycles)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if cycles[i][0] & cycles[j][0]:
                adj[i][j] = 1
                adj[j][i] = 1
    return adj


def count_independent_sets(adj, n):
    """Count all independent sets by size (alpha_j) using inclusion-exclusion on bitmask.
    Only feasible for small n (say n <= 25)."""
    alpha = [0] * (n + 1)
    # Use backtracking with bitmask for speed
    # independent set = subset of vertices with no edges between them

    # For each vertex, precompute its neighbor bitmask
    nbr = [0] * n
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                nbr[i] |= (1 << j)

    # Enumerate all independent sets via backtracking
    def backtrack(v, mask, size):
        alpha[size] += 1
        for w in range(v + 1, n):
            if not (mask & (1 << w)):  # w is not excluded
                # Add w to independent set
                new_mask = mask | nbr[w]
                backtrack(w + 1, new_mask, size + 1)

    backtrack(-1, 0, 0)
    return alpha


def count_independent_sets_large(adj, n, max_size=4):
    """Count independent sets up to size max_size for larger graphs.
    Uses explicit enumeration."""
    # alpha_0 = 1
    alpha = [0] * (max_size + 1)
    alpha[0] = 1

    # alpha_1 = n
    alpha[1] = n

    # alpha_2: count non-adjacent pairs
    if max_size >= 2:
        count2 = 0
        for i in range(n):
            for j in range(i+1, n):
                if not adj[i][j]:
                    count2 += 1
        alpha[2] = count2

    # alpha_3: count non-adjacent triples
    if max_size >= 3:
        # Build neighbor set for each vertex
        nbr_set = [set() for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if adj[i][j]:
                    nbr_set[i].add(j)

        count3 = 0
        for i in range(n):
            for j in range(i+1, n):
                if adj[i][j]:
                    continue
                for k in range(j+1, n):
                    if not adj[i][k] and not adj[j][k]:
                        count3 += 1
        alpha[3] = count3

    # alpha_4: count non-adjacent quadruples
    if max_size >= 4:
        count4 = 0
        for i in range(n):
            for j in range(i+1, n):
                if adj[i][j]:
                    continue
                for k in range(j+1, n):
                    if adj[i][k] or adj[j][k]:
                        continue
                    for l in range(k+1, n):
                        if not adj[i][l] and not adj[j][l] and not adj[k][l]:
                            count4 += 1
        alpha[4] = count4

    return alpha


def three_cycle_gap_analysis(A, p, S):
    """Analyze the "gap structure" of 3-cycles.

    A 3-cycle {a,b,c} with a->b->c->a has gaps (b-a, c-b, a-c) mod p.
    These gaps all lie in S. This gives a triple of elements of S summing to 0.
    """
    S_set = set(S)
    m = (p - 1) // 2

    # Find all 3-cycles as ordered triples (a,b,c) with a<b<c as vertex set
    cycles_3 = []
    gap_triples = []

    for a, b, c in combinations(range(p), 3):
        # Check a->b->c->a
        if A[a][b] and A[b][c] and A[c][a]:
            g1 = (b - a) % p
            g2 = (c - b) % p
            g3 = (a - c) % p
            cycles_3.append(((a, b, c), 'fwd'))
            gap_triples.append(tuple(sorted([g1, g2, g3])))
        # Check a->c->b->a
        if A[a][c] and A[c][b] and A[b][a]:
            g1 = (c - a) % p
            g2 = (b - c) % p
            g3 = (a - b) % p
            cycles_3.append(((a, c, b), 'bwd'))
            gap_triples.append(tuple(sorted([g1, g2, g3])))

    # Count gap triple types
    gap_counts = defaultdict(int)
    for gt in gap_triples:
        gap_counts[gt] += 1

    return cycles_3, gap_triples, gap_counts


def disjoint_pair_gap_analysis(A, p, S, cycles):
    """For each pair of disjoint 3-cycles, analyze the "complementary gap" structure.

    Two 3-cycles {a,b,c} and {d,e,f} are disjoint iff they use 6 distinct vertices.
    At p=7, this leaves exactly 1 unused vertex.

    Question: what is the relationship between the gap structure of the two cycles
    and whether they are disjoint?
    """
    c3_cycles = [(fs, k) for fs, k in cycles if k == 3]

    disjoint_pairs = []
    adjacent_pairs = []

    for i in range(len(c3_cycles)):
        for j in range(i+1, len(c3_cycles)):
            if not (c3_cycles[i][0] & c3_cycles[j][0]):
                disjoint_pairs.append((i, j))
            else:
                ov = len(c3_cycles[i][0] & c3_cycles[j][0])
                adjacent_pairs.append((i, j, ov))

    return c3_cycles, disjoint_pairs, adjacent_pairs


def vertex_degree_in_omega(adj, n):
    """Compute degree distribution of Omega(T)."""
    degrees = [sum(adj[i]) for i in range(n)]
    return degrees


def omega_clique_analysis(adj, n, max_clique=5):
    """Analyze clique structure of Omega(T).
    Cliques in Omega = sets of mutually conflicting cycles."""
    cliques = {1: n}

    if max_clique >= 2:
        c2 = 0
        for i in range(n):
            for j in range(i+1, n):
                if adj[i][j]:
                    c2 += 1
        cliques[2] = c2

    if max_clique >= 3 and n <= 200:
        c3 = 0
        for i in range(n):
            for j in range(i+1, n):
                if not adj[i][j]:
                    continue
                for k in range(j+1, n):
                    if adj[i][k] and adj[j][k]:
                        c3 += 1
        cliques[3] = c3

    return cliques


def analyze_overlap_by_vertex_pair(cycles, p):
    """For each pair of Z_p vertices (u,v), count how many cycles contain both.

    This gives a "co-occurrence matrix" on Z_p. For circulant tournaments,
    co-occurrence(u,v) depends only on (v-u) mod p.
    """
    co_occur = [[0]*p for _ in range(p)]
    for fs, k in cycles:
        verts = list(fs)
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                co_occur[verts[i]][verts[j]] += 1
                co_occur[verts[j]][verts[i]] += 1

    # For circulant: co-occurrence should depend only on gap
    gap_co = [0] * p
    for d in range(1, p):
        gap_co[d] = co_occur[0][d]  # by circulant symmetry

    return co_occur, gap_co


def sum_free_analysis(S, p):
    """Analyze the sum-free structure of S.

    A set A is sum-free if (A+A) ∩ A = ∅.
    For tournament connection sets, S has size m=(p-1)/2.

    The number of 3-cycles relates to |{(a,b) in S x S : a+b in S}|.
    """
    S_set = set(S)
    m = len(S)

    # Count sum pairs: (a,b) in S x S with a+b in S
    sum_pairs = 0
    sum_triples = []  # (a, b, a+b) triples
    for a in S:
        for b in S:
            c = (a + b) % p
            if c in S_set:
                sum_pairs += 1
                sum_triples.append((a, b, c))

    # Sum triples give directed 3-cycles: vertex v, v+a, v+a+b = v+c'
    # where the cycle goes v -> v+a -> v+(a+b) -> v if a,b,(p-a-b) all in S
    # Actually: a 3-cycle {v, v+g1, v+g1+g2} needs g1 in S, g2 in S, and
    # -(g1+g2) = p-g1-g2 in S.

    # Count zero-sum triples in S: (g1, g2, g3) with g1+g2+g3 = 0, all in S
    zero_sum = 0
    for a in S:
        for b in S:
            c = (p - a - b) % p
            if c in S_set:
                zero_sum += 1

    # Additive energy
    energy = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    energy += 1

    # Fourier concentration (simplified)
    import cmath
    omega = cmath.exp(2j * cmath.pi / p)
    fourier = []
    for t in range(p):
        val = sum(omega ** (t * s) for s in S)
        fourier.append(abs(val)**2)

    return {
        'sum_pairs': sum_pairs,
        'zero_sum_triples': zero_sum,
        'additive_energy': energy,
        'fourier_spectrum': fourier,
        'sum_triples': sum_triples
    }


def shared_difference_analysis(cycles, p):
    """For each pair of cycles, compute the set of "differences" used.

    A cycle C = {v0, v1, ..., v_{k-1}} uses arcs with gaps
    {v1-v0, v2-v1, ..., v0-v_{k-1}} mod p.

    Two cycles are more likely to conflict if they use similar gaps
    (because similar gaps => similar vertex placement on Z_p).
    """
    # For 3-cycles, the gap set is a 3-element subset of S
    # For two 3-cycles to be disjoint, their gap sets must be "compatible"
    # in the sense that the vertex placements don't overlap

    c3_cycles = [(fs, k) for fs, k in cycles if k == 3]

    # For each 3-cycle, extract its gap set
    gap_sets = []
    for fs, k in c3_cycles:
        verts = sorted(fs)
        # All pairwise differences
        diffs = set()
        for a in verts:
            for b in verts:
                if a != b:
                    diffs.add((b - a) % p)
        gap_sets.append(diffs)

    # Count shared differences between pairs
    shared_diff_vs_overlap = defaultdict(list)
    for i in range(len(c3_cycles)):
        for j in range(i+1, len(c3_cycles)):
            ov = len(c3_cycles[i][0] & c3_cycles[j][0])
            sd = len(gap_sets[i] & gap_sets[j])
            shared_diff_vs_overlap[ov].append(sd)

    return shared_diff_vs_overlap


def c3_only_overlap_analysis(A, p, S, name):
    """Detailed analysis restricted to 3-cycles only — tractable for any p."""
    print(f"\n  --- 3-CYCLE FOCUSED ANALYSIS ---")

    S_set = set(S)
    m = (p - 1) // 2

    # Enumerate 3-cycles
    c3_list = []
    for a, b, c in combinations(range(p), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            c3_list.append(frozenset([a, b, c]))
        if A[a][c] and A[c][b] and A[b][a]:
            c3_list.append(frozenset([a, b, c]))

    # Deduplicate (each vertex set can have 0, 1, or 2 directed 3-cycles)
    # Actually for tournaments, a 3-vertex set has EXACTLY 1 directed 3-cycle or 0
    # No: 3-vertex tournament has either 1 Ham cycle (both dirs) or 0
    # For our listing: each {a,b,c} with a cycle contributes 2 directed cycles (fwd+bwd)
    # But we only care about vertex sets for overlap. Let's count directed and vertex sets.
    c3_vertex_sets = list(set(c3_list))
    c3_directed = len(c3_list)

    print(f"    3-cycle vertex sets: {len(c3_vertex_sets)}")
    print(f"    Directed 3-cycles: {c3_directed}")

    # Overlap matrix for 3-cycle vertex sets
    n3 = len(c3_vertex_sets)
    disjoint_3_pairs = 0
    overlap_1_pairs = 0
    overlap_2_pairs = 0
    overlap_3_pairs = 0

    for i in range(n3):
        for j in range(i + 1, n3):
            ov = len(c3_vertex_sets[i] & c3_vertex_sets[j])
            if ov == 0:
                disjoint_3_pairs += 1
            elif ov == 1:
                overlap_1_pairs += 1
            elif ov == 2:
                overlap_2_pairs += 1
            else:
                overlap_3_pairs += 1

    total_3_pairs = n3 * (n3 - 1) // 2
    print(f"    C(c3, 2) = {total_3_pairs}")
    print(f"    Disjoint 3-3 pairs: {disjoint_3_pairs} ({100*disjoint_3_pairs/total_3_pairs:.2f}%)")
    print(f"    Overlap=1 pairs:    {overlap_1_pairs} ({100*overlap_1_pairs/total_3_pairs:.2f}%)")
    print(f"    Overlap=2 pairs:    {overlap_2_pairs} ({100*overlap_2_pairs/total_3_pairs:.2f}%)")
    print(f"    Overlap=3 pairs:    {overlap_3_pairs} ({100*overlap_3_pairs/total_3_pairs:.2f}%)")

    # Disjoint triple count (3 mutually disjoint 3-cycles = 9 vertices needed)
    if p >= 9:
        disjoint_triples = 0
        for i in range(n3):
            for j in range(i + 1, n3):
                if c3_vertex_sets[i] & c3_vertex_sets[j]:
                    continue
                for k in range(j + 1, n3):
                    if not (c3_vertex_sets[i] & c3_vertex_sets[k]) and \
                       not (c3_vertex_sets[j] & c3_vertex_sets[k]):
                        disjoint_triples += 1
        print(f"    Mutually disjoint 3-3-3 triples: {disjoint_triples}")

    return c3_vertex_sets, disjoint_3_pairs


def main():
    print("=" * 70)
    print("DEEP OVERLAP WEIGHT ANALYSIS")
    print("=" * 70)

    for p in [7, 11, 13, 19]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))
        # Always compare both Paley and Interval
        tournaments = [("Paley", S_qr), ("Interval", S_int)]

        for name, S in tournaments:
            A = build_adj(p, S)

            print(f"\n{'='*70}")
            print(f"p={p}, {name}, S={S}")
            print(f"{'='*70}")

            # ====== CYCLE ENUMERATION ======
            t0 = time.time()
            # For p >= 13, limit to k=7 to avoid huge cycle counts
            max_k = p if p <= 11 else 7
            cycles = enumerate_directed_cycles(A, p, max_k=max_k)
            t1 = time.time()
            print(f"\n  Total cycles (k<={max_k}): {len(cycles)} (enumerated in {t1-t0:.1f}s)")

            # Count by length
            by_length = defaultdict(int)
            for fs, k in cycles:
                by_length[k] += 1
            for k in sorted(by_length):
                print(f"    k={k}: {by_length[k]} directed cycles")

            n_cyc = len(cycles)

            # ====== FULL ALPHA DECOMPOSITION (only feasible for small n_cyc) ======
            if n_cyc <= 200:
                print(f"\n  --- Alpha decomposition ---")
                adj = omega_adjacency(cycles)

                if n_cyc <= 22:
                    alpha = count_independent_sets(adj, n_cyc)
                    max_j = max(j for j in range(len(alpha)) if alpha[j] > 0)
                    for j in range(max_j + 1):
                        print(f"    alpha_{j} = {alpha[j]}")
                    H_check = sum(alpha[j] * (2**j) for j in range(len(alpha)))
                    print(f"    I(Omega, 2) = {H_check}")
                else:
                    alpha = count_independent_sets_large(adj, n_cyc, max_size=4)
                    for j in range(len(alpha)):
                        if alpha[j] > 0:
                            print(f"    alpha_{j} = {alpha[j]}")

                # ====== OVERLAP WEIGHT SPECTRUM ======
                print(f"\n  --- Overlap weight spectrum ---")
                overlap_dist = defaultdict(int)
                for i in range(n_cyc):
                    for j in range(i+1, n_cyc):
                        ov = len(cycles[i][0] & cycles[j][0])
                        overlap_dist[ov] += 1

                total_pairs = n_cyc * (n_cyc - 1) // 2
                for ov in sorted(overlap_dist):
                    pct = 100 * overlap_dist[ov] / total_pairs
                    print(f"    overlap={ov}: {overlap_dist[ov]:>8d} pairs ({pct:6.2f}%)")

                mean_ov = sum(ov * cnt for ov, cnt in overlap_dist.items()) / total_pairs
                var_ov = sum((ov - mean_ov)**2 * cnt for ov, cnt in overlap_dist.items()) / total_pairs
                print(f"    Mean overlap = {mean_ov:.3f}")
                print(f"    Var overlap = {var_ov:.3f}")
                print(f"    Std overlap = {var_ov**0.5:.3f}")

                # ====== OMEGA DEGREE DISTRIBUTION ======
                print(f"\n  --- Omega degree distribution ---")
                degrees = vertex_degree_in_omega(adj, n_cyc)
                deg_dist = defaultdict(int)
                for d in degrees:
                    deg_dist[d] += 1
                avg_deg = sum(degrees) / n_cyc if n_cyc > 0 else 0
                max_deg = max(degrees) if degrees else 0
                min_deg = min(degrees) if degrees else 0
                print(f"    Avg degree = {avg_deg:.1f}")
                print(f"    Min degree = {min_deg}")
                print(f"    Max degree = {max_deg}")
                print(f"    Degree values: {dict(sorted(deg_dist.items()))}")

                # ====== CLIQUE ANALYSIS ======
                if n_cyc <= 500:
                    print(f"\n  --- Omega clique counts ---")
                    cliques = omega_clique_analysis(adj, n_cyc, max_clique=3)
                    for k_cl, cnt in sorted(cliques.items()):
                        print(f"    {k_cl}-cliques: {cnt}")

                # ====== SHARED DIFFERENCES vs OVERLAP ======
                if n_cyc <= 1000:
                    print(f"\n  --- Shared differences vs vertex overlap ---")
                    sd_vs_ov = shared_difference_analysis(cycles, p)
                    for ov in sorted(sd_vs_ov):
                        vals = sd_vs_ov[ov]
                        if vals:
                            avg_sd = sum(vals) / len(vals)
                            print(f"    overlap={ov}: avg shared diffs = {avg_sd:.2f} (n={len(vals)})")

            else:
                print(f"\n  [Skipping full Omega analysis: {n_cyc} cycles too large]")

            # ====== 3-CYCLE FOCUSED ANALYSIS (always feasible) ======
            c3_sets, disj_3 = c3_only_overlap_analysis(A, p, S, name)

            # ====== 3-CYCLE GAP ANATOMY ======
            print(f"\n  --- 3-cycle gap anatomy ---")
            c3_info, gap_triples, gap_counts = three_cycle_gap_analysis(A, p, S)
            print(f"    {len(c3_info)} directed 3-cycles")
            print(f"    Gap triple types (sorted):")
            for gt, cnt in sorted(gap_counts.items(), key=lambda x: -x[1]):
                print(f"      {gt}: {cnt} cycles")

            # ====== CO-OCCURRENCE ON Z_p (using 3-cycles only for large p) ======
            if n_cyc <= 200:
                print(f"\n  --- Co-occurrence by gap (ALL cycles) ---")
                co_occur, gap_co = analyze_overlap_by_vertex_pair(cycles, p)
            else:
                print(f"\n  --- Co-occurrence by gap (3-cycles only) ---")
                c3_as_cycles = [(fs, 3) for fs in c3_sets]
                co_occur, gap_co = analyze_overlap_by_vertex_pair(c3_as_cycles, p)

            for d in range(1, (p+1)//2 + 1):
                print(f"    gap {d}: {gap_co[d]} shared cycles (gap {p-d}: {gap_co[p-d]})")

            # Co-occurrence variance
            co_vals = [gap_co[d] for d in range(1, p)]
            co_mean = sum(co_vals) / len(co_vals)
            co_var = sum((c - co_mean)**2 for c in co_vals) / len(co_vals)
            print(f"    Co-occurrence mean = {co_mean:.2f}, variance = {co_var:.2f}")

            # ====== SUM-FREE ANALYSIS ======
            print(f"\n  --- Additive structure of S ---")
            sf = sum_free_analysis(S, p)
            print(f"    Sum pairs in S: {sf['sum_pairs']}")
            print(f"    Zero-sum triples: {sf['zero_sum_triples']}")
            print(f"    Additive energy: {sf['additive_energy']}")
            print(f"    c3 = zero_sum / 3 = {sf['zero_sum_triples'] // 3}")

            # Fourier spectrum
            fourier = sf['fourier_spectrum']
            print(f"    Fourier |S_hat|^2: ", end="")
            for t in range(p):
                print(f"{fourier[t]:.1f}", end=" ")
            print()
            L4 = sum(f**2 for f in fourier)
            L2 = sum(f for f in fourier)
            print(f"    Fourier L4/L2^2 = {L4/L2**2:.6f}")

            # ====== DISJOINT PAIR VERTEX ANALYSIS ======
            if p <= 11 and c3_sets:
                print(f"\n  --- Disjoint pair vertex analysis ---")
                # Count disjoint pairs from c3_sets
                n3 = len(c3_sets)
                disj_pairs_list = []
                for i in range(n3):
                    for j in range(i+1, n3):
                        if not (c3_sets[i] & c3_sets[j]):
                            disj_pairs_list.append((i, j))

                if p == 7 and disj_pairs_list:
                    print(f"\n    Free vertex distribution (disjoint 3-3 pairs at p=7):")
                    free_verts = []
                    for i, j in disj_pairs_list:
                        used = c3_sets[i] | c3_sets[j]
                        free = set(range(p)) - used
                        free_verts.extend(list(free))
                    free_dist = defaultdict(int)
                    for v in free_verts:
                        free_dist[v] += 1
                    print(f"    {dict(sorted(free_dist.items()))}")
                    if len(set(free_dist.values())) == 1:
                        print(f"    UNIFORM: each vertex is free in {list(free_dist.values())[0]} pairs")

    # ====== COMPARATIVE SUMMARY ======
    print(f"\n{'='*70}")
    print("COMPARATIVE SUMMARY")
    print("=" * 70)

    for p in [7, 11, 13, 19]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))
        tourn_list = [("Paley", S_qr), ("Interval", S_int)]

        for name, S in tourn_list:
            A = build_adj(p, S)

            # 3-cycle vertex sets
            c3_sets_loc = []
            for a, b, c in combinations(range(p), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3_sets_loc.append(frozenset([a, b, c]))
            c3_sets_loc = list(set(c3_sets_loc))
            n3 = len(c3_sets_loc)

            disj_3 = 0
            for i in range(n3):
                for j in range(i + 1, n3):
                    if not (c3_sets_loc[i] & c3_sets_loc[j]):
                        disj_3 += 1
            total_3_pairs = n3 * (n3 - 1) // 2

            # Additive analysis
            sf = sum_free_analysis(S, p)
            fourier = sf['fourier_spectrum']
            L4 = sum(f**2 for f in fourier)
            L2 = sum(f for f in fourier)

            # Co-occurrence
            c3_as_cycles = [(fs, 3) for fs in c3_sets_loc]
            _, gap_co = analyze_overlap_by_vertex_pair(c3_as_cycles, p)
            co_vals = [gap_co[d] for d in range(1, p)]
            co_var = sum((c - sum(co_vals)/len(co_vals))**2 for c in co_vals) / len(co_vals)

            print(f"\n  p={p}, {name}, S={S}:")
            print(f"    c3 vertex sets = {n3}")
            print(f"    Disjoint 3-3 pairs = {disj_3}")
            if total_3_pairs > 0:
                print(f"    Disjointness ratio = {disj_3/total_3_pairs:.6f}")
            print(f"    Zero-sum triples = {sf['zero_sum_triples']}")
            print(f"    Additive energy = {sf['additive_energy']}")
            print(f"    Fourier L4/L2^2 = {L4/L2**2:.6f}")
            print(f"    Co-occurrence var = {co_var:.2f}")


if __name__ == '__main__':
    main()
