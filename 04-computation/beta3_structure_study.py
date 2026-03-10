"""
beta3_structure_study.py — Deep structural analysis of beta_3=2 tournaments at n=8

Extracts the beta_3=2 tournaments found in claims_n8_extended.py (seed 12345, 5000 trials)
and studies their structural properties:
1. Score sequences, cycle counts, regularity
2. Omega(T) graph structure (independence polynomial, chromatic number)
3. Deletion patterns — which vertices are good/bad?
4. H_3 generator structure — are the two generators related?
5. What distinguishes beta_3=2 from beta_3=1 tournaments?

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp, adj_to_bits,
    enumerate_all_allowed, boundary_faces,
    _build_constraint_matrix, _gauss_nullbasis_modp, _gauss_rank_np,
    RANK_PRIME
)

PRIME = RANK_PRIME


def count_directed_cycles(A, n, max_len=None):
    """Count directed cycles of each length in tournament A."""
    if max_len is None:
        max_len = n
    cycle_counts = {}
    for k in range(3, max_len + 1):
        count = 0
        for subset in __import__('itertools').combinations(range(n), k):
            sub = [[A[subset[i]][subset[j]] for j in range(k)] for i in range(k)]
            # Count directed cycles on this vertex set using DFS
            cycles = count_cycles_on_subset(sub, k)
            count += cycles
        cycle_counts[k] = count
    return cycle_counts


def count_cycles_on_subset(sub, k):
    """Count directed Hamiltonian cycles on a k-vertex tournament (adjacency sub)."""
    if k < 3:
        return 0
    count = 0
    # Fix vertex 0 as start to avoid counting each cycle k times
    # Enumerate permutations of remaining vertices
    from itertools import permutations
    remaining = list(range(1, k))
    for perm in permutations(remaining):
        path = [0] + list(perm)
        valid = True
        for i in range(k):
            if sub[path[i]][path[(i+1) % k]] != 1:
                valid = False
                break
        if valid:
            count += 1
    # Each cycle is counted (k-1)! / k times by fixing vertex 0
    # Actually fixing vertex 0 counts each cycle exactly once
    # Wait: fixing vertex 0 and going through permutations of others
    # counts each Hamiltonian cycle once (since we broke the cycle at 0)
    return count


def omega_graph_info(A, n):
    """Compute Omega(T) structure: cycle list, adjacency, independence polynomial."""
    from itertools import combinations

    cycles_by_len = {}
    all_cycles = []

    for k in range(3, n + 1, 2):  # odd cycles only
        for subset in combinations(range(n), k):
            sub = [[A[subset[i]][subset[j]] for j in range(k)] for i in range(k)]
            nc = count_cycles_on_subset(sub, k)
            if nc > 0:
                for _ in range(nc):
                    all_cycles.append(set(subset))
                cycles_by_len.setdefault(k, []).append(set(subset))

    # Wait, this double-counts: multiple directed cycles on same vertex set
    # For Omega(T), each directed odd cycle is a separate vertex
    # But two cycles on the SAME vertex set are always adjacent (share all vertices)
    # So for independence polynomial, we need to track vertex sets + cycle count

    # Actually for Omega(T), the conflict graph has one vertex per directed odd cycle
    # and adjacency = share a vertex. Two directed cycles on the same vertex set
    # are trivially adjacent. So the IP treats each cycle as a separate node.

    # For simplicity, let's compute alpha_k (max independent sets of k cycles)
    # by vertex-disjoint cycle collections

    # Collect all vertex sets that support odd cycles, with multiplicity
    cycle_vertex_sets = []
    cycle_multiplicities = {}

    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            sub = [[A[subset[i]][subset[j]] for j in range(k)] for i in range(k)]
            nc = count_cycles_on_subset(sub, k)
            if nc > 0:
                cycle_vertex_sets.append((frozenset(subset), nc))
                cycle_multiplicities[frozenset(subset)] = nc

    # Independent sets = vertex-disjoint collections of cycle vertex sets
    # Alpha_k = sum over such collections of product of multiplicities
    alpha = [1]  # alpha_0 = 1

    # alpha_1 = sum of multiplicities
    alpha1 = sum(m for _, m in cycle_vertex_sets)
    alpha.append(alpha1)

    # alpha_2 = sum over disjoint pairs
    unique_sets = list(cycle_multiplicities.keys())
    alpha2 = 0
    for i in range(len(unique_sets)):
        for j in range(i + 1, len(unique_sets)):
            if unique_sets[i].isdisjoint(unique_sets[j]):
                alpha2 += cycle_multiplicities[unique_sets[i]] * cycle_multiplicities[unique_sets[j]]
    alpha.append(alpha2)

    # alpha_3 for completeness
    alpha3 = 0
    for i in range(len(unique_sets)):
        for j in range(i + 1, len(unique_sets)):
            if not unique_sets[i].isdisjoint(unique_sets[j]):
                continue
            for kk in range(j + 1, len(unique_sets)):
                if unique_sets[kk].isdisjoint(unique_sets[i]) and unique_sets[kk].isdisjoint(unique_sets[j]):
                    alpha3 += (cycle_multiplicities[unique_sets[i]] *
                              cycle_multiplicities[unique_sets[j]] *
                              cycle_multiplicities[unique_sets[kk]])
    alpha.append(alpha3)

    return {
        'cycle_vertex_sets': cycle_multiplicities,
        'alpha': alpha,
        'total_odd_cycles': alpha1,
        'I_at_2': sum(a * (2**k) for k, a in enumerate(alpha)),
    }


def analyze_tournament(A, n, label=""):
    """Full structural analysis of one tournament."""
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")

    # Score sequence
    scores = sorted([int(sum(A[i])) for i in range(n)])
    print(f"  Score sequence: {scores}")
    print(f"  Score variance: {np.var(scores):.3f}")

    # Betti numbers (full)
    res = full_chain_complex_modp(A, n, max_p=7)
    bettis = res['bettis']
    profile = tuple(bettis.get(p, 0) for p in range(n))
    print(f"  Betti profile: {profile}")

    # Hamiltonian path count via OCF
    # (Quick: just get H from Omega)

    # Cycle counts
    from itertools import combinations
    c3, c5, c7 = 0, 0, 0
    for subset in combinations(range(n), 3):
        sub = [[A[subset[i]][subset[j]] for j in range(3)] for i in range(3)]
        if count_cycles_on_subset(sub, 3) > 0:
            c3 += 1
    for subset in combinations(range(n), 5):
        sub = [[A[subset[i]][subset[j]] for j in range(5)] for i in range(5)]
        c5 += count_cycles_on_subset(sub, 5)
    for subset in combinations(range(n), 7):
        sub = [[A[subset[i]][subset[j]] for j in range(7)] for i in range(7)]
        c7 += count_cycles_on_subset(sub, 7)

    print(f"  Cycle counts: c3={c3}, c5={c5}, c7={c7}")

    # Omega structure
    omega = omega_graph_info(A, n)
    print(f"  Omega: alpha = {omega['alpha']}")
    H_ocf = omega['I_at_2']
    print(f"  H(T) = I(Omega, 2) = {H_ocf}")

    # Deletion analysis
    print(f"\n  Vertex deletion analysis:")
    good_vertices = []
    bad_vertices = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
        res_v = full_chain_complex_modp(A_sub, n-1, max_p=6)
        b3_v = res_v['bettis'].get(3, 0)
        b4_v = res_v['bettis'].get(4, 0)
        bv_profile = tuple(res_v['bettis'].get(p, 0) for p in range(n-1))
        out_deg = int(sum(A[v]))

        tag = "GOOD" if b3_v == 0 else "BAD"
        if b3_v == 0:
            good_vertices.append(v)
        else:
            bad_vertices.append(v)

        print(f"    v={v} (out={out_deg}): b3={b3_v}, b4={b4_v}, betti={bv_profile} [{tag}]")

    print(f"\n  Good vertices: {good_vertices} ({len(good_vertices)}/{n})")
    print(f"  Bad vertices:  {bad_vertices} ({len(bad_vertices)}/{n})")

    # Bits encoding
    bits = adj_to_bits(A, n)
    print(f"  Bits encoding: {bits}")

    return {
        'scores': scores,
        'profile': profile,
        'c3': c3, 'c5': c5, 'c7': c7,
        'omega': omega,
        'H': H_ocf,
        'good_vertices': good_vertices,
        'bad_vertices': bad_vertices,
        'bits': bits,
    }


def main():
    print("=" * 70)
    print("BETA_3=2 STRUCTURAL ANALYSIS AT n=8")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(12345)

    # First pass: find all beta_3=2 tournaments
    print("\n--- Phase 1: Extracting beta_3=2 tournaments ---")
    b3_2_tours = []
    b3_1_examples = []  # Keep a few beta_3=1 for comparison

    t0 = time.time()
    for trial in range(5000):
        A = random_tournament(n, rng)
        res = full_chain_complex_modp(A, n, max_p=5)
        b3 = res['bettis'].get(3, 0)

        if b3 == 2:
            b3_2_tours.append((trial, A.copy()))
            print(f"  Found beta_3=2 at trial {trial}")
        elif b3 == 1 and len(b3_1_examples) < 3:
            b3_1_examples.append((trial, A.copy()))

        if (trial + 1) % 1000 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/5000, {elapsed:.1f}s, found {len(b3_2_tours)} beta_3=2")

    print(f"\nFound {len(b3_2_tours)} beta_3=2 tournaments in 5000 trials")

    # Phase 2: Deep analysis of each beta_3=2 tournament
    print(f"\n--- Phase 2: Deep structural analysis ---")

    results_b3_2 = []
    for trial, A in b3_2_tours:
        r = analyze_tournament(A, n, label=f"BETA_3=2 (trial {trial})")
        results_b3_2.append(r)

    # Phase 3: Compare with beta_3=1 examples
    print(f"\n--- Phase 3: Comparison with beta_3=1 tournaments ---")

    results_b3_1 = []
    for trial, A in b3_1_examples[:2]:
        r = analyze_tournament(A, n, label=f"BETA_3=1 (trial {trial})")
        results_b3_1.append(r)

    # Phase 4: Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    print(f"\nbeta_3=2 tournaments ({len(results_b3_2)}):")
    for i, r in enumerate(results_b3_2):
        print(f"  [{i}] scores={r['scores']}, c3={r['c3']}, c5={r['c5']}, "
              f"H={r['H']}, alpha={r['omega']['alpha']}, "
              f"#good={len(r['good_vertices'])}, #bad={len(r['bad_vertices'])}")

    print(f"\nbeta_3=1 comparison ({len(results_b3_1)}):")
    for i, r in enumerate(results_b3_1):
        print(f"  [{i}] scores={r['scores']}, c3={r['c3']}, c5={r['c5']}, "
              f"H={r['H']}, alpha={r['omega']['alpha']}, "
              f"#good={len(r['good_vertices'])}, #bad={len(r['bad_vertices'])}")

    # Check for common structural patterns in beta_3=2
    if results_b3_2:
        print("\n--- Structural patterns in beta_3=2 ---")
        all_scores = [tuple(r['scores']) for r in results_b3_2]
        score_counts = {}
        for s in all_scores:
            score_counts[s] = score_counts.get(s, 0) + 1
        print(f"  Score distributions: {score_counts}")

        all_c3 = [r['c3'] for r in results_b3_2]
        all_c5 = [r['c5'] for r in results_b3_2]
        print(f"  c3 range: {min(all_c3)}-{max(all_c3)}")
        print(f"  c5 range: {min(all_c5)}-{max(all_c5)}")

        all_good = [len(r['good_vertices']) for r in results_b3_2]
        all_bad = [len(r['bad_vertices']) for r in results_b3_2]
        print(f"  Good vertices: {min(all_good)}-{max(all_good)}")
        print(f"  Bad vertices:  {min(all_bad)}-{max(all_bad)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
