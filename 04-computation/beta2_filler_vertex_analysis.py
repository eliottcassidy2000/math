#!/usr/bin/env python3
"""
Filler Vertex Analysis for β₁ deletion bound.

Claim: If β₁(T)=0, then Σ_v β₁(T\v) ≤ 3.

For each tournament T with β₁(T)=0:
  - Find critical vertices v where β₁(T\v)=1
  - For each such v, find the 1-cycle z in T\v and the filling 2-chain c in T
  - Analyze which filling paths pass through v
  - Check if critical sets have size ≤ 3 and share structure (e.g., 3-cycles)
"""

import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
import sys

# ============================================================
# Path homology core (from path_homology_v2.py, no validation)
# ============================================================

def enumerate_allowed_paths(A, n, p):
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))
    if p == 0:
        return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space

def compute_homology_detailed(A, n, dim_target):
    """Compute path homology at dimension dim_target, returning
    kernel basis, image basis, and betti number, plus allowed paths and omega basis."""
    max_dim = dim_target + 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    omega = {}
    for p in range(max_dim + 2):
        if p <= max_dim + 1:
            omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    p = dim_target
    dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
    if dim_omega_p == 0:
        return {'betti': 0, 'ker_basis': None, 'im_basis': None,
                'allowed': allowed, 'omega': omega}

    # kernel of ∂_p restricted to Ω_p
    bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
    bd_p_omega = bd_p @ omega[p]

    if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
        U_p, S_p, Vt_p = np.linalg.svd(bd_p_omega, full_matrices=True)
        rank_p = sum(s > 1e-8 for s in S_p)
    else:
        rank_p = 0
        Vt_p = np.eye(dim_omega_p)

    ker_dim = dim_omega_p - rank_p
    # Kernel basis in Ω_p coordinates
    ker_omega_coords = Vt_p[rank_p:].T  # columns
    # Kernel basis in A_p coordinates
    ker_Ap = omega[p] @ ker_omega_coords if ker_dim > 0 else None

    # image of ∂_{p+1}
    p1 = p + 1
    dim_omega_p1 = omega[p1].shape[1] if p1 in omega and omega[p1].ndim == 2 else 0
    if dim_omega_p1 > 0:
        bd_p1 = build_full_boundary_matrix(allowed[p1], allowed[p])
        bd_p1_omega = bd_p1 @ omega[p1]
        U_p1, S_p1, Vt_p1 = np.linalg.svd(bd_p1_omega, full_matrices=True)
        im_rank = sum(s > 1e-8 for s in S_p1)
        # Image basis in A_p coordinates
        im_Ap = U_p1[:, :im_rank] if im_rank > 0 else None
    else:
        im_rank = 0
        im_Ap = None

    betti = max(0, ker_dim - im_rank)

    return {'betti': betti, 'ker_basis': ker_Ap, 'im_basis': im_Ap,
            'allowed': allowed, 'omega': omega, 'ker_dim': ker_dim, 'im_rank': im_rank}

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None:
        max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti.append(0)
            continue
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0
        betti.append(max(0, ker_dim - im_dim))
    return betti


# ============================================================
# Tournament utilities
# ============================================================

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

def delete_vertex(A, n, v):
    """Return adjacency matrix of T\v on n-1 vertices."""
    verts = [i for i in range(n) if i != v]
    nm1 = n - 1
    B = [[0]*nm1 for _ in range(nm1)]
    for i_new, i_old in enumerate(verts):
        for j_new, j_old in enumerate(verts):
            B[i_new][j_new] = A[i_old][j_old]
    return B, verts

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def find_3cycles_through(A, n, vertices):
    """Find 3-cycles that pass through ALL vertices in the given set."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                triple = {i, j, k}
                if not set(vertices).issubset(triple):
                    continue
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i, j, k))
                elif A[j][i] and A[k][j] and A[i][k]:
                    cycles.append((k, j, i))
    return cycles

def score_sequence(A, n):
    return sorted([sum(A[i]) for i in range(n)])

def tournament_id(A, n):
    """Simple canonical ID (not isomorphism-invariant, just for dedup within enumeration)."""
    bits = []
    for i in range(n):
        for j in range(i+1, n):
            bits.append(A[i][j])
    return tuple(bits)


def find_filling_chain(A_full, n, A_del, n_del, v, verts_del, cycle_vec, allowed_1_del, allowed_2_full):
    """Given a 1-cycle z in T\v (as vector over allowed 1-paths of T\v),
    find a 2-chain c in T such that ∂₂(c) restricted to T\v edges equals z.

    Returns: the 2-chain c as a dict {2-path: coefficient}, and info about
    which 2-paths pass through v.
    """
    # Map allowed 1-paths of T\v into the full tournament's indexing
    # verts_del[i] = original vertex index for vertex i in T\v

    # Build allowed 1-paths and 2-paths for full tournament
    allowed_1_full = enumerate_allowed_paths(A_full, n, 1)
    idx_1_full = {path: i for i, path in enumerate(allowed_1_full)}

    # The cycle z is over allowed_1_del; translate to full indices
    z_full = np.zeros(len(allowed_1_full))
    for i, path_del in enumerate(allowed_1_del):
        if abs(cycle_vec[i]) < 1e-10:
            continue
        # Map to original vertices
        orig_path = tuple(verts_del[x] for x in path_del)
        if orig_path in idx_1_full:
            z_full[idx_1_full[orig_path]] = cycle_vec[i]

    # Now solve: find c (over allowed_2_full) such that ∂₂(c) = z_full
    # Build boundary matrix ∂₂: allowed_2_full → allowed_1_full
    bd2 = build_full_boundary_matrix(allowed_2_full, allowed_1_full)

    # Solve bd2 @ c = z_full (least-norm solution)
    if bd2.shape[1] == 0:
        return None, None

    # Use least squares
    result, residuals, rank, sv = np.linalg.lstsq(bd2, z_full, rcond=None)

    # Check if solution is valid
    reconstruction = bd2 @ result
    error = np.linalg.norm(reconstruction - z_full)
    if error > 1e-6:
        return None, None  # z is not a boundary in full T (shouldn't happen if β₁(T)=0)

    # Analyze which 2-paths pass through v
    fillers_through_v = {}
    fillers_not_through_v = {}
    for i, path in enumerate(allowed_2_full):
        if abs(result[i]) < 1e-10:
            continue
        if v in path:
            fillers_through_v[path] = result[i]
        else:
            fillers_not_through_v[path] = result[i]

    return fillers_through_v, fillers_not_through_v


# ============================================================
# Main analysis
# ============================================================

def analyze_tournament(A, n, verbose=False):
    """Full analysis of one tournament."""
    # Compute β₁(T)
    betti_T = path_betti_numbers(A, n, max_dim=2)
    beta1_T = betti_T[1] if len(betti_T) > 1 else 0

    if beta1_T != 0:
        return None  # Skip: β₁(T) ≠ 0

    # For each vertex, compute β₁(T\v)
    critical_verts = []
    deletion_info = {}

    for v in range(n):
        B, verts = delete_vertex(A, n, v)
        nm1 = n - 1
        betti_del = path_betti_numbers(B, nm1, max_dim=2)
        beta1_del = betti_del[1] if len(betti_del) > 1 else 0

        if beta1_del > 0:
            critical_verts.append(v)

            # Get detailed H₁ info for T\v
            h1_info = compute_homology_detailed(B, nm1, 1)

            deletion_info[v] = {
                'beta1': beta1_del,
                'betti': betti_del,
                'h1_info': h1_info,
                'sub_adj': B,
                'verts': verts,
            }

    sum_beta1 = sum(deletion_info[v]['beta1'] for v in critical_verts)

    result = {
        'n': n,
        'betti_T': betti_T,
        'critical_verts': critical_verts,
        'sum_beta1': sum_beta1,
        'deletion_info': deletion_info,
        'score_seq': score_sequence(A, n),
        't3': count_3cycles(A, n),
    }

    # Filling analysis: for each critical v, find filling chain
    if len(critical_verts) > 0:
        allowed_2_full = enumerate_allowed_paths(A, n, 2)

        for v in critical_verts:
            info = deletion_info[v]
            h1 = info['h1_info']

            if h1['ker_basis'] is not None and h1['betti'] > 0:
                # Get the cycle generator
                # ker_basis columns span ker(∂₁), im_basis columns span im(∂₂)
                # We want a representative of H₁ = ker/im
                # Take first kernel vector and project out image
                allowed_1_del = h1['allowed'][1]

                cycle_vec = h1['ker_basis'][:, 0]  # first kernel vector

                if h1['im_basis'] is not None:
                    # Project out image component
                    im = h1['im_basis']
                    proj = im @ np.linalg.lstsq(im, cycle_vec, rcond=None)[0]
                    cycle_vec = cycle_vec - proj

                # Normalize
                norm = np.linalg.norm(cycle_vec)
                if norm > 1e-10:
                    cycle_vec = cycle_vec / norm

                # Find filling in full T
                through_v, not_through_v = find_filling_chain(
                    A, n, info['sub_adj'], n-1, v, info['verts'],
                    cycle_vec, allowed_1_del, allowed_2_full
                )

                info['cycle_vec'] = cycle_vec
                info['allowed_1_del'] = allowed_1_del
                info['fillers_through_v'] = through_v
                info['fillers_not_through_v'] = not_through_v

    # Check 3-cycle structure of critical set
    if len(critical_verts) >= 2:
        # Check if any pair of critical vertices share a 3-cycle
        shared_3cycles = []
        for i, u in enumerate(critical_verts):
            for w in critical_verts[i+1:]:
                cycles = find_3cycles_through(A, n, [u, w])
                if cycles:
                    shared_3cycles.append((u, w, cycles))
        result['shared_3cycles'] = shared_3cycles

        if len(critical_verts) >= 3:
            # Check if all three are in a common 3-cycle
            for triple in combinations(critical_verts, 3):
                cycles = find_3cycles_through(A, n, list(triple))
                if cycles:
                    result['triple_3cycle'] = (triple, cycles)

    return result


def run_analysis(n_val, verbose_threshold=10):
    """Run full analysis for all tournaments on n vertices."""
    print(f"\n{'='*70}")
    print(f"FILLER VERTEX ANALYSIS: n = {n_val}")
    print(f"{'='*70}")

    total = 0
    beta1_zero_count = 0
    sum_dist = Counter()  # distribution of Σ β₁(T\v)
    max_sum = 0
    max_sum_examples = []
    critical_size_dist = Counter()

    # Track 3-cycle overlap stats
    triple_in_3cycle_count = 0
    pair_in_3cycle_count = 0
    total_with_3_critical = 0

    # Track filling chain stats
    filling_stats = {'through_v': 0, 'not_through_v': 0, 'total_analyzed': 0}

    # Detailed examples for verbose output
    interesting_examples = []

    for A in all_tournaments(n_val):
        total += 1
        result = analyze_tournament(A, n_val)

        if result is None:
            continue  # β₁(T) ≠ 0

        beta1_zero_count += 1
        s = result['sum_beta1']
        sum_dist[s] += 1
        critical_size_dist[len(result['critical_verts'])] += 1

        if s > max_sum:
            max_sum = s
            max_sum_examples = [result]
        elif s == max_sum and len(max_sum_examples) < 5:
            max_sum_examples.append(result)

        # Track structure of critical sets
        if len(result['critical_verts']) >= 3:
            total_with_3_critical += 1
            if 'triple_3cycle' in result:
                triple_in_3cycle_count += 1

        if len(result['critical_verts']) >= 2:
            if 'shared_3cycles' in result and result['shared_3cycles']:
                pair_in_3cycle_count += 1

        # Track filling chain stats
        for v in result['critical_verts']:
            info = result['deletion_info'][v]
            if 'fillers_through_v' in info and info['fillers_through_v'] is not None:
                filling_stats['total_analyzed'] += 1
                if info['fillers_through_v']:
                    filling_stats['through_v'] += 1
                if info['fillers_not_through_v']:
                    filling_stats['not_through_v'] += 1

        # Save interesting examples
        if s >= 2 and len(interesting_examples) < 20:
            interesting_examples.append(result)

    # Print results
    print(f"\nTotal tournaments: {total}")
    print(f"With β₁(T) = 0: {beta1_zero_count}")
    print(f"\nDistribution of Σ_v β₁(T\\v) among β₁(T)=0 tournaments:")
    for k in sorted(sum_dist.keys()):
        pct = 100 * sum_dist[k] / beta1_zero_count if beta1_zero_count > 0 else 0
        print(f"  Σ = {k}: {sum_dist[k]} ({pct:.1f}%)")

    print(f"\n|Critical set| distribution:")
    for k in sorted(critical_size_dist.keys()):
        pct = 100 * critical_size_dist[k] / beta1_zero_count if beta1_zero_count > 0 else 0
        print(f"  |crit| = {k}: {critical_size_dist[k]} ({pct:.1f}%)")

    print(f"\nMAXIMUM Σ_v β₁(T\\v) = {max_sum}")
    claim_holds = max_sum <= 3
    print(f"Claim (Σ ≤ 3): {'HOLDS' if claim_holds else 'FAILS!'}")

    print(f"\n3-Cycle structure of critical sets:")
    if total_with_3_critical > 0:
        print(f"  Tournaments with |crit| ≥ 3: {total_with_3_critical}")
        print(f"  ...with all 3 in a common 3-cycle: {triple_in_3cycle_count} " +
              f"({100*triple_in_3cycle_count/total_with_3_critical:.1f}%)")
    else:
        print(f"  No tournaments with |crit| ≥ 3")

    with_2_critical = sum(1 for k in critical_size_dist if k >= 2)
    total_with_2plus = sum(critical_size_dist[k] for k in critical_size_dist if k >= 2)
    if total_with_2plus > 0:
        print(f"  Tournaments with |crit| ≥ 2: {total_with_2plus}")
        print(f"  ...with some pair in a common 3-cycle: {pair_in_3cycle_count} " +
              f"({100*pair_in_3cycle_count/total_with_2plus:.1f}%)")

    print(f"\nFilling chain analysis:")
    print(f"  Total critical deletions analyzed: {filling_stats['total_analyzed']}")
    print(f"  Fillings with 2-paths through v: {filling_stats['through_v']}")
    print(f"  Fillings with 2-paths NOT through v: {filling_stats['not_through_v']}")

    # Detailed examples for max sum
    if max_sum_examples:
        print(f"\n{'='*70}")
        print(f"DETAILED EXAMPLES (Σ = {max_sum})")
        print(f"{'='*70}")
        for idx, ex in enumerate(max_sum_examples[:5]):
            print(f"\n--- Example {idx+1} ---")
            print(f"  Score seq: {ex['score_seq']}, t₃ = {ex['t3']}")
            print(f"  β(T) = {ex['betti_T']}")
            print(f"  Critical vertices: {ex['critical_verts']}")

            for v in ex['critical_verts']:
                info = ex['deletion_info'][v]
                print(f"  Vertex {v}: β₁(T\\{v}) = {info['beta1']}, β(T\\{v}) = {info['betti']}")
                if 'fillers_through_v' in info and info['fillers_through_v'] is not None:
                    n_through = len(info['fillers_through_v'])
                    n_other = len(info['fillers_not_through_v']) if info['fillers_not_through_v'] else 0
                    print(f"    Filling: {n_through} paths through v, {n_other} paths not through v")
                    if n_through > 0 and n_through <= 10:
                        for path, coeff in info['fillers_through_v'].items():
                            print(f"      {path}: coeff = {coeff:.4f}")

            if 'shared_3cycles' in ex:
                print(f"  Shared 3-cycles among critical pairs: {ex['shared_3cycles']}")
            if 'triple_3cycle' in ex:
                print(f"  Triple in common 3-cycle: {ex['triple_3cycle']}")

    # Additional: check if critical verts always form specific substructure
    if interesting_examples:
        print(f"\n{'='*70}")
        print(f"STRUCTURAL ANALYSIS OF CRITICAL SETS (Σ ≥ 2)")
        print(f"{'='*70}")

        # Check: in-degree / out-degree of critical vertices within {crit}
        for idx, ex in enumerate(interesting_examples[:10]):
            crit = ex['critical_verts']
            if len(crit) < 2:
                continue
            A = None
            # Reconstruct A from deletion info
            # Actually we need to pass A through... let me store it
            pass  # We'll handle this in v2 if needed

        # Check how many have critical set = exactly a 3-cycle
        triple_is_3cycle = 0
        triple_total = 0
        for ex in interesting_examples:
            crit = ex['critical_verts']
            if len(crit) == 3:
                triple_total += 1
                if 'triple_3cycle' in ex:
                    triple_is_3cycle += 1

        if triple_total > 0:
            print(f"  Among examples with |crit|=3 and Σ≥2:")
            print(f"    Total: {triple_total}")
            print(f"    Critical set forms a 3-cycle: {triple_is_3cycle} " +
                  f"({100*triple_is_3cycle/triple_total:.1f}%)")

    return sum_dist, max_sum, critical_size_dist


def run_deeper_analysis(n_val):
    """Deeper structural analysis: store adjacency for detailed inspection."""
    print(f"\n{'='*70}")
    print(f"DEEP STRUCTURAL ANALYSIS: n = {n_val}")
    print(f"{'='*70}")

    # Focus on tournaments with β₁(T)=0 and Σ β₁(T\v) = max value
    results_by_sum = defaultdict(list)

    for A in all_tournaments(n_val):
        betti_T = path_betti_numbers(A, n_val, max_dim=1)
        if betti_T[1] != 0:
            continue

        crit = []
        for v in range(n_val):
            B, _ = delete_vertex(A, n_val, v)
            betti_del = path_betti_numbers(B, n_val - 1, max_dim=1)
            if betti_del[1] > 0:
                crit.append(v)

        s = len(crit)  # since β₁(T\v) ∈ {0,1} for tournaments on small n
        if s >= 2:
            results_by_sum[s].append((A, crit))

    for s in sorted(results_by_sum.keys(), reverse=True):
        examples = results_by_sum[s]
        print(f"\nΣ = {s}: {len(examples)} tournaments")

        for A, crit in examples[:3]:
            print(f"  Critical set: {crit}")
            print(f"  Score seq: {score_sequence(A, n_val)}")
            print(f"  t₃ = {count_3cycles(A, n_val)}")

            # Print adjacency among critical vertices
            print(f"  Edges among critical vertices:")
            for u in crit:
                for w in crit:
                    if u < w:
                        direction = f"{u}→{w}" if A[u][w] else f"{w}→{u}"
                        print(f"    {direction}")

            # Check: which 3-cycles contain ≥2 critical vertices?
            for i in range(n_val):
                for j in range(i+1, n_val):
                    for k in range(j+1, n_val):
                        triple = {i,j,k}
                        n_crit_in = len(triple & set(crit))
                        if n_crit_in >= 2:
                            is_cycle = (A[i][j] and A[j][k] and A[k][i]) or \
                                       (A[j][i] and A[k][j] and A[i][k])
                            if is_cycle:
                                print(f"    3-cycle {{{i},{j},{k}}} contains {n_crit_in} critical verts")

            # Score of each critical vertex
            for v in crit:
                out_deg = sum(A[v])
                print(f"    Vertex {v}: out-degree = {out_deg}")
            print()


# ============================================================
# RUN
# ============================================================

if __name__ == '__main__':
    print("FILLER VERTEX ANALYSIS")
    print("Claim: β₁(T)=0 ⟹ Σ_v β₁(T\\v) ≤ 3")
    print("="*70)

    for n_val in [5, 6]:
        run_analysis(n_val)
        run_deeper_analysis(n_val)

    print("\n\nDONE.")
