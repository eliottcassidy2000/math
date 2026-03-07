#!/usr/bin/env python3
"""
Skeleton (GS flip graph) eigenvalue analysis at n=7.

Builds the 88x88 adjacency matrix for SC tournament isomorphism classes
connected by GS (Grinberg-Stanley) tiling flips, then computes eigenvalues.

At n=5: 8 SC classes, eigenvalues +/-(1+sqrt2), +/-1, +/-1, +/-(sqrt2-1).
At n=7: 88 SC classes forming an 88x88 skeleton.

Optimization: use refined invariant (score seq + A^2 row sums + A^3 diagonal)
to bucket tournaments, then compute canonical form only within each bucket.

kind-pasteur-2026-03-07
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
import time
import sys

# ============================================================
# Core tournament functions
# ============================================================

def tournament_from_tiling(n, tiling_bits):
    """Build tournament adjacency matrix from tiling bits.
    Edges indexed as: for i in range(n), for j in range(i+2, n).
    The (i, i+1) edges are always directed i -> i+1 (Hamiltonian path).
    """
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def refined_invariant(A, n):
    """Strong invariant: score seq + A^2 row sums + A^3 diagonal sums."""
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    A2 = [[sum(A[i][k]*A[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    a2_scores = tuple(sorted(sum(A2[i]) for i in range(n)))
    # A^3 diagonal = sum_j A2[i][j]*A[j][i] for each i
    a3_diag = tuple(sorted(sum(A2[i][j]*A[j][i] for j in range(n)) for i in range(n)))
    return (scores, a2_scores, a3_diag)

def canonical(A, n):
    """Canonical form via full permutation search."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def is_isomorphic(A, B, n):
    """Check if A and B are isomorphic (faster than full canonical)."""
    target = tuple(B[i][j] for i in range(n) for j in range(n) if i != j)
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if form == target:
            return True
    return False

def tiling_transpose_pairs(n):
    """Find pairs of tiling bit positions that are swapped by the
    vertex reversal i -> n-1-i (which maps T to its converse)."""
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs, fixed, seen = [], [], set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    """Generate all GS tilings (those invariant under transpose pairing)."""
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            t3 += 1
    return t3

def count_ham_paths(A, n):
    """Count Hamiltonian paths using bitmask DP."""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            if not (mask & (1 << last)):
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if A[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[full])

# ============================================================
# Phase 1: Build isomorphism class database
# ============================================================

def build_class_database(n):
    """Build class database for ALL tilings using refined invariant bucketing."""
    m = num_tiling_bits(n)
    total = 2**m

    print(f"Phase 1: Building class database for n={n}")
    print(f"  Tiling bits: {m}, Total tilings: {total}")
    t0 = time.time()

    # Step 1: Group all tilings by refined invariant
    ri_groups = defaultdict(list)
    for bits in range(total):
        A = tournament_from_tiling(n, bits)
        ri = refined_invariant(A, n)
        ri_groups[ri].append(bits)

    t1 = time.time()
    print(f"  Refined invariant grouping: {t1-t0:.1f}s, {len(ri_groups)} groups")

    # Step 2: Within each group, find isomorphism classes
    canon_db = {}      # canonical form -> class_idx
    class_list = []    # list of class info dicts
    bits_to_class = {} # tiling bits -> class_idx

    groups_done = 0
    for ri, bits_list in ri_groups.items():
        groups_done += 1
        if groups_done % 50 == 0:
            print(f"    Processing group {groups_done}/{len(ri_groups)}...", flush=True)

        # Within this group, cluster by isomorphism
        group_classes = []  # list of (canon, class_idx, rep_A)

        for bits in bits_list:
            A = tournament_from_tiling(n, bits)

            # Try to match against existing classes in this group
            matched = False
            for gc_canon, gc_idx, gc_rep in group_classes:
                if is_isomorphic(A, gc_rep, n):
                    bits_to_class[bits] = gc_idx
                    class_list[gc_idx]['tilings'].add(bits)
                    matched = True
                    break

            if not matched:
                # New class - compute canonical form
                c = canonical(A, n)
                if c in canon_db:
                    idx = canon_db[c]
                    bits_to_class[bits] = idx
                    class_list[idx]['tilings'].add(bits)
                    group_classes.append((c, idx, class_list[idx]['rep']))
                else:
                    idx = len(class_list)
                    canon_db[c] = idx
                    # Check SC: is converse isomorphic to original?
                    A_op = converse(A, n)
                    c_op = canonical(A_op, n)
                    sc = (c == c_op)
                    scores = score_seq(A, n)
                    class_list.append({
                        'canon': c, 'rep': A, 'tilings': {bits}, 'sc': sc,
                        'scores': scores, 'gs_tilings': set(),
                    })
                    bits_to_class[bits] = idx
                    group_classes.append((c, idx, A))

    t2 = time.time()
    num_classes = len(class_list)
    num_sc = sum(1 for c in class_list if c['sc'])
    print(f"  Canonical form computation: {t2-t1:.1f}s")
    print(f"  {num_classes} total classes, {num_sc} SC, {num_classes - num_sc} NSC")
    print(f"  Total phase 1 time: {t2-t0:.1f}s")

    return class_list, bits_to_class, canon_db

# ============================================================
# Phase 2: Build skeleton adjacency matrix
# ============================================================

def build_skeleton(n, class_list, bits_to_class):
    """Build the GS flip skeleton adjacency matrix for SC classes."""
    m = num_tiling_bits(n)
    pairs, fixed = tiling_transpose_pairs(n)
    gs_tilings = gen_gs_tilings(n, pairs, fixed)
    gs_dof = len(pairs) + len(fixed)

    print(f"\nPhase 2: Building GS flip skeleton")
    print(f"  GS DOF: {gs_dof}, GS tilings: {len(gs_tilings)}")

    # Mark GS tilings
    for bits in gs_tilings:
        idx = bits_to_class[bits]
        class_list[idx]['gs_tilings'].add(bits)

    # SC classes only
    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
    sc_to_local = {v: i for i, v in enumerate(sc_indices)}
    nsc = len(sc_indices)
    print(f"  SC classes: {nsc}")

    # Build adjacency matrix (weighted by number of GS tilings connecting them)
    K_weighted = np.zeros((nsc, nsc))
    K_binary = np.zeros((nsc, nsc))  # unweighted (0/1)

    same_class = 0
    cross_class = 0

    for bits in gs_tilings:
        c_from = bits_to_class[bits]
        flipped = flip_tiling(bits, m)
        c_to = bits_to_class[flipped]

        if c_from == c_to:
            same_class += 1
        else:
            cross_class += 1
            if c_from in sc_to_local and c_to in sc_to_local:
                i, j = sc_to_local[c_from], sc_to_local[c_to]
                K_weighted[i][j] += 1
                K_weighted[j][i] = K_weighted[i][j]  # symmetric
                K_binary[i][j] = 1
                K_binary[j][i] = 1

    print(f"  Same-class flips: {same_class}, Cross-class flips: {cross_class}")

    # Compute class properties
    for idx in sc_indices:
        if 't3' not in class_list[idx]:
            A = class_list[idx]['rep']
            class_list[idx]['t3'] = count_3cycles(A, n)
            class_list[idx]['H'] = count_ham_paths(A, n)

    return K_weighted, K_binary, sc_indices, sc_to_local

# ============================================================
# Phase 3: Spectral analysis
# ============================================================

def spectral_analysis(K, sc_indices, class_list, n):
    """Full spectral analysis of the skeleton adjacency matrix."""
    nsc = len(sc_indices)

    print(f"\n{'='*70}")
    print(f"SPECTRAL ANALYSIS: n={n}, {nsc}x{nsc} skeleton")
    print(f"{'='*70}")

    # Compute t3 and H for all SC classes
    print(f"\nComputing t3 and H for {nsc} SC classes...")
    t0 = time.time()
    for i, idx in enumerate(sc_indices):
        if 't3' not in class_list[idx]:
            A = class_list[idx]['rep']
            class_list[idx]['t3'] = count_3cycles(A, n)
        if 'H' not in class_list[idx]:
            A = class_list[idx]['rep']
            class_list[idx]['H'] = count_ham_paths(A, n)
        if (i+1) % 20 == 0:
            print(f"  {i+1}/{nsc} done...", flush=True)
    t1 = time.time()
    print(f"  Done in {t1-t0:.1f}s")

    # 1. Eigenvalues
    print(f"\n--- EIGENVALUES ---")
    evals = np.linalg.eigvalsh(K)
    evals_sorted = np.sort(evals)[::-1]  # descending

    # Round and find distinct
    evals_rounded = np.round(evals_sorted, 6)
    distinct = []
    multiplicities = []
    for ev in evals_rounded:
        found = False
        for i, d in enumerate(distinct):
            if abs(ev - d) < 1e-5:
                multiplicities[i] += 1
                found = True
                break
        if not found:
            distinct.append(ev)
            multiplicities.append(1)

    print(f"\nAll {nsc} eigenvalues (descending):")
    for i in range(0, len(evals_sorted), 10):
        chunk = evals_sorted[i:i+10]
        print(f"  [{i:3d}-{i+len(chunk)-1:3d}]: {[f'{e:.6f}' for e in chunk]}")

    print(f"\nDistinct eigenvalues (with multiplicities):")
    for d, m in zip(distinct, multiplicities):
        print(f"  {d:+12.6f}  (mult {m})")

    # 2. Bipartite check: eigenvalues in +/- pairs
    print(f"\n--- BIPARTITE SIGNATURE CHECK ---")
    pos_evals = sorted([e for e in evals_sorted if e > 1e-8], reverse=True)
    neg_evals = sorted([-e for e in evals_sorted if e < -1e-8], reverse=True)
    zero_evals = [e for e in evals_sorted if abs(e) < 1e-8]

    print(f"  Positive: {len(pos_evals)}, Negative: {len(neg_evals)}, Zero: {len(zero_evals)}")

    is_bipartite = True
    if len(pos_evals) != len(neg_evals):
        is_bipartite = False
        print(f"  NOT bipartite: unequal +/- counts")
    else:
        for p, n_val in zip(pos_evals, neg_evals):
            if abs(p - n_val) > 1e-4:
                is_bipartite = False
                print(f"  NOT bipartite: {p:.6f} vs {n_val:.6f}")
                break
        if is_bipartite:
            print(f"  BIPARTITE confirmed: all eigenvalues in +/- pairs")

    # Also check via t3 parity coloring
    t3_vals = [class_list[idx]['t3'] for idx in sc_indices]
    t3_parities = [t % 2 for t in t3_vals]

    # Check if K only connects different parities
    bipartite_by_t3 = True
    for i in range(nsc):
        for j in range(nsc):
            if K[i][j] > 0 and t3_parities[i] == t3_parities[j]:
                bipartite_by_t3 = False
                break
        if not bipartite_by_t3:
            break
    print(f"  Bipartite by t3 parity: {bipartite_by_t3}")
    if bipartite_by_t3:
        side_A = sum(1 for p in t3_parities if p == 1)
        side_B = sum(1 for p in t3_parities if p == 0)
        print(f"  Side A (odd t3): {side_A}, Side B (even t3): {side_B}")

    # 3. Algebraic number identification
    print(f"\n--- ALGEBRAIC NUMBER IDENTIFICATION ---")
    sqrt2 = np.sqrt(2)
    sqrt3 = np.sqrt(3)
    sqrt5 = np.sqrt(5)
    phi = (1 + sqrt5) / 2  # golden ratio
    silver = 1 + sqrt2      # silver ratio

    algebraic_candidates = {
        'sqrt(2)': sqrt2, '-sqrt(2)': -sqrt2,
        'sqrt(3)': sqrt3, '-sqrt(3)': -sqrt3,
        'sqrt(5)': sqrt5, '-sqrt(5)': -sqrt5,
        'phi': phi, '-phi': -phi,
        '1/phi': 1/phi, '-1/phi': -1/phi,
        'silver': silver, '-silver': -silver,
        'sqrt(2)-1': sqrt2-1, '-(sqrt(2)-1)': -(sqrt2-1),
    }

    # Also check small integer + sqrt combinations
    for a in range(-10, 11):
        for b in [0, 1, -1, 2, -2, 3, -3]:
            for r in [sqrt2, sqrt3, sqrt5]:
                val = a + b*r
                name = f"{a}+{b}*{['sqrt2','sqrt3','sqrt5'][[sqrt2,sqrt3,sqrt5].index(r)]}"
                if abs(val) < max(abs(evals_sorted)) + 1:
                    algebraic_candidates[name] = val

    # Also check a/b * sqrt(c)
    for a in range(1, 11):
        for b in range(1, 11):
            for r in [sqrt2, sqrt3, sqrt5]:
                val = a/b * r
                rname = ['sqrt2','sqrt3','sqrt5'][[sqrt2,sqrt3,sqrt5].index(r)]
                algebraic_candidates[f"{a}/{b}*{rname}"] = val
                algebraic_candidates[f"-{a}/{b}*{rname}"] = -val

    matched = {}
    for ev in evals_sorted:
        best_name = None
        best_diff = 1e-3  # threshold
        for name, val in algebraic_candidates.items():
            if abs(ev - val) < best_diff:
                best_diff = abs(ev - val)
                best_name = name
        if best_name:
            matched[round(ev, 6)] = best_name
        # Also check if it's a small rational
        for num in range(-50, 51):
            for den in range(1, 21):
                if abs(ev - num/den) < 1e-5:
                    if round(ev, 6) not in matched:
                        matched[round(ev, 6)] = f"{num}/{den}"

    print(f"  Identified algebraic eigenvalues:")
    for ev_val, name in sorted(matched.items(), reverse=True):
        print(f"    {ev_val:+12.6f} = {name}")

    # 4. Try to find minimal polynomials
    print(f"\n--- MINIMAL POLYNOMIAL SEARCH ---")
    for d, mult in zip(distinct, multiplicities):
        # Try to find integer polynomial of degree <= 4 that d satisfies
        x = d
        found = False
        for deg in range(1, 7):
            # Try p(x) = 0 with integer coefficients
            # Use numpy to check
            if deg == 1:
                # ax + b = 0 => x = -b/a
                for a in range(1, 20):
                    b = -round(a * x)
                    if abs(a*x + b) < 1e-4:
                        print(f"  {d:+12.6f} (mult {mult}): root of {a}x + {b}")
                        found = True
                        break
            elif deg == 2:
                # ax^2 + bx + c = 0
                for a in range(1, 20):
                    for b in range(-30, 31):
                        c = -round(a*x**2 + b*x)
                        if abs(a*x**2 + b*x + c) < 1e-3:
                            disc = b*b - 4*a*c
                            print(f"  {d:+12.6f} (mult {mult}): root of {a}x^2 + {b}x + {c}  [disc={disc}]")
                            found = True
                            break
                    if found:
                        break
            elif deg == 3:
                for a in range(1, 10):
                    for b in range(-15, 16):
                        for c in range(-30, 31):
                            dd = -round(a*x**3 + b*x**2 + c*x)
                            if abs(a*x**3 + b*x**2 + c*x + dd) < 1e-2:
                                print(f"  {d:+12.6f} (mult {mult}): root of {a}x^3 + {b}x^2 + {c}x + {dd}")
                                found = True
                                break
                        if found:
                            break
                    if found:
                        break
            if found:
                break
        if not found:
            print(f"  {d:+12.6f} (mult {mult}): no small polynomial found (deg <= 6)")

    # 5. Spectral gap
    print(f"\n--- SPECTRAL GAP ---")
    lambda1 = evals_sorted[0]
    lambda2 = evals_sorted[1]
    gap = lambda1 - lambda2
    print(f"  Largest eigenvalue:  {lambda1:.6f}")
    print(f"  Second largest:      {lambda2:.6f}")
    print(f"  Spectral gap:        {gap:.6f}")
    print(f"  Ratio lambda2/lambda1: {lambda2/lambda1:.6f}")

    # 6. Degree and t3 projections onto eigenspaces
    print(f"\n--- EIGENVECTOR PROJECTIONS ---")
    evals_full, evecs = np.linalg.eigh(K)
    idx_sorted = np.argsort(-evals_full)

    # Degree vector
    degree = K.sum(axis=1)
    print(f"\n  Degree vector stats: min={degree.min():.0f}, max={degree.max():.0f}, "
          f"mean={degree.mean():.1f}, std={degree.std():.1f}")

    # t3 vector
    t3_vec = np.array([class_list[idx]['t3'] for idx in sc_indices], dtype=float)
    print(f"  t3 vector stats: min={t3_vec.min():.0f}, max={t3_vec.max():.0f}, "
          f"mean={t3_vec.mean():.1f}, std={t3_vec.std():.1f}")

    # H vector
    H_vec = np.array([class_list[idx]['H'] for idx in sc_indices], dtype=float)
    print(f"  H vector stats: min={H_vec.min():.0f}, max={H_vec.max():.0f}, "
          f"mean={H_vec.mean():.1f}, std={H_vec.std():.1f}")

    # Project onto eigenspaces
    print(f"\n  Degree vector projection (top 10 eigenspaces):")
    deg_centered = degree - degree.mean()
    deg_norm = np.linalg.norm(deg_centered)
    if deg_norm > 0:
        for k in range(min(10, nsc)):
            j = idx_sorted[k]
            comp = np.dot(degree, evecs[:, j])
            frac = comp**2 / np.dot(degree, degree)
            print(f"    ev={evals_full[j]:+10.6f}: component={comp:+10.4f}, "
                  f"fraction={frac:.4f}")

    print(f"\n  t3 vector projection (top 10 eigenspaces):")
    for k in range(min(10, nsc)):
        j = idx_sorted[k]
        comp = np.dot(t3_vec, evecs[:, j])
        frac = comp**2 / np.dot(t3_vec, t3_vec)
        print(f"    ev={evals_full[j]:+10.6f}: component={comp:+10.4f}, "
              f"fraction={frac:.4f}")

    print(f"\n  H vector projection (top 10 eigenspaces):")
    for k in range(min(10, nsc)):
        j = idx_sorted[k]
        comp = np.dot(H_vec, evecs[:, j])
        frac = comp**2 / np.dot(H_vec, H_vec)
        print(f"    ev={evals_full[j]:+10.6f}: component={comp:+10.4f}, "
              f"fraction={frac:.4f}")

    # 7. Check if degree vector is an eigenvector
    print(f"\n--- DEGREE AS EIGENVECTOR CHECK ---")
    Kd = K @ degree
    # Is Kd proportional to degree?
    if degree.min() > 0:
        ratios = Kd / degree
        print(f"  K*degree / degree: min={ratios.min():.4f}, max={ratios.max():.4f}")
        if abs(ratios.max() - ratios.min()) < 0.01:
            print(f"  YES: degree is an eigenvector with eigenvalue {ratios.mean():.6f}")
        else:
            print(f"  NO: degree is not an eigenvector (ratio varies)")

    # 8. Connected components
    print(f"\n--- CONNECTIVITY ---")
    visited = set()
    components = []
    for start in range(nsc):
        if start in visited:
            continue
        comp = set()
        stack = [start]
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            visited.add(u)
            for v in range(nsc):
                if K[u][v] > 0 and v not in comp:
                    stack.append(v)
        components.append(comp)

    print(f"  Connected components: {len(components)}")
    for i, comp in enumerate(components):
        print(f"    Component {i}: {len(comp)} vertices")

    # 9. Characteristic polynomial (if small enough)
    if nsc <= 100:
        print(f"\n--- CHARACTERISTIC POLYNOMIAL ---")
        charpoly = np.round(np.poly(evals_sorted), 2)
        print(f"  Degree: {len(charpoly)-1}")
        # Just show first and last few coefficients
        if len(charpoly) > 20:
            print(f"  First 10 coeffs: {charpoly[:10].tolist()}")
            print(f"  Last 10 coeffs:  {charpoly[-10:].tolist()}")
        else:
            print(f"  Coefficients: {charpoly.tolist()}")

    return evals_sorted, evecs, evals_full

# ============================================================
# Main
# ============================================================

def main():
    n = 7
    t_start = time.time()

    print(f"SKELETON EIGENVALUE ANALYSIS at n={n}")
    print(f"{'='*70}\n")

    # Build class database
    class_list, bits_to_class, canon_db = build_class_database(n)

    # Build skeleton
    K_weighted, K_binary, sc_indices, sc_to_local = build_skeleton(n, class_list, bits_to_class)

    # Use the WEIGHTED adjacency matrix for spectral analysis
    # (K[i,j] = number of GS tilings in class i that flip to class j)
    print(f"\nUsing WEIGHTED adjacency matrix for spectral analysis")
    print(f"  (K[i,j] = # GS tilings in class i that flip to class j)")

    evals_w, evecs_w, evals_full_w = spectral_analysis(
        K_weighted, sc_indices, class_list, n)

    # Also analyze the BINARY (0/1) adjacency matrix
    print(f"\n\n{'='*70}")
    print(f"BINARY (0/1) SKELETON ANALYSIS")
    print(f"{'='*70}")

    evals_b, evecs_b, evals_full_b = spectral_analysis(
        K_binary, sc_indices, class_list, n)

    t_end = time.time()
    print(f"\n\nTotal computation time: {t_end - t_start:.1f}s")

    # Save data for later use
    print(f"\n--- SAVING DATA ---")
    np.savez('C:/Users/Eliott/Documents/GitHub/math/04-computation/skeleton_n7_data.npz',
             K_weighted=K_weighted, K_binary=K_binary,
             evals_weighted=evals_w, evals_binary=evals_b,
             sc_indices=np.array(sc_indices),
             t3_vals=np.array([class_list[idx]['t3'] for idx in sc_indices]),
             H_vals=np.array([class_list[idx]['H'] for idx in sc_indices]),
             scores=[str(class_list[idx]['scores']) for idx in sc_indices])
    print(f"  Saved to skeleton_n7_data.npz")

if __name__ == '__main__':
    main()
