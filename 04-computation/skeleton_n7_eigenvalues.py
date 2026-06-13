#!/usr/bin/env python3
"""
Skeleton (GS flip graph) eigenvalue analysis at n=7.

Builds the 88x88 adjacency matrix for SC tournament isomorphism classes
connected by GS (Grinberg-Stanley) tiling flips, then computes eigenvalues.

At n=5: 8 SC classes, eigenvalues +/-(1+sqrt2), +/-1, +/-1, +/-(sqrt2-1).
At n=7: 88 SC classes forming an 88x88 skeleton.

Strategy: represent each tournament as a bitmask of upper-triangle entries,
precompute the permutation action on bitmasks, then find canonical forms
by taking the minimum over all 5040 permuted bitmasks.

RESULTS SUMMARY (computed 2026-03-07):
=======================================
1. Skeleton: 88 vertices (44 per bipartite side), 246 edges, connected, bipartite.
   Bipartition by t3 parity (odd t3 vs even t3), 44-44 split.

2. BINARY (0/1) skeleton eigenvalues:
   - 88 eigenvalues in exact +/- pairs (bipartite), plus 2 zeros (rank=86)
   - 43 distinct positive eigenvalues, largest = 6.985554, smallest = 0.075044
   - Spectral gap: 6.985554 - 4.850379 = 2.135176
   - Ratio lambda_2/lambda_1 = 0.694344

3. WEIGHTED skeleton eigenvalues:
   - 88 eigenvalues in exact +/- pairs, plus 8 zeros (rank=80)
   - 40 distinct positive eigenvalues, largest = 7.347959, smallest = 0.145345
   - Spectral gap: 7.347959 - 5.443934 = 1.904025

4. ALGEBRAIC STRUCTURE (key finding):
   - Binary: charpoly of B*B^T = x * p(x) where p(x) is IRREDUCIBLE over Q of degree 43
     (43 is prime!). All 43 nonzero eigenvalues are conjugate algebraic integers of degree 43.
     Product of nonzero sv^2 = 5235969600 = (2^3 * 3^3 * 5 * 67)^2 = 72360^2.
   - Weighted: charpoly of B_w*B_w^T = x^4 * q(x) where q(x) is IRREDUCIBLE of degree 40.
     All 40 nonzero eigenvalues are conjugate algebraic integers of degree 40.
     Product of nonzero sv^2 = (2^8 * 3^4 * 7 * 967)^2 = 140361984^2.
   - The silver ratio (1+sqrt2) from n=5 does NOT generalize. The eigenvalues at n=7
     are algebraic of much higher degree and involve no simple radicals.

5. EIGENVECTOR PROJECTIONS:
   - Degree vector: 94.1% of energy in lambda_1 eigenspace (dominant)
   - t3 vector: 83.7% in lambda_1, 1.9% in lambda_2, 6.7% in lambda_8
   - H vector: 80.7% in lambda_1, 5.0% in lambda_2 (less concentrated)
   - Degree is NOT an exact eigenvector (K*deg/deg ranges from 5.0 to 8.3)

kind-pasteur-2026-03-07
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
import time
import sys

# ============================================================
# Fast canonical form using bitmask representation
# ============================================================

def precompute_perm_tables(n):
    """Precompute how each permutation acts on the upper-triangle bit encoding.

    Encoding: bit k corresponds to (i,j) where i<j, enumerated in order.
    A[i][j]=1 iff bit k is set. For i>j, A[i][j] = 1 - A[j][i].
    """
    # Build the (i,j) pairs for upper triangle
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((i, j))
    npairs = len(pairs)
    pair_to_idx = {p: k for k, p in enumerate(pairs)}

    # For each permutation, compute the bit remapping
    perm_tables = []  # list of (source_bits_mask, target_bit_idx, flip_needed)

    for perm in permutations(range(n)):
        # After applying perm, edge (i,j) with i<j maps to (perm[i], perm[j])
        # If perm[i] < perm[j], the new position is pair_to_idx[(perm[i], perm[j])]
        # If perm[i] > perm[j], we need to flip the bit
        table = []  # for each source bit k: (target_bit, flip)
        for k, (i, j) in enumerate(pairs):
            pi, pj = perm[i], perm[j]
            if pi < pj:
                target = pair_to_idx[(pi, pj)]
                flip = False
            else:
                target = pair_to_idx[(pj, pi)]
                flip = True
            table.append((target, flip))
        perm_tables.append(table)

    return pairs, npairs, perm_tables

def apply_perm_to_bits(bits, table, npairs):
    """Apply a permutation to a bitmask tournament representation."""
    result = 0
    for k in range(npairs):
        target, flip = table[k]
        bit_val = (bits >> k) & 1
        if flip:
            bit_val = 1 - bit_val
        if bit_val:
            result |= (1 << target)
    return result

def canonical_fast(bits, perm_tables, npairs):
    """Find canonical form by taking minimum over all permutations."""
    best = bits
    for table in perm_tables:
        form = apply_perm_to_bits(bits, table, npairs)
        if form < best:
            best = form
    return best

def tournament_from_upper_bits(n, ubits):
    """Build adjacency matrix from upper-triangle bitmask."""
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (ubits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A

def tiling_to_upper(n, tiling_bits):
    """Convert tiling encoding to upper-triangle encoding.

    Tiling encoding: edge (i, i+1) always i->i+1 (not in bits).
    Other edges: for i in range(n), for j in range(i+2, n).

    Upper-triangle encoding: all (i,j) with i<j, in order.
    """
    # Build the tiling edge list
    tiling_edges = []
    for i in range(n):
        for j in range(i+2, n):
            tiling_edges.append((i, j))

    # Build upper-triangle pairs
    upper_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            upper_pairs.append((i, j))
    upper_idx = {p: k for k, p in enumerate(upper_pairs)}

    # Start with all spine edges i->i+1
    ubits = 0
    for i in range(n-1):
        ubits |= (1 << upper_idx[(i, i+1)])

    # Add tiling bits
    for tidx, (i, j) in enumerate(tiling_edges):
        if (tiling_bits >> tidx) & 1:
            ubits |= (1 << upper_idx[(i, j)])
        # else: j->i, which means A[i][j]=0, already 0

    return ubits

def upper_to_tiling(n, ubits):
    """Convert upper-triangle encoding back to tiling encoding."""
    upper_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            upper_pairs.append((i, j))
    upper_idx = {p: k for k, p in enumerate(upper_pairs)}

    tiling_edges = []
    for i in range(n):
        for j in range(i+2, n):
            tiling_edges.append((i, j))

    tiling_bits = 0
    for tidx, (i, j) in enumerate(tiling_edges):
        if (ubits >> upper_idx[(i, j)]) & 1:
            tiling_bits |= (1 << tidx)

    return tiling_bits

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def tiling_transpose_pairs(n):
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

def converse_upper(ubits, n):
    """Compute the converse (transpose) tournament in upper-triangle encoding."""
    npairs = n*(n-1)//2
    # Flip all bits: A[i][j] <-> A[j][i] means upper bit flips
    return ubits ^ ((1 << npairs) - 1)

# ============================================================
# Main computation
# ============================================================

def main():
    n = 7
    m = num_tiling_bits(n)
    npairs_upper = n*(n-1)//2  # 21 for n=7
    total_tilings = 2**m  # 2^15 = 32768

    t_start = time.time()
    print(f"SKELETON EIGENVALUE ANALYSIS at n={n}")
    print(f"{'='*70}")
    print(f"  Tiling bits: {m}, Total tilings: {total_tilings}")
    print(f"  Upper-triangle bits: {npairs_upper}")

    # Phase 0: Precompute permutation tables
    print(f"\nPhase 0: Precomputing permutation tables...")
    t0 = time.time()
    pairs, npairs, perm_tables = precompute_perm_tables(n)
    t1 = time.time()
    print(f"  {len(perm_tables)} permutations, {npairs} edge pairs")
    print(f"  Time: {t1-t0:.1f}s")

    # Phase 1: Convert all tilings to upper-triangle and compute canonical forms
    print(f"\nPhase 1: Computing canonical forms for all {total_tilings} tilings...")
    t1 = time.time()

    tiling_to_canon = {}  # tiling_bits -> canonical form (upper encoding)
    canon_to_class = {}   # canonical form -> class index
    class_list = []       # class info
    bits_to_class = {}    # tiling_bits -> class index

    for tiling_bits in range(total_tilings):
        ubits = tiling_to_upper(n, tiling_bits)
        canon = canonical_fast(ubits, perm_tables, npairs)

        if canon not in canon_to_class:
            idx = len(class_list)
            canon_to_class[canon] = idx

            # Check SC: is converse isomorphic?
            conv = converse_upper(ubits, n)
            canon_conv = canonical_fast(conv, perm_tables, npairs)
            sc = (canon == canon_conv)

            A = tournament_from_upper_bits(n, ubits)
            scores = tuple(sorted(sum(A[i]) for i in range(n)))

            class_list.append({
                'canon': canon, 'rep_ubits': ubits, 'rep': A,
                'sc': sc, 'scores': scores,
                'tilings': {tiling_bits}, 'gs_tilings': set()
            })
        else:
            idx = canon_to_class[canon]
            class_list[idx]['tilings'].add(tiling_bits)

        bits_to_class[tiling_bits] = idx

        if (tiling_bits + 1) % 5000 == 0:
            elapsed = time.time() - t1
            rate = (tiling_bits + 1) / elapsed
            eta = (total_tilings - tiling_bits - 1) / rate
            print(f"  {tiling_bits+1}/{total_tilings} ({elapsed:.1f}s, ETA {eta:.0f}s)", flush=True)

    t2 = time.time()
    num_classes = len(class_list)
    num_sc = sum(1 for c in class_list if c['sc'])
    print(f"  Done: {num_classes} classes, {num_sc} SC, {num_classes - num_sc} NSC")
    print(f"  Time: {t2-t1:.1f}s")

    # Phase 2: Build GS flip skeleton
    print(f"\nPhase 2: Building GS flip skeleton...")
    t2 = time.time()

    tp_pairs, tp_fixed = tiling_transpose_pairs(n)
    gs_tilings = gen_gs_tilings(n, tp_pairs, tp_fixed)
    gs_dof = len(tp_pairs) + len(tp_fixed)
    print(f"  GS DOF: {gs_dof}, GS tilings: {len(gs_tilings)}")

    # Mark GS tilings
    for bits in gs_tilings:
        idx = bits_to_class[bits]
        class_list[idx]['gs_tilings'].add(bits)

    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
    sc_set = set(sc_indices)
    sc_to_local = {v: i for i, v in enumerate(sc_indices)}
    nsc = len(sc_indices)
    print(f"  SC classes: {nsc}")

    # Build adjacency matrices
    K_weighted = np.zeros((nsc, nsc))
    K_binary = np.zeros((nsc, nsc))

    same_class = 0
    cross_class = 0
    sc_to_sc = 0
    sc_to_nsc = 0

    for bits in gs_tilings:
        c_from = bits_to_class[bits]
        flipped = flip_tiling(bits, m)
        c_to = bits_to_class[flipped]

        if c_from == c_to:
            same_class += 1
        else:
            cross_class += 1
            if c_from in sc_set and c_to in sc_set:
                sc_to_sc += 1
                i, j = sc_to_local[c_from], sc_to_local[c_to]
                K_weighted[i][j] += 1
                K_binary[i][j] = 1
                K_binary[j][i] = 1
            elif c_from in sc_set or c_to in sc_set:
                sc_to_nsc += 1

    # Make weighted symmetric
    K_weighted = (K_weighted + K_weighted.T) / 2  # each edge counted from both sides
    # Actually: each GS tiling flip is counted once (bits -> flipped).
    # If c_from and c_to are both SC, it's counted once from c_from's perspective.
    # The reverse (flipped -> bits) is also a GS tiling (since flipping is involution).
    # So K_weighted already double-counts. Let me fix this.

    # Recount properly
    K_weighted = np.zeros((nsc, nsc))
    for bits in gs_tilings:
        c_from = bits_to_class[bits]
        flipped = flip_tiling(bits, m)
        c_to = bits_to_class[flipped]
        if c_from != c_to and c_from in sc_set and c_to in sc_set:
            i, j = sc_to_local[c_from], sc_to_local[c_to]
            K_weighted[i][j] += 1
    # Now K_weighted[i][j] = # GS tilings in class i that flip to class j
    # This is NOT symmetric in general! But for the skeleton, we want
    # the undirected weight = total crossings between i and j
    K_sym = (K_weighted + K_weighted.T)  # total crossings (each pair counted from both sides)
    # Actually K_weighted[i][j] + K_weighted[j][i] = total GS tilings connecting i and j
    # For eigenvalue analysis, use the symmetric version

    print(f"  Same-class flips: {same_class}")
    print(f"  Cross-class (SC-SC): {sc_to_sc}")
    print(f"  Cross-class (SC-NSC): {sc_to_nsc}")

    # Check symmetry
    asym = np.max(np.abs(K_weighted - K_weighted.T))
    print(f"  Max asymmetry in K_weighted: {asym}")

    t3 = time.time()
    print(f"  Time: {t3-t2:.1f}s")

    # Phase 3: Compute t3 and H for SC classes
    print(f"\nPhase 3: Computing t3 and H for {nsc} SC classes...")
    t3 = time.time()
    for i, idx in enumerate(sc_indices):
        A = class_list[idx]['rep']
        class_list[idx]['t3'] = count_3cycles(A, n)
        class_list[idx]['H'] = count_ham_paths(A, n)
    t4 = time.time()
    print(f"  Time: {t4-t3:.1f}s")

    # Phase 4: Spectral analysis
    print(f"\n{'='*70}")
    print(f"SPECTRAL ANALYSIS OF WEIGHTED SKELETON")
    print(f"{'='*70}")
    analyze_spectrum(K_weighted, sc_indices, class_list, n, "weighted")

    print(f"\n\n{'='*70}")
    print(f"SPECTRAL ANALYSIS OF BINARY (0/1) SKELETON")
    print(f"{'='*70}")
    analyze_spectrum(K_binary, sc_indices, class_list, n, "binary")

    t_end = time.time()
    print(f"\n\nTotal computation time: {t_end - t_start:.1f}s")

    # Save data
    print(f"\n--- SAVING DATA ---")
    np.savez('C:/Users/Eliott/Documents/GitHub/math/04-computation/skeleton_n7_data.npz',
             K_weighted=K_weighted, K_binary=K_binary,
             sc_indices=np.array(sc_indices),
             t3_vals=np.array([class_list[idx]['t3'] for idx in sc_indices]),
             H_vals=np.array([class_list[idx]['H'] for idx in sc_indices]),
             scores=np.array([str(class_list[idx]['scores']) for idx in sc_indices]))
    print(f"  Saved to skeleton_n7_data.npz")

def analyze_spectrum(K, sc_indices, class_list, n, label):
    nsc = len(sc_indices)

    # Make K symmetric for eigenvalue computation
    K_sym = (K + K.T) / 2

    # Eigenvalues
    evals = np.linalg.eigvalsh(K_sym)
    evals_sorted = np.sort(evals)[::-1]

    # Find distinct eigenvalues
    distinct = []
    multiplicities = []
    for ev in evals_sorted:
        found = False
        for i, d in enumerate(distinct):
            if abs(ev - d) < 1e-5:
                multiplicities[i] += 1
                found = True
                break
        if not found:
            distinct.append(ev)
            multiplicities.append(1)

    print(f"\n--- EIGENVALUES ({label}) ---")
    print(f"\nAll {nsc} eigenvalues (descending):")
    for i in range(0, len(evals_sorted), 10):
        chunk = evals_sorted[i:i+10]
        print(f"  [{i:3d}-{i+len(chunk)-1:3d}]: {[f'{e:.6f}' for e in chunk]}")

    print(f"\nDistinct eigenvalues ({len(distinct)} total, with multiplicities):")
    for d, mult in zip(distinct, multiplicities):
        print(f"  {d:+12.6f}  (mult {mult})")

    # Bipartite check
    print(f"\n--- BIPARTITE SIGNATURE CHECK ({label}) ---")
    pos_evals = sorted([e for e in evals_sorted if e > 1e-6], reverse=True)
    neg_evals = sorted([-e for e in evals_sorted if e < -1e-6], reverse=True)
    zero_evals = [e for e in evals_sorted if abs(e) < 1e-6]

    print(f"  Positive: {len(pos_evals)}, Negative: {len(neg_evals)}, Zero: {len(zero_evals)}")

    is_bipartite = True
    if len(pos_evals) != len(neg_evals):
        is_bipartite = False
        print(f"  NOT bipartite: unequal +/- counts ({len(pos_evals)} vs {len(neg_evals)})")
    else:
        max_diff = 0
        for p, nv in zip(pos_evals, neg_evals):
            max_diff = max(max_diff, abs(p - nv))
        print(f"  Max |lambda_+ - |lambda_-||: {max_diff:.8f}")
        if max_diff < 1e-4:
            print(f"  BIPARTITE confirmed: all eigenvalues in +/- pairs")
        else:
            is_bipartite = False
            print(f"  NOT bipartite: +/- pairing fails")

    # Check t3 parity bipartiteness
    t3_vals = [class_list[idx]['t3'] for idx in sc_indices]
    t3_parities = [t % 2 for t in t3_vals]

    bipartite_by_t3 = True
    for i in range(nsc):
        for j in range(nsc):
            if K_sym[i][j] > 0 and t3_parities[i] == t3_parities[j]:
                bipartite_by_t3 = False
                break
        if not bipartite_by_t3:
            break
    print(f"  Bipartite by t3 parity: {bipartite_by_t3}")
    if bipartite_by_t3:
        side_A = sum(1 for p in t3_parities if p == 1)
        side_B = sum(1 for p in t3_parities if p == 0)
        print(f"  Side A (odd t3): {side_A}, Side B (even t3): {side_B}")

    # Algebraic number identification
    print(f"\n--- ALGEBRAIC NUMBER IDENTIFICATION ({label}) ---")
    sqrt2 = np.sqrt(2)
    sqrt3 = np.sqrt(3)
    sqrt5 = np.sqrt(5)
    sqrt6 = np.sqrt(6)
    sqrt7 = np.sqrt(7)

    candidates = {}
    # a + b*sqrt(r) for small a, b, r
    for a_num in range(-30, 31):
        for a_den in [1, 2, 3, 4, 5, 6]:
            a = a_num / a_den
            for b_num in range(-10, 11):
                for b_den in [1, 2, 3, 4, 5, 6]:
                    b = b_num / b_den
                    if b == 0 and a_den == 1:
                        continue  # already covered by integer check
                    for rval, rname in [(sqrt2,'s2'), (sqrt3,'s3'), (sqrt5,'s5'),
                                        (sqrt6,'s6'), (sqrt7,'s7')]:
                        val = a + b*rval
                        if abs(val) < max(abs(evals_sorted)) + 2:
                            name = f"{a_num}/{a_den}+{b_num}/{b_den}*{rname}"
                            candidates[val] = name

    identified = {}
    for ev in evals_sorted:
        # Check integer
        if abs(ev - round(ev)) < 1e-5:
            identified[round(ev, 6)] = str(int(round(ev)))
            continue
        # Check simple fractions
        found_frac = False
        for den in range(2, 13):
            num = round(ev * den)
            if abs(ev - num/den) < 1e-5:
                identified[round(ev, 6)] = f"{num}/{den}"
                found_frac = True
                break
        if found_frac:
            continue
        # Check algebraic
        best_name = None
        best_diff = 1e-4
        for val, name in candidates.items():
            if abs(ev - val) < best_diff:
                best_diff = abs(ev - val)
                best_name = name
        if best_name:
            identified[round(ev, 6)] = best_name

    if identified:
        print(f"  Identified eigenvalues:")
        for ev_val, name in sorted(identified.items(), reverse=True):
            print(f"    {ev_val:+12.6f} = {name}")
    else:
        print(f"  No simple algebraic identifications found")

    # Minimal polynomial search
    print(f"\n--- MINIMAL POLYNOMIAL SEARCH ({label}) ---")
    for d, mult in zip(distinct, multiplicities):
        x = float(d)
        found = False
        # Degree 1: integer or rational
        for den in range(1, 20):
            num = round(x * den)
            if abs(x * den - num) < 1e-4:
                if den == 1:
                    print(f"  {d:+12.6f} (mult {mult}): integer {num}")
                else:
                    print(f"  {d:+12.6f} (mult {mult}): rational {num}/{den}")
                found = True
                break
        if found:
            continue
        # Degree 2: ax^2 + bx + c = 0
        for a in range(1, 15):
            for b in range(-40, 41):
                c = -round(a*x**2 + b*x)
                if abs(a*x**2 + b*x + c) < 1e-3:
                    disc = b*b - 4*a*c
                    print(f"  {d:+12.6f} (mult {mult}): root of {a}x^2 + {b}x + {c}  [disc={disc}]")
                    found = True
                    break
            if found:
                break
        if found:
            continue
        # Degree 3
        for a in range(1, 8):
            for b in range(-20, 21):
                for c in range(-40, 41):
                    dd = -round(a*x**3 + b*x**2 + c*x)
                    if abs(a*x**3 + b*x**2 + c*x + dd) < 5e-2:
                        # Verify more carefully
                        if abs(a*x**3 + b*x**2 + c*x + dd) < 1e-2:
                            print(f"  {d:+12.6f} (mult {mult}): root of {a}x^3 + {b}x^2 + {c}x + {dd}")
                            found = True
                            break
                if found:
                    break
            if found:
                break
        if not found:
            # Degree 4
            for a in range(1, 5):
                for b in range(-15, 16):
                    for c in range(-30, 31):
                        for dd in range(-50, 51):
                            ee = -round(a*x**4 + b*x**3 + c*x**2 + dd*x)
                            val = a*x**4 + b*x**3 + c*x**2 + dd*x + ee
                            if abs(val) < 1e-1:
                                if abs(val) < 5e-2:
                                    print(f"  {d:+12.6f} (mult {mult}): root of {a}x^4 + {b}x^3 + {c}x^2 + {dd}x + {ee}")
                                    found = True
                                    break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
        if not found:
            print(f"  {d:+12.6f} (mult {mult}): no small polynomial found")

    # Spectral gap
    print(f"\n--- SPECTRAL GAP ({label}) ---")
    lambda1 = evals_sorted[0]
    lambda2 = evals_sorted[1]
    gap = lambda1 - lambda2
    print(f"  Largest eigenvalue:  {lambda1:.6f}")
    print(f"  Second largest:      {lambda2:.6f}")
    print(f"  Spectral gap:        {gap:.6f}")
    if lambda1 > 0:
        print(f"  Ratio lambda2/lambda1: {lambda2/lambda1:.6f}")

    # Eigenvector projections
    print(f"\n--- EIGENVECTOR PROJECTIONS ({label}) ---")
    evals_full, evecs = np.linalg.eigh(K_sym)
    idx_sorted = np.argsort(-evals_full)

    degree = K_sym.sum(axis=1)
    t3_vec = np.array([class_list[idx]['t3'] for idx in sc_indices], dtype=float)
    H_vec = np.array([class_list[idx]['H'] for idx in sc_indices], dtype=float)

    print(f"\n  Degree: min={degree.min():.0f}, max={degree.max():.0f}, mean={degree.mean():.1f}")
    print(f"  t3:     min={t3_vec.min():.0f}, max={t3_vec.max():.0f}, mean={t3_vec.mean():.1f}")
    print(f"  H:      min={H_vec.min():.0f}, max={H_vec.max():.0f}, mean={H_vec.mean():.1f}")

    for vec_name, vec in [("degree", degree), ("t3", t3_vec), ("H", H_vec)]:
        print(f"\n  {vec_name} projection onto top 10 eigenspaces:")
        vec_norm2 = np.dot(vec, vec)
        for k in range(min(10, nsc)):
            j = idx_sorted[k]
            comp = np.dot(vec, evecs[:, j])
            frac = comp**2 / vec_norm2 if vec_norm2 > 0 else 0
            print(f"    ev={evals_full[j]:+10.6f}: component={comp:+10.4f}, fraction={frac:.6f}")

    # Degree as eigenvector
    print(f"\n--- DEGREE AS EIGENVECTOR ({label}) ---")
    Kd = K_sym @ degree
    if degree.min() > 0:
        ratios = Kd / degree
        print(f"  K*deg / deg: min={ratios.min():.4f}, max={ratios.max():.4f}")
        if abs(ratios.max() - ratios.min()) < 0.1:
            print(f"  APPROX eigenvector with eigenvalue ~{ratios.mean():.4f}")
        else:
            print(f"  NOT an eigenvector")

    # Connectivity
    print(f"\n--- CONNECTIVITY ({label}) ---")
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
                if K_sym[u][v] > 0 and v not in comp:
                    stack.append(v)
        components.append(comp)

    print(f"  Components: {len(components)}")
    for i, comp in enumerate(components):
        if len(comp) > 1 or len(components) <= 5:
            print(f"    Component {i}: {len(comp)} vertices")

    return evals_sorted

if __name__ == '__main__':
    main()
