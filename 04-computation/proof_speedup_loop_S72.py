#!/usr/bin/env python3
"""
proof_speedup_loop_S72.py — opus-2026-03-15-S72

PROOF ↔ SPEEDUP FEEDBACK LOOP EXPLORER
========================================

Systematically hunts for formulas that, if proved, enable massive speedups.
Each discovery feeds back: formula → skip computation → reach larger n → find more formulas.

TARGETS:
  A. Can Ω_k be computed from score sequence alone? (skip constraint matrix)
  B. Does β₁ satisfy deletion-contraction? (O(2^n) instead of O(n³))
  C. Is rank(∂_k|Ω_k) determined by (n, β₁, k)? (skip Gaussian elimination)
  D. Does eigenspace rank constancy hold for ALL tournaments? (n× speedup)
  E. Is dim(Ω_k) a polynomial in cycle counts? (instant Ω computation)
  F. Can β₃ be predicted from β₁ and score sequence? (skip β₃ computation)
"""

import numpy as np
from math import comb, factorial
from itertools import combinations
from collections import Counter, defaultdict
import sys
sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# Core tournament utilities (self-contained for speed)
# ============================================================

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def score_seq(A, n):
    return tuple(sorted(int(np.sum(A[i])) for i in range(n)))

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i,j] and A[j,k] and A[k,i]) or \
                   (A[i,k] and A[k,j] and A[j,i]):
                    c3 += 1
    return c3

def ham_path_count(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u] == 1:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# ============================================================
# GLMY chain complex (mod-p, fast)
# ============================================================

PRIME = 2**31 - 1

def enum_allowed(A, n, length):
    if length == 0:
        return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)
    paths = []
    for start in range(n):
        stack = [([start], 1 << start)]
        while stack:
            path, visited = stack.pop()
            if len(path) == length + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def gauss_rank_np(M, prime):
    M = M.copy() % prime
    nrows, ncols = M.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[rank])) % prime
        rank += 1
    return rank

def gauss_nullbasis_np(M, prime):
    M = M.copy() % prime
    nrows, ncols = M.shape
    pivot_cols = []
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[rank])) % prime
        pivot_cols.append(col)
        rank += 1
    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    if not free_cols:
        return rank, np.zeros((0, ncols), dtype=np.int64)
    basis = np.zeros((len(free_cols), ncols), dtype=np.int64)
    for i, fc in enumerate(free_cols):
        basis[i, fc] = 1
        for r, pc in enumerate(pivot_cols):
            if r < M.shape[0]:
                basis[i, pc] = (-M[r, fc]) % prime
    return rank, basis

def full_chain_data(A, n, max_p=None):
    """Compute full chain complex data: allowed counts, Ω dims, boundary ranks, Betti."""
    if max_p is None:
        max_p = n - 1
    max_p = min(max_p, n - 1)

    allowed = {}
    for p in range(max_p + 2):
        if p > n - 1:
            allowed[p] = []
        else:
            allowed[p] = enum_allowed(A, n, p)

    # Omega dims and bases
    omega_dims = {}
    omega_bases = {}
    for p in range(max_p + 1):
        paths = allowed.get(p, [])
        if not paths:
            omega_dims[p] = 0
            omega_bases[p] = None
            continue
        if p == 0:
            omega_dims[p] = len(paths)
            omega_bases[p] = None
            continue

        set_prev = set(allowed.get(p-1, []))
        non_allowed = {}
        na_count = 0
        for path in paths:
            for i in range(len(path)):
                face = path[:i] + path[i+1:]
                if len(set(face)) == len(face) and face not in set_prev:
                    if face not in non_allowed:
                        non_allowed[face] = na_count
                        na_count += 1

        if na_count == 0:
            omega_dims[p] = len(paths)
            omega_bases[p] = None
            continue

        P = np.zeros((na_count, len(paths)), dtype=np.int64)
        for j, path in enumerate(paths):
            for i in range(len(path)):
                face = path[:i] + path[i+1:]
                if face in non_allowed:
                    P[non_allowed[face], j] = (P[non_allowed[face], j] + (-1)**i) % PRIME

        rank_P, nbasis = gauss_nullbasis_np(P % PRIME, PRIME)
        omega_dims[p] = len(paths) - rank_P
        omega_bases[p] = nbasis if len(nbasis) > 0 else None

    # Boundary ranks
    boundary_ranks = {}
    for p in range(1, max_p + 1):
        if omega_dims.get(p, 0) == 0 or not allowed.get(p-1, []):
            boundary_ranks[p] = 0
            continue
        idx_prev = {s: i for i, s in enumerate(allowed[p-1])}
        bd = np.zeros((len(allowed[p-1]), len(allowed[p])), dtype=np.int64)
        for j, path in enumerate(allowed[p]):
            for i in range(len(path)):
                face = path[:i] + path[i+1:]
                if face in idx_prev:
                    bd[idx_prev[face], j] = (bd[idx_prev[face], j] + (-1)**i) % PRIME

        N_p = omega_bases.get(p)
        if N_p is not None and len(N_p) > 0:
            composed = bd @ N_p.T % PRIME
            boundary_ranks[p] = gauss_rank_np(composed, PRIME)
        else:
            boundary_ranks[p] = gauss_rank_np(bd % PRIME, PRIME)

    # Betti
    betti = {}
    for p in range(max_p + 1):
        ker = omega_dims.get(p, 0) - boundary_ranks.get(p, 0)
        im_next = boundary_ranks.get(p + 1, 0)
        betti[p] = max(0, ker - im_next)

    return {
        'allowed_counts': {p: len(allowed.get(p, [])) for p in range(max_p + 1)},
        'omega_dims': omega_dims,
        'boundary_ranks': boundary_ranks,
        'betti': betti,
    }

# ============================================================
# TARGET A: Ω_k from score sequence?
# ============================================================

def test_omega_from_scores():
    """Test if Ω_k is determined by score sequence alone."""
    print("\n" + "=" * 72)
    print("TARGET A: Is Ω_k determined by score sequence?")
    print("If YES → skip constraint matrix entirely → O(n²) Ω computation")
    print("=" * 72)

    for n in [4, 5, 6]:
        print(f"\n  n={n}:")
        # Group tournaments by score sequence
        score_groups = defaultdict(list)
        m = n * (n - 1) // 2
        total = 2**m

        for bits in range(total):
            A = bits_to_adj(bits, n)
            ss = score_seq(A, n)
            data = full_chain_data(A, n, max_p=min(n-1, 4))
            omega_tuple = tuple(data['omega_dims'].get(k, 0) for k in range(min(n, 5)))
            score_groups[ss].append(omega_tuple)

        # Check: same score seq → same Ω?
        determined = True
        violations = 0
        for ss, omegas in sorted(score_groups.items()):
            unique = set(omegas)
            if len(unique) > 1:
                determined = False
                violations += 1
                print(f"    VIOLATION: score_seq={ss} → {len(unique)} different Ω vectors:")
                for u in unique:
                    count = omegas.count(u)
                    print(f"      Ω={list(u)}: {count} tournaments")

        if determined:
            print(f"    ✓ Ω IS determined by score sequence at n={n}")
            # Print the mapping
            for ss, omegas in sorted(score_groups.items()):
                print(f"    {ss} → Ω={list(omegas[0])}")
        else:
            print(f"    ✗ {violations} score sequences have multiple Ω vectors")

# ============================================================
# TARGET B: β₁ via deletion-contraction
# ============================================================

def test_beta1_deletion():
    """Test if β₁ satisfies a deletion-contraction recurrence."""
    print("\n" + "=" * 72)
    print("TARGET B: Does β₁ satisfy deletion-contraction?")
    print("If YES → O(2^n) β₁ instead of O(n³) matrix rank")
    print("=" * 72)

    for n in [4, 5, 6]:
        print(f"\n  n={n}:")
        m = n * (n - 1) // 2
        total = 2**m
        successes = 0
        failures = 0

        # For each tournament, delete each vertex, check recurrence
        # β₁(T) = ? from {β₁(T\v) : v ∈ V(T)}
        beta1_data = []

        for bits in range(total):
            A = bits_to_adj(bits, n)
            data = full_chain_data(A, n, max_p=2)
            b1 = data['betti'].get(1, 0)

            # Delete each vertex
            del_betas = []
            for v in range(n):
                vertices = [i for i in range(n) if i != v]
                A_del = np.zeros((n-1, n-1), dtype=np.int8)
                for i, vi in enumerate(vertices):
                    for j, vj in enumerate(vertices):
                        A_del[i, j] = A[vi, vj]
                d = full_chain_data(A_del, n-1, max_p=2)
                del_betas.append(d['betti'].get(1, 0))

            beta1_data.append((b1, tuple(del_betas), score_seq(A, n)))

        # Check: is β₁(T) determined by {β₁(T\v)}?
        # Group by deletion signature
        sig_groups = defaultdict(list)
        for b1, del_sig, ss in beta1_data:
            sig_groups[del_sig].append(b1)

        determined = True
        for sig, b1s in sig_groups.items():
            if len(set(b1s)) > 1:
                determined = False
                break

        if determined:
            print(f"    ✓ β₁ IS determined by deletion signature at n={n}")
            # Try to find the formula
            for sig, b1s in sorted(sig_groups.items()):
                b1 = b1s[0]
                sum_del = sum(sig)
                max_del = max(sig)
                print(f"    del_betas={list(sig)} → β₁={b1}  (sum={sum_del}, max={max_del})")
        else:
            print(f"    ✗ Deletion signature does NOT determine β₁ uniquely")

        # Try simpler: β₁(T) = max(β₁(T\v)) or sum or ...
        formulas = {
            'max': lambda sigs: max(sigs),
            'sum': lambda sigs: sum(sigs),
            'sum - (n-2)': lambda sigs: sum(sigs) - (n - 2),
            'sum - n + 1': lambda sigs: sum(sigs) - n + 1,
            'max + min': lambda sigs: max(sigs) + min(sigs),
        }

        for fname, formula in formulas.items():
            match = sum(1 for b1, del_sig, _ in beta1_data
                       if formula(del_sig) == b1)
            if match == len(beta1_data):
                print(f"    ★ FORMULA: β₁(T) = {fname}(β₁(T\\v₁), ..., β₁(T\\vₙ)) — PERFECT")
            elif match > 0.95 * len(beta1_data):
                print(f"    ~ {fname}: {match}/{len(beta1_data)} match ({100*match/len(beta1_data):.1f}%)")

# ============================================================
# TARGET C: rank(∂_k|Ω_k) formula
# ============================================================

def test_rank_formula():
    """Test if rank(∂_k|Ω_k) has a closed form."""
    print("\n" + "=" * 72)
    print("TARGET C: Is rank(∂_k|Ω_k) determined by (n, k, β₁, t3)?")
    print("If YES → skip Gaussian elimination entirely")
    print("=" * 72)

    for n in [4, 5, 6]:
        print(f"\n  n={n}:")
        m = n * (n - 1) // 2
        total = 2**m

        rank_data = defaultdict(list)  # (b1, t3) -> [rank_by_degree]

        for bits in range(total):
            A = bits_to_adj(bits, n)
            b1 = None
            t3 = count_3cycles(A, n)
            data = full_chain_data(A, n, max_p=min(n-1, 5))
            b1 = data['betti'].get(1, 0)

            ranks = tuple(data['boundary_ranks'].get(k, 0) for k in range(1, min(n, 6)))
            omegas = tuple(data['omega_dims'].get(k, 0) for k in range(min(n, 6)))

            rank_data[(b1, t3)].append((ranks, omegas))

        # Check: same (β₁, t3) → same ranks?
        for key in sorted(rank_data.keys()):
            entries = rank_data[key]
            unique_ranks = set(e[0] for e in entries)
            unique_omegas = set(e[1] for e in entries)
            b1, t3 = key
            if len(unique_ranks) > 1:
                print(f"    (β₁={b1}, t3={t3}): {len(unique_ranks)} different rank vectors ({len(entries)} tournaments)")
            else:
                r = entries[0][0]
                o = entries[0][1]
                print(f"    (β₁={b1}, t3={t3}): ranks={list(r)}, Ω={list(o)} [{len(entries)} tournaments]")

# ============================================================
# TARGET D: Eigenspace rank constancy for general tournaments
# ============================================================

def test_eigenspace_constancy():
    """Test HYP-549: is boundary rank constant across eigenspaces for ALL tournaments?"""
    print("\n" + "=" * 72)
    print("TARGET D: Eigenspace rank constancy (HYP-549)")
    print("Currently proved for interval tournaments only.")
    print("If true for ALL → n× speedup for every tournament")
    print("=" * 72)

    # For non-circulant tournaments, "eigenspace" means something different.
    # But we can test: for circulant tournaments, is rank constant?
    # And for general tournaments, does the AVERAGE rank predict the total?

    # Test on Paley tournaments
    for p in [5, 7]:
        print(f"\n  Paley P_{p}:")
        QR = set()
        for k in range(1, p):
            QR.add((k * k) % p)
        A = np.zeros((p, p), dtype=np.int8)
        for i in range(p):
            for j in range(p):
                if i != j and ((j - i) % p) in QR:
                    A[i, j] = 1

        data = full_chain_data(A, p, max_p=p-1)
        print(f"    Ω dims: {[data['omega_dims'].get(k, 0) for k in range(p)]}")
        print(f"    Boundary ranks: {[data['boundary_ranks'].get(k, 0) for k in range(1, p)]}")
        print(f"    Betti: {[data['betti'].get(k, 0) for k in range(p)]}")

        # Check: is Ω_k / p always an integer? (constancy implies this)
        for k in range(p):
            dim = data['omega_dims'].get(k, 0)
            if dim > 0 and k > 0:
                # For k=0 eigenspace, dim = C(p, k+1) / p (?)
                # For k≠0 eigenspaces: dim_total = dim_k0 + (p-1)*dim_kneq0
                # If constant: dim_total = dim_k0 + (p-1)*dim_k1
                pass

        # What matters: Ω_k (total) = dim_k0 + (p-1) * dim_k1
        # If constant across k≠0, then:
        #   dim_k1 = (Ω_k_total - dim_k0) / (p-1)
        # Test: is Ω_k_total ≡ dim_k0 mod (p-1)?
        for k in range(1, p):
            dim = data['omega_dims'].get(k, 0)
            if dim > 0:
                # The k=0 part: for circulant of prime order,
                # allowed paths that are shift-invariant
                remainder = dim % p
                print(f"    Ω_{k} = {dim},  mod {p} = {remainder}")

# ============================================================
# TARGET E: Ω_k as polynomial in cycle counts
# ============================================================

def test_omega_polynomial():
    """Test if Ω_k = polynomial in (n, t3, t5, ...)."""
    print("\n" + "=" * 72)
    print("TARGET E: Is Ω_k a polynomial in cycle counts?")
    print("If YES → instant Ω computation from cycle census")
    print("=" * 72)

    for n in [5, 6]:
        print(f"\n  n={n}:")
        m = n * (n - 1) // 2
        total = 2**m

        # Collect (t3, Ω_k) for each tournament
        data_by_t3 = defaultdict(list)

        for bits in range(total):
            A = bits_to_adj(bits, n)
            t3 = count_3cycles(A, n)
            chain = full_chain_data(A, n, max_p=min(n-1, 4))
            omegas = tuple(chain['omega_dims'].get(k, 0) for k in range(min(n, 5)))
            data_by_t3[t3].append(omegas)

        # For each degree k, check if Ω_k is determined by t3
        for k in range(min(n, 5)):
            determined = True
            for t3, omegas_list in data_by_t3.items():
                vals = set(o[k] for o in omegas_list)
                if len(vals) > 1:
                    determined = False
                    break

            if determined:
                # Build the function t3 → Ω_k
                mapping = {t3: omegas_list[0][k] for t3, omegas_list in data_by_t3.items()}
                t3_vals = sorted(mapping.keys())
                omega_vals = [mapping[t] for t in t3_vals]
                print(f"    Ω_{k} IS determined by t3: {dict(zip(t3_vals, omega_vals))}")

                # Try polynomial fit
                if len(t3_vals) >= 2:
                    coeffs = np.polyfit(t3_vals, omega_vals, min(3, len(t3_vals)-1))
                    pred = np.polyval(coeffs, t3_vals)
                    residual = max(abs(p - a) for p, a in zip(pred, omega_vals))
                    if residual < 0.01:
                        # Clean up coefficients
                        clean = [round(c, 6) for c in coeffs]
                        print(f"      Polynomial fit: Ω_{k} ≈ {clean} (max residual {residual:.6f})")
            else:
                # How many t3 values have unique Ω_k?
                unique_count = sum(1 for t3, ol in data_by_t3.items()
                                  if len(set(o[k] for o in ol)) == 1)
                print(f"    Ω_{k} NOT determined by t3 alone ({unique_count}/{len(data_by_t3)} t3-values are unique)")

# ============================================================
# TARGET F: β₃ from β₁ and structure
# ============================================================

def test_beta3_prediction():
    """Test if β₃ can be predicted from β₁ and simple invariants."""
    print("\n" + "=" * 72)
    print("TARGET F: Can β₃ be predicted from β₁ + score sequence?")
    print("If YES + β₁·β₃=0 → skip half the Betti computation")
    print("=" * 72)

    for n in [5, 6, 7]:
        print(f"\n  n={n}:")
        m_edges = n * (n - 1) // 2

        if n <= 6:
            total = 2**m_edges
            sample = range(total)
            exhaustive = True
        else:
            total = 2**m_edges
            np.random.seed(42)
            sample = [np.random.randint(0, total) for _ in range(500)]
            exhaustive = False

        beta_pairs = []
        for bits in sample:
            A = bits_to_adj(bits, n)
            data = full_chain_data(A, n, max_p=min(n-1, 5))
            b1 = data['betti'].get(1, 0)
            b3 = data['betti'].get(3, 0)
            t3 = count_3cycles(A, n)
            beta_pairs.append((b1, b3, t3))

        # Check mutual exclusivity
        both_nonzero = sum(1 for b1, b3, _ in beta_pairs if b1 > 0 and b3 > 0)
        total_checked = len(beta_pairs)
        print(f"    β₁·β₃ = 0: {'CONFIRMED' if both_nonzero == 0 else f'VIOLATED {both_nonzero} times'} "
              f"({'exhaustive' if exhaustive else 'sampled'}, {total_checked} tournaments)")

        # When β₃ > 0, what determines it?
        b3_nonzero = [(b1, b3, t3) for b1, b3, t3 in beta_pairs if b3 > 0]
        if b3_nonzero:
            print(f"    β₃ > 0: {len(b3_nonzero)} tournaments")
            t3_vals = set(t3 for _, _, t3 in b3_nonzero)
            print(f"    t3 values when β₃>0: {sorted(t3_vals)}")

            # Is β₃ determined by t3?
            t3_to_b3 = defaultdict(set)
            for _, b3, t3 in b3_nonzero:
                t3_to_b3[t3].add(b3)
            all_unique = all(len(v) == 1 for v in t3_to_b3.values())
            if all_unique:
                print(f"    β₃ IS determined by t3 when β₃>0")
            else:
                print(f"    β₃ NOT determined by t3 alone")

        # β₁ distribution
        b1_dist = Counter(b1 for b1, _, _ in beta_pairs)
        b3_dist = Counter(b3 for _, b3, _ in beta_pairs)
        print(f"    β₁ distribution: {dict(sorted(b1_dist.items()))}")
        print(f"    β₃ distribution: {dict(sorted(b3_dist.items()))}")

# ============================================================
# TARGET G: H from deletion values (new DC for H)
# ============================================================

def test_H_deletion_formula():
    """Test if H(T) has a simple formula from H(T\\v) values."""
    print("\n" + "=" * 72)
    print("TARGET G: Is H(T) = f(H(T\\v₁), ..., H(T\\vₙ))?")
    print("Known: H(T) = Σ H(T\\v_last) over vertices with in-degree pattern")
    print("If exact formula exists → O(n · 2^(n-1)) H computation")
    print("=" * 72)

    for n in [4, 5]:
        print(f"\n  n={n}:")
        m = n * (n - 1) // 2
        total = 2**m

        # Check: H(T) = sum of H(T\v) / something?
        relations = []
        for bits in range(total):
            A = bits_to_adj(bits, n)
            H = ham_path_count(A, n)

            del_H = []
            for v in range(n):
                vertices = [i for i in range(n) if i != v]
                A_del = np.zeros((n-1, n-1), dtype=np.int8)
                for i, vi in enumerate(vertices):
                    for j, vj in enumerate(vertices):
                        A_del[i, j] = A[vi, vj]
                del_H.append(ham_path_count(A_del, n-1))

            relations.append((H, tuple(del_H)))

        # Test: H = sum(del_H)?
        sum_match = sum(1 for H, dH in relations if H == sum(dH))
        print(f"    H = Σ H(T\\v): {sum_match}/{total} match")

        # Test: H = sum(del_H) / (n-1)?
        div_match = sum(1 for H, dH in relations if H * (n-1) == sum(dH) * n)
        print(f"    H·(n-1) = n·Σ H(T\\v): {div_match}/{total} match")

        # What IS the actual ratio?
        ratios = [(H, sum(dH), sum(dH)/H if H > 0 else 0) for H, dH in relations]
        ratio_vals = set(r for _, _, r in ratios if r > 0)
        if len(ratio_vals) == 1:
            print(f"    ★ UNIVERSAL RATIO: Σ H(T\\v) / H(T) = {ratio_vals.pop()}")
        else:
            print(f"    Ratio Σ H(T\\v)/H(T) varies: {sorted(ratio_vals)[:5]}...")
            # Check if ratio is always n
            n_match = sum(1 for H, sH, _ in ratios if sH == n * H)
            print(f"    Σ H(T\\v) = n · H(T): {n_match}/{total} match")

# ============================================================
# MAIN
# ============================================================

print("=" * 72)
print("PROOF ↔ SPEEDUP FEEDBACK LOOP EXPLORER")
print("opus-2026-03-15-S72")
print("=" * 72)
print("""
Strategy: Find formulas that skip expensive computation.
Each proved formula → larger n reachable → more patterns → more proofs.

The loop: PROVE formula → SKIP matrix op → COMPUTE at 10× larger n
           → SEE new pattern → PROVE next formula → ...
""")

test_omega_from_scores()  # Target A
test_beta1_deletion()     # Target B
test_rank_formula()       # Target C
test_omega_polynomial()   # Target E
test_beta3_prediction()   # Target F
test_H_deletion_formula() # Target G

# Target D done inline above (eigenspace)

print("\n" + "=" * 72)
print("SYNTHESIS: WHICH TARGETS HIT?")
print("=" * 72)
print("""
Results will show which formulas hold and which break.
Any CONFIRMED formula is a proof target that unlocks speedup.
Any REFUTED formula tells us which invariants are NOT sufficient.

The feedback loop works when:
  FORMULA CONFIRMED → prove it → skip computation → reach larger n
  FORMULA REFUTED → tells us WHAT extra info is needed → refine → retry
""")

print("\nDONE — proof_speedup_loop_S72.py complete")
