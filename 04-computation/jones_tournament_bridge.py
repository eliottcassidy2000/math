#!/usr/bin/env python3
"""
jones_tournament_bridge.py
==========================
Investigates topological polynomial invariants for tournaments, inspired by
knot theory (Jones polynomial / Kauffman bracket).

For n=3,4,5,6,7:
1. Compute GLMY path homology Betti numbers beta_0..beta_{n-1}
2. Euler characteristic chi = sum (-1)^k beta_k
3. H(T) = number of Hamiltonian paths
4. "Jones-like polynomial" J_T(q) = sum_k beta_k q^k, evaluated at q=-1 (=chi), q=1 (=total Betti)
5. Relationship between chi and H(T)
6. "Kauffman state sum" K(A) for small n=3,4,5: generating function of H over all tournaments

Uses GLMY path homology code from glmy_walsh_bridge.py.
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import sys
import random
import time

# ============================================================
# CORE: Path Homology (from glmy_walsh_bridge.py)
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
        beta_p = ker_dim - im_dim
        betti.append(max(0, beta_p))
    return betti

# ============================================================
# Tournament utilities
# ============================================================

def all_tournaments(n):
    """Enumerate all 2^C(n,2) tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A, mask

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def ham_path_count(A, n):
    """Count Hamiltonian paths via DP (Bellman-Held-Karp)."""
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
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def score_sequence(A, n):
    scores = sorted([sum(A[i]) for i in range(n)])
    return tuple(scores)

# ============================================================
# Analysis functions
# ============================================================

def euler_char(betti):
    return sum((-1)**k * b for k, b in enumerate(betti))

def jones_poly_eval(betti, q):
    """J_T(q) = sum_k beta_k * q^k, evaluated at given q."""
    return sum(b * q**k for k, b in enumerate(betti))

def jones_poly_str(betti):
    """String representation of J_T(q) = sum beta_k q^k."""
    terms = []
    for k, b in enumerate(betti):
        if b == 0:
            continue
        if k == 0:
            terms.append(str(b))
        elif k == 1:
            terms.append(f"{b}q" if b != 1 else "q")
        else:
            terms.append(f"{b}q^{k}" if b != 1 else f"q^{k}")
    return " + ".join(terms) if terms else "0"

# ============================================================
# PART 1-3: Betti numbers, Euler char, Jones-like polynomial
# ============================================================

def analyze_tournaments(n, sample_size=None):
    """Analyze all (or sampled) tournaments at given n."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    total = 1 << m

    if sample_size is not None and sample_size < total:
        print(f"\n  Sampling {sample_size} tournaments (out of {total})")
        is_sampled = True
    else:
        sample_size = total
        print(f"\n  Enumerating all {total} tournaments")
        is_sampled = False

    records = []
    t0 = time.time()

    if is_sampled:
        for idx in range(sample_size):
            if (idx + 1) % 500 == 0:
                elapsed = time.time() - t0
                rate = (idx + 1) / elapsed
                print(f"    ... {idx+1}/{sample_size} ({rate:.1f}/s)", file=sys.stderr)
            A = random_tournament(n)
            betti = path_betti_numbers(A, n)
            H = ham_path_count(A, n)
            chi = euler_char(betti)
            total_betti = sum(betti)
            score = score_sequence(A, n)
            records.append({
                'betti': tuple(betti),
                'H': H,
                'chi': chi,
                'total_betti': total_betti,
                'score': score,
            })
    else:
        for idx, (A, mask) in enumerate(all_tournaments(n)):
            if (idx + 1) % 1000 == 0:
                elapsed = time.time() - t0
                rate = (idx + 1) / elapsed
                print(f"    ... {idx+1}/{total} ({rate:.1f}/s)", file=sys.stderr)
            betti = path_betti_numbers(A, n)
            H = ham_path_count(A, n)
            chi = euler_char(betti)
            total_betti = sum(betti)
            score = score_sequence(A, n)
            records.append({
                'betti': tuple(betti),
                'H': H,
                'chi': chi,
                'total_betti': total_betti,
                'score': score,
            })

    elapsed = time.time() - t0
    print(f"  Completed in {elapsed:.1f}s")
    return records

# ============================================================
# PART 4: Kauffman bracket / state sum
# ============================================================

def kauffman_state_sum(n):
    """
    Compute K(A) = sum_{T'} A^{H(T')} over ALL tournaments T' on n vertices.
    This is the generating function for H values.
    Returns dict: {H_value: count}.
    """
    print(f"\n  Computing Kauffman state sum for n={n}...")
    H_dist = Counter()
    total = 1 << (n*(n-1)//2)
    for A, mask in all_tournaments(n):
        H = ham_path_count(A, n)
        H_dist[H] += 1
    return H_dist

def format_polynomial(coeffs_dict, var='A'):
    """Format a polynomial from {exponent: coeff} dict."""
    if not coeffs_dict:
        return "0"
    terms = []
    for exp in sorted(coeffs_dict.keys()):
        c = coeffs_dict[exp]
        if c == 0:
            continue
        if exp == 0:
            terms.append(str(c))
        elif exp == 1:
            terms.append(f"{c}{var}" if c != 1 else var)
        else:
            terms.append(f"{c}{var}^{exp}" if c != 1 else f"{var}^{exp}")
    return " + ".join(terms) if terms else "0"

# ============================================================
# Main
# ============================================================

def main():
    print("=" * 72)
    print("JONES-TOURNAMENT BRIDGE: Topological Polynomial Invariants")
    print("=" * 72)

    all_results = {}

    # =============================================
    # PART 1-3: Betti, chi, Jones-like for n=3..7
    # =============================================
    for n in [3, 4, 5, 6, 7]:
        print(f"\n{'='*72}")
        print(f"  n = {n}")
        print(f"{'='*72}")

        sample = 5000 if n == 7 else None
        records = analyze_tournaments(n, sample_size=sample)
        all_results[n] = records

        # --- Chi distribution ---
        chi_counts = Counter(r['chi'] for r in records)
        print(f"\n  EULER CHARACTERISTIC (chi = sum (-1)^k beta_k):")
        print(f"    Distribution: {dict(sorted(chi_counts.items()))}")
        if len(chi_counts) == 1:
            print(f"    *** chi is CONSTANT = {list(chi_counts.keys())[0]} for all tournaments ***")
        else:
            chis = [r['chi'] for r in records]
            print(f"    Range: [{min(chis)}, {max(chis)}]")
            print(f"    Mean: {np.mean(chis):.4f}")

        # --- Betti distribution ---
        betti_counts = Counter(r['betti'] for r in records)
        print(f"\n  BETTI NUMBER DISTRIBUTION (top 15 most common):")
        for betti_tuple, cnt in betti_counts.most_common(15):
            pct = 100 * cnt / len(records)
            print(f"    beta={str(list(betti_tuple)):30s}  count={cnt:6d}  ({pct:5.1f}%)")
        if len(betti_counts) > 15:
            print(f"    ... ({len(betti_counts)} distinct Betti vectors total)")

        # --- H distribution ---
        H_counts = Counter(r['H'] for r in records)
        print(f"\n  H(T) DISTRIBUTION (distinct values):")
        H_vals = sorted(H_counts.keys())
        print(f"    Values: {H_vals}")
        print(f"    Range: [{min(H_vals)}, {max(H_vals)}]")

        # --- Total Betti distribution ---
        tb_counts = Counter(r['total_betti'] for r in records)
        print(f"\n  TOTAL BETTI (J_T(1) = sum beta_k):")
        print(f"    Distribution: {dict(sorted(tb_counts.items()))}")

        # --- Jones polynomial J_T(q) examples ---
        print(f"\n  JONES-LIKE POLYNOMIAL J_T(q) = sum beta_k q^k (examples):")
        seen_betti = set()
        examples_shown = 0
        for r in records:
            b = r['betti']
            if b not in seen_betti:
                seen_betti.add(b)
                jp = jones_poly_str(list(b))
                print(f"    beta={str(list(b)):30s} => J_T(q) = {jp}")
                print(f"       J_T(-1) = chi = {r['chi']},  J_T(1) = {r['total_betti']},  H = {r['H']}")
                examples_shown += 1
                if examples_shown >= 20:
                    break

        # --- Chi vs H relationship ---
        print(f"\n  CHI vs H(T) RELATIONSHIP:")
        chi_H_pairs = defaultdict(list)
        for r in records:
            chi_H_pairs[r['chi']].append(r['H'])
        for chi_val in sorted(chi_H_pairs.keys()):
            Hs = chi_H_pairs[chi_val]
            print(f"    chi={chi_val:3d}: H values = {sorted(set(Hs))}, count={len(Hs)}")

        # --- Check if chi == 1 always ---
        all_chi_one = all(r['chi'] == 1 for r in records)
        print(f"\n  chi == 1 for all tournaments? {all_chi_one}")

        # --- Check chi mod relationships ---
        print(f"\n  CHI mod PATTERNS:")
        for mod in [2, 3, 4]:
            chi_mod = Counter(r['chi'] % mod for r in records)
            print(f"    chi mod {mod}: {dict(sorted(chi_mod.items()))}")

        # --- H mod chi ---
        if not all_chi_one:
            print(f"\n  H(T) mod chi(T) (when chi != 0):")
            for r in records[:20]:
                if r['chi'] != 0:
                    print(f"    H={r['H']}, chi={r['chi']}, H mod chi = {r['H'] % r['chi']}")

        # --- Correlation ---
        if len(set(r['chi'] for r in records)) > 1 and len(set(r['H'] for r in records)) > 1:
            chis = np.array([r['chi'] for r in records], dtype=float)
            Hs = np.array([r['H'] for r in records], dtype=float)
            corr = np.corrcoef(chis, Hs)[0, 1]
            print(f"\n  Pearson correlation(chi, H) = {corr:.6f}")

            tbs = np.array([r['total_betti'] for r in records], dtype=float)
            corr2 = np.corrcoef(tbs, Hs)[0, 1]
            print(f"  Pearson correlation(total_betti, H) = {corr2:.6f}")

        # --- Beta_1 vs H ---
        if n >= 3:
            b1_H = defaultdict(list)
            for r in records:
                b1_H[r['betti'][1] if len(r['betti']) > 1 else 0].append(r['H'])
            print(f"\n  BETA_1 vs H:")
            for b1 in sorted(b1_H.keys()):
                Hs = b1_H[b1]
                print(f"    beta_1={b1:3d}: mean H = {np.mean(Hs):8.2f}, "
                      f"range [{min(Hs)}, {max(Hs)}], count={len(Hs)}")

    # =============================================
    # PART 4: Kauffman state sum for n=3,4,5
    # =============================================
    print(f"\n{'='*72}")
    print("PART 4: KAUFFMAN STATE SUM K(A) = sum_{{T'}} A^{{H(T')}}")
    print("=" * 72)

    for n in [3, 4, 5]:
        print(f"\n  --- n = {n} ---")
        H_dist = kauffman_state_sum(n)
        total_tours = sum(H_dist.values())
        print(f"  Total tournaments: {total_tours}")
        print(f"  H distribution: {dict(sorted(H_dist.items()))}")
        poly_str = format_polynomial(H_dist, 'A')
        print(f"  K(A) = {poly_str}")

        # Check for factorization patterns
        # K(A) at A=1: number of tournaments
        K_at_1 = sum(H_dist.values())
        K_at_neg1 = sum(cnt * ((-1)**H) for H, cnt in H_dist.items())
        K_at_2 = sum(cnt * (2**H) for H, cnt in H_dist.items())
        print(f"  K(1) = {K_at_1} (= 2^C({n},2) = {1 << (n*(n-1)//2)})")
        print(f"  K(-1) = {K_at_neg1}")
        print(f"  K(2) = {K_at_2}")

        # Check if all H values are odd
        all_odd_H = all(H % 2 == 1 for H in H_dist.keys())
        print(f"  All H values odd? {all_odd_H}")

        # Parity of K
        even_count = sum(cnt for H, cnt in H_dist.items() if H % 2 == 0)
        odd_count = sum(cnt for H, cnt in H_dist.items() if H % 2 == 1)
        print(f"  Tournaments with even H: {even_count}, odd H: {odd_count}")

        # Try to express K(A) in factored form
        # K(A) = product over some simple factors?
        # Check if K(A) = C * A^{min_H} * P(A) where P is nice
        min_H = min(H_dist.keys())
        shifted = {H - min_H: cnt for H, cnt in H_dist.items()}
        print(f"  K(A) / A^{min_H} = {format_polynomial(shifted, 'A')}")

        # Check symmetry: is K(A) = A^S * K(1/A) for some S?
        max_H = max(H_dist.keys())
        is_palindrome = True
        for H in H_dist:
            mirror = max_H + min_H - H
            if H_dist.get(mirror, 0) != H_dist[H]:
                is_palindrome = False
                break
        print(f"  Palindromic around (min+max)/2? {is_palindrome}")

        # GCD of coefficients
        from math import gcd
        from functools import reduce
        all_coeffs = list(H_dist.values())
        g = reduce(gcd, all_coeffs)
        print(f"  GCD of all coefficients: {g}")
        if g > 1:
            reduced = {H: cnt // g for H, cnt in H_dist.items()}
            print(f"  K(A) / {g} = {format_polynomial(reduced, 'A')}")

    # =============================================
    # PART 5: Deep chi analysis across all n
    # =============================================
    print(f"\n{'='*72}")
    print("PART 5: DEEP CHI ANALYSIS")
    print("=" * 72)

    print(f"\n  Summary table:")
    print(f"  {'n':>3s} {'#tours':>8s} {'chi values':>30s} {'chi constant?':>15s}")
    print(f"  {'---':>3s} {'--------':>8s} {'------------------------------':>30s} {'---------------':>15s}")
    for n in [3, 4, 5, 6, 7]:
        records = all_results[n]
        chi_set = sorted(set(r['chi'] for r in records))
        is_const = "YES" if len(chi_set) == 1 else "NO"
        chi_str = str(chi_set) if len(chi_set) <= 10 else f"[{min(chi_set)}..{max(chi_set)}], {len(chi_set)} values"
        num = len(records)
        print(f"  {n:3d} {num:8d} {chi_str:>30s} {is_const:>15s}")

    # =============================================
    # PART 6: Jones polynomial type distribution
    # =============================================
    print(f"\n{'='*72}")
    print("PART 6: JONES POLYNOMIAL TYPE DISTRIBUTION")
    print("=" * 72)

    for n in [3, 4, 5, 6]:
        records = all_results[n]
        jp_types = Counter()
        for r in records:
            jp_types[r['betti']] += 1
        print(f"\n  n={n}: {len(jp_types)} distinct J_T(q) types")
        for betti, cnt in jp_types.most_common():
            jp = jones_poly_str(list(betti))
            chi = euler_char(list(betti))
            tb = sum(betti)
            pct = 100 * cnt / len(records)
            print(f"    J_T(q) = {jp:35s}  chi={chi:3d}  sum={tb:3d}  count={cnt:5d} ({pct:5.1f}%)")

    # =============================================
    # PART 7: Score sequence vs chi
    # =============================================
    print(f"\n{'='*72}")
    print("PART 7: SCORE SEQUENCE vs CHI")
    print("=" * 72)

    for n in [3, 4, 5, 6]:
        records = all_results[n]
        score_chi = defaultdict(set)
        score_H = defaultdict(set)
        for r in records:
            score_chi[r['score']].add(r['chi'])
            score_H[r['score']].add(r['H'])
        print(f"\n  n={n}:")
        for score in sorted(score_chi.keys()):
            chis = sorted(score_chi[score])
            Hs = sorted(score_H[score])
            print(f"    score={list(score)}: chi={chis}, H={Hs}")

    # =============================================
    # PART 8: Check beta_0 constancy
    # =============================================
    print(f"\n{'='*72}")
    print("PART 8: BETA_0 (connected components)")
    print("=" * 72)

    for n in [3, 4, 5, 6, 7]:
        records = all_results[n]
        b0_counts = Counter(r['betti'][0] for r in records)
        print(f"  n={n}: beta_0 distribution = {dict(sorted(b0_counts.items()))}")

    # =============================================
    # PART 9: Alternating Betti sum formulas
    # =============================================
    print(f"\n{'='*72}")
    print("PART 9: SEARCHING FOR chi-H FORMULAS")
    print("=" * 72)

    for n in [3, 4, 5, 6]:
        records = all_results[n]
        # Check: is chi related to H mod 2? mod 4?
        # Check: is chi = f(H) for some function?
        H_to_chis = defaultdict(set)
        for r in records:
            H_to_chis[r['H']].add(r['chi'])

        # Is chi a function of H?
        is_function = all(len(v) == 1 for v in H_to_chis.values())
        print(f"\n  n={n}: chi is a FUNCTION of H? {is_function}")
        if is_function:
            for H in sorted(H_to_chis.keys()):
                chi_val = list(H_to_chis[H])[0]
                print(f"    H={H:3d} => chi={chi_val}")
        else:
            for H in sorted(H_to_chis.keys()):
                if len(H_to_chis[H]) > 1:
                    print(f"    H={H:3d} => chi in {sorted(H_to_chis[H])}")

        # Check: is H congruent to chi mod 2?
        all_match_mod2 = all(r['H'] % 2 == r['chi'] % 2 for r in records)
        print(f"  n={n}: H == chi (mod 2)? {all_match_mod2}")

        # Check: is (H - chi) divisible by something?
        diffs = [r['H'] - r['chi'] for r in records]
        if diffs:
            g = reduce(gcd, [abs(d) for d in diffs if d != 0]) if any(d != 0 for d in diffs) else 0
            print(f"  n={n}: gcd(H - chi) = {g}")

    print(f"\n{'='*72}")
    print("DONE")
    print("=" * 72)

if __name__ == '__main__':
    main()
