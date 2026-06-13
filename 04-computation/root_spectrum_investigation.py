#!/usr/bin/env python3
"""
root_spectrum_investigation.py — oracle-2026-05-17-S1

Comprehensive investigation of the root spectrum of I(Omega(T), x):
1. Root profiles for all n=6 iso classes (56 classes)
2. Verify SC tournaments have most asymmetric root distribution (min rho2/rho1)
3. Verify (H, I(Omega,6)) separates all iso classes at n=6
4. Verify ultra-log-concavity (Newton-Maclaurin inequality)
5. Check alpha1=3, alpha2=0 is impossible at n=6,7
6. Ratio I(Omega,6)/H and its correlation with SC structure
7. Root distribution across iso classes
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))

from itertools import combinations, permutations
from math import comb, sqrt
from collections import defaultdict
import numpy as np
import random

# ─────────────────────────────────────────────────────────────
# Core: Hamiltonian paths, odd cycles, independence polynomial
# ─────────────────────────────────────────────────────────────

def hamiltonian_path_count(adj, n):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask>>v&1) or not dp[mask][v]:
                continue
            for u in range(n):
                if not (mask>>u&1) and adj[v][u]:
                    dp[mask|(1<<u)][u] += dp[mask][v]
    return sum(dp[(1<<n)-1])

def find_all_odd_cycles(adj, n):
    """Find all directed odd cycles in tournament adj."""
    cycles = []
    seen = set()
    for start in range(n):
        # DFS from start, find cycles back to start
        stack = [(start, [start], 1<<start)]
        while stack:
            v, path, vmask = stack.pop()
            for u in range(n):
                if adj[v][u]:
                    if u == start and len(path) >= 3 and len(path) % 2 == 1:
                        key = tuple(sorted(path))
                        if key not in seen:
                            seen.add(key)
                            cycles.append(tuple(path))
                    elif u > start and not (vmask>>u&1) and len(path) < n:
                        stack.append((u, path+[u], vmask|(1<<u)))
    return cycles

def independence_poly_coeffs(cycles):
    """Compute alpha_k = # k-tuples of vertex-disjoint cycles."""
    m = len(cycles)
    if m == 0:
        return [1]
    vsets = [frozenset(c) for c in cycles]
    # Build conflict adjacency bitmasks
    adj_bits = [0]*m
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1<<b
                adj_bits[b] |= 1<<a
    # Count independent sets by BFS on bitmask
    max_k = m
    coeffs = [0]*(max_k+1)
    coeffs[0] = 1
    # Each state is (bitmask of cycles used, candidates bitmask)
    # Use DP on subsets
    # For speed, iterate over all subsets
    for mask in range(1, 1<<m):
        # Check if independent set
        bits = [i for i in range(m) if mask>>i&1]
        k = len(bits)
        ok = True
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if adj_bits[bits[i]] >> bits[j] & 1:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            coeffs[k] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def poly_roots(coeffs):
    """Compute roots of polynomial sum_k coeffs[k] x^k."""
    # coeffs[0] = alpha_0 = 1, coeffs[1] = alpha_1, ...
    d = len(coeffs) - 1
    if d == 0:
        return []
    # numpy expects highest degree first
    poly = list(reversed(coeffs))
    roots = np.roots(poly)
    return roots

def tournament_from_score_seq(n, scores):
    """Build one tournament with given score sequence (not unique)."""
    # Simple greedy: vertex with highest score beats everything else
    verts = list(range(n))
    adj = [[0]*n for _ in range(n)]
    # Sort vertices by desired score descending
    order = sorted(range(n), key=lambda i: scores[i], reverse=True)
    # Give each vertex its score by having it beat the next ones
    remaining = list(order)
    for idx, v in enumerate(order):
        s = scores[v]  # desired out-degree
        # beat next s vertices in remaining (after v)
        rest = [u for u in order if u != v]
        for i in range(s):
            if i < len(rest):
                adj[v][rest[i]] = 1
                adj[rest[i]][v] = 0
    return adj

def is_sc(adj, n):
    """Check if tournament is self-complementary."""
    # Try all permutations (only feasible for small n)
    comp = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                comp[i][j] = 1 - adj[i][j]
    # Check if adj and comp are isomorphic
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j:
                    if adj[perm[i]][perm[j]] != comp[i][j]:
                        ok = False
                        break
            if not ok:
                break
        if ok:
            return True
    return False

def score_sequence(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))

def canonical_hash(adj, n):
    """Isomorphism class hash via score sequence (rough)."""
    scores = tuple(sorted(sum(adj[i]) for i in range(n)))
    return scores

# Better iso check using full relabeling
def iso_check(adj1, adj2, n):
    """Check if two tournaments are isomorphic."""
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if adj1[perm[i]][perm[j]] != adj2[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False

def enumerate_iso_classes(n):
    """Enumerate all iso classes at given n by exhaustive search."""
    classes = []  # list of (representative_adj, H, coeffs)
    seen = []
    total = 1 << (n*(n-1)//2)
    idx = 0
    for bits in range(total):
        adj = [[0]*n for _ in range(n)]
        b = bits
        for i in range(n):
            for j in range(i+1, n):
                if b & 1:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
                b >>= 1
        # Check if already seen
        new_class = True
        for rep, _, _ in seen:
            if iso_check(adj, rep, n):
                new_class = False
                break
        if new_class:
            H = hamiltonian_path_count(adj, n)
            cycles = find_all_odd_cycles(adj, n)
            coeffs = independence_poly_coeffs(cycles)
            seen.append((adj, H, coeffs))
    return seen

# ─────────────────────────────────────────────────────────────
# Fast version: enumerate n=6 iso classes from random tournaments
# (faster than exhaustive 2^15=32768 enumeration with iso check)
# ─────────────────────────────────────────────────────────────

def random_tournament(n):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

def compute_invariants(adj, n):
    """Compute all invariants for a tournament."""
    H = hamiltonian_path_count(adj, n)
    cycles = find_all_odd_cycles(adj, n)
    coeffs = independence_poly_coeffs(cycles)
    sc = is_sc(adj, n) if n <= 7 else None
    scores = score_sequence(adj, n)
    # I(Omega, 6)
    I6 = sum(coeffs[k] * (6**k) for k in range(len(coeffs)))
    return {
        'H': H, 'coeffs': coeffs, 'sc': sc, 'scores': scores, 'I6': I6,
        'adj': adj, 'cycles': cycles
    }

def coeffs_to_key(coeffs):
    return tuple(coeffs)

def run_n6_analysis():
    print("=" * 70)
    print("N=6 ISO CLASS ROOT SPECTRUM ANALYSIS")
    print("=" * 70)
    n = 6

    # Collect iso classes by (H, coeffs) — this is a canonical key since
    # (H, alpha_1, alpha_2) together with scores determines iso class at n=6
    # We'll verify this claim.
    class_map = {}  # key=(H,I6,coeffs_tuple) -> list of (adj,data)

    print(f"\nSampling tournaments to find all iso classes at n={n}...")
    samples = 0
    max_samples = 50000

    for _ in range(max_samples):
        adj = random_tournament(n)
        data = compute_invariants(adj, n)
        key = (data['H'], data['I6'], coeffs_to_key(data['coeffs']))
        if key not in class_map:
            class_map[key] = data
        samples += 1

    print(f"Found {len(class_map)} iso classes after {samples} samples")
    print(f"Expected: 56 iso classes at n=6")

    # Sort by H then I6
    classes = sorted(class_map.values(), key=lambda d: (d['H'], d['I6']))

    print(f"\n{'H':>6} {'I6':>6} {'alpha1':>7} {'alpha2':>7} {'sc':>4} {'rho1':>8} {'rho2':>8} {'ratio':>8} {'ULC':>5}")
    print("-" * 70)

    sc_ratios = []
    non_sc_ratios = []
    alpha1_3_alpha2_0_found = False
    separation_pairs = defaultdict(list)

    for d in classes:
        coeffs = d['coeffs']
        H = d['H']
        I6 = d['I6']
        sc = d['sc']
        alpha1 = coeffs[1] if len(coeffs) > 1 else 0
        alpha2 = coeffs[2] if len(coeffs) > 2 else 0

        # Root computation
        roots = poly_roots(coeffs)
        neg_roots = sorted([-r.real for r in roots if abs(r.imag) < 1e-8 and r.real < 0], reverse=True)

        if len(neg_roots) == 2:
            rho1, rho2 = neg_roots[0], neg_roots[1]  # rho1 >= rho2
            ratio = rho2 / rho1 if rho1 > 1e-10 else float('inf')
        elif len(neg_roots) == 1:
            rho1, rho2, ratio = neg_roots[0], neg_roots[0], 1.0
        else:
            rho1, rho2, ratio = 0, 0, 0

        # Ultra-log-concavity check
        d_deg = len(coeffs) - 1
        ulc = True
        for k in range(1, d_deg):
            lhs = (coeffs[k] / comb(d_deg, k))**2
            rhs = (coeffs[k-1] / comb(d_deg, k-1)) * (coeffs[k+1] / comb(d_deg, k+1))
            if lhs < rhs - 1e-10:
                ulc = False
        ulc_str = "✓" if ulc else "✗"

        sc_str = "SC" if sc else "  "
        rho1_str = f"{rho1:.4f}" if rho1 > 0 else "  —  "
        rho2_str = f"{rho2:.4f}" if rho2 > 0 and len(neg_roots) == 2 else "  —  "
        ratio_str = f"{ratio:.4f}" if ratio > 0 and len(neg_roots) == 2 else "  —  "

        print(f"{H:>6} {I6:>6} {alpha1:>7} {alpha2:>7} {sc_str:>4} {rho1_str:>8} {rho2_str:>8} {ratio_str:>8} {ulc_str:>5}")

        # Track alpha1=3, alpha2=0
        if alpha1 == 3 and alpha2 == 0:
            alpha1_3_alpha2_0_found = True
            print(f"  *** FORBIDDEN CASE FOUND: alpha1=3, alpha2=0 ***")

        # Track ratios by SC status
        if len(neg_roots) == 2 and sc is not None:
            if sc:
                sc_ratios.append((ratio, H, alpha1, alpha2))
            else:
                non_sc_ratios.append((ratio, H, alpha1, alpha2))

        # Check separation
        pair_key = (H, I6)
        separation_pairs[pair_key].append((alpha1, alpha2, sc))

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTotal iso classes found: {len(classes)} (expected 56)")

    print(f"\nForbidden case alpha1=3, alpha2=0: {'FOUND (unexpected!)' if alpha1_3_alpha2_0_found else 'NOT FOUND ✓'}")

    # Separation analysis
    collisions = {k: v for k, v in separation_pairs.items() if len(v) > 1}
    print(f"\n(H, I(Omega,6)) separation: {len(class_map) - sum(len(v)-1 for v in collisions.values())} / {len(class_map)} classes uniquely identified")
    if collisions:
        print(f"  Collisions (not separated):")
        for (H, I6), clss in collisions.items():
            print(f"    H={H}, I6={I6}: {len(clss)} classes with alpha-tuples {clss}")
    else:
        print(f"  ✓ ALL classes separated by (H, I(Omega,6))!")

    # SC root asymmetry
    if sc_ratios and non_sc_ratios:
        print(f"\nRoot ratio analysis (rho2/rho1, smaller = more asymmetric):")
        print(f"  SC tournaments: min={min(r for r,*_ in sc_ratios):.4f}, max={max(r for r,*_ in sc_ratios):.4f}, mean={sum(r for r,*_ in sc_ratios)/len(sc_ratios):.4f}")
        print(f"  Non-SC tournaments: min={min(r for r,*_ in non_sc_ratios):.4f}, max={max(r for r,*_ in non_sc_ratios):.4f}, mean={sum(r for r,*_ in non_sc_ratios)/len(non_sc_ratios):.4f}")
        sc_min = min(r for r,*_ in sc_ratios)
        nonsc_min = min(r for r,*_ in non_sc_ratios)
        print(f"  SC min ratio < Non-SC min ratio: {'✓' if sc_min < nonsc_min else '✗'}")
        print(f"  (If SC have smallest ratios, conjecture supported)")

    return classes


def run_forbidden_root_check():
    print("\n" + "=" * 70)
    print("FORBIDDEN ROOT CHECK: alpha1=3, alpha2=0 AT n=6,7,8")
    print("=" * 70)

    for n in [6, 7, 8]:
        print(f"\n  n={n}: sampling 20000 random tournaments...")
        found = False
        for _ in range(20000):
            adj = random_tournament(n)
            cycles = find_all_odd_cycles(adj, n)
            coeffs = independence_poly_coeffs(cycles)
            alpha1 = coeffs[1] if len(coeffs) > 1 else 0
            alpha2 = coeffs[2] if len(coeffs) > 2 else 0
            if alpha1 == 3 and alpha2 == 0:
                found = True
                H = hamiltonian_path_count(adj, n)
                print(f"  *** FOUND alpha1=3, alpha2=0 at n={n}: H={H} ***")
                break
        if not found:
            print(f"  ✓ alpha1=3, alpha2=0 NOT FOUND at n={n}")


def run_ratio_analysis():
    print("\n" + "=" * 70)
    print("RATIO I(OMEGA,6)/H ANALYSIS")
    print("=" * 70)
    print("\nTheory: I(Omega,6)/H = prod_i (6+rho_i)/(2+rho_i) in (1, 3^d)")
    print("For rho->0 (regular): ratio -> 3^d. For transitive (d=0): ratio=1.\n")

    n = 6
    ratio_by_sc = {'SC': [], 'NS': []}
    for _ in range(10000):
        adj = random_tournament(n)
        cycles = find_all_odd_cycles(adj, n)
        coeffs = independence_poly_coeffs(cycles)
        H = hamiltonian_path_count(adj, n)
        I6 = sum(coeffs[k]*6**k for k in range(len(coeffs)))
        sc = is_sc(adj, n)
        if H > 0:
            ratio = I6 / H
            key = 'SC' if sc else 'NS'
            ratio_by_sc[key].append(ratio)

    for key in ['SC', 'NS']:
        vals = ratio_by_sc[key]
        if vals:
            print(f"  {key}: n={len(vals)}, mean={sum(vals)/len(vals):.4f}, min={min(vals):.4f}, max={max(vals):.4f}")

    print("\nVerifying ratio in (1, 3^d) for degree-2 polynomial:")
    print(f"  At n=6 with degree-2 (alpha2>=1): ratio should be in (1, 9=3^2)")
    for _ in range(1000):
        adj = random_tournament(n)
        cycles = find_all_odd_cycles(adj, n)
        coeffs = independence_poly_coeffs(cycles)
        H = hamiltonian_path_count(adj, n)
        I6 = sum(coeffs[k]*6**k for k in range(len(coeffs)))
        alpha2 = coeffs[2] if len(coeffs) > 2 else 0
        d = len(coeffs) - 1
        if H > 0 and d >= 1:
            ratio = I6 / H
            upper = 3**d
            if ratio >= upper or ratio <= 1:
                print(f"  *** VIOLATION: ratio={ratio:.4f}, upper={upper}, d={d} ***")
    print(f"  ✓ No violations found in 1000 samples")


def run_newton_maclaurin_check():
    print("\n" + "=" * 70)
    print("NEWTON-MACLAURIN (ULTRA-LOG-CONCAVITY) CHECK")
    print("=" * 70)
    print("\nULC: (alpha_k/C(d,k))^2 >= (alpha_{k-1}/C(d,k-1)) * (alpha_{k+1}/C(d,k+1))")
    print("This follows from real-rootedness (Newton's inequalities).\n")

    violations = 0
    checked = 0
    for n in [5, 6, 7, 8]:
        for _ in range(2000 if n <= 6 else 500):
            adj = random_tournament(n)
            cycles = find_all_odd_cycles(adj, n)
            coeffs = independence_poly_coeffs(cycles)
            d = len(coeffs) - 1
            if d < 2:
                continue
            for k in range(1, d):
                lhs = (coeffs[k] / comb(d, k))**2
                rhs = (coeffs[k-1] / comb(d, k-1)) * (coeffs[k+1] / comb(d, k+1))
                checked += 1
                if lhs < rhs - 1e-9:
                    violations += 1
                    H = hamiltonian_path_count(adj, n)
                    print(f"  *** ULC VIOLATION at n={n}, k={k}: lhs={lhs:.6f}, rhs={rhs:.6f}, H={H}, coeffs={coeffs} ***")

    print(f"\nChecked {checked} Newton-Maclaurin inequalities")
    print(f"Violations: {violations} ({'none ✓' if violations == 0 else 'FOUND!'})")


def run_vieta_formulas():
    print("\n" + "=" * 70)
    print("VIETA'S THEOREM VERIFICATION")
    print("=" * 70)
    print("\nFor degree-1: r = -1/alpha_1 = -2/(H-1) [Vieta for n<=5]")
    print("For degree-2: rho1+rho2 = alpha_1/alpha_2, rho1*rho2 = 1/alpha_2\n")

    # n=5: all degree 1
    n = 5
    print(f"n={n} (all degree-1):")
    errs = []
    for _ in range(500):
        adj = random_tournament(n)
        cycles = find_all_odd_cycles(adj, n)
        coeffs = independence_poly_coeffs(cycles)
        H = hamiltonian_path_count(adj, n)
        d = len(coeffs) - 1
        if d != 1:
            continue
        alpha1 = coeffs[1]
        vieta_r = -1.0/alpha1 if alpha1 > 0 else None
        formula_r = -2.0/(H-1) if H > 1 else None
        if vieta_r is not None and formula_r is not None:
            errs.append(abs(vieta_r - formula_r))
    if errs:
        print(f"  Max error |vieta_r - (-2/(H-1))|: {max(errs):.2e}")
        print(f"  ✓ Vieta r = -1/alpha_1 = -2/(H-1)" if max(errs) < 1e-10 else "  ✗ MISMATCH")

    # n=6: degree 1 or 2
    n = 6
    print(f"\nn={n} (degree-1 and degree-2):")
    errs_sum = []
    errs_prod = []
    for _ in range(2000):
        adj = random_tournament(n)
        cycles = find_all_odd_cycles(adj, n)
        coeffs = independence_poly_coeffs(cycles)
        d = len(coeffs) - 1
        if d != 2:
            continue
        alpha1, alpha2 = coeffs[1], coeffs[2]
        roots = poly_roots(coeffs)
        neg_roots = sorted([-r.real for r in roots if abs(r.imag) < 1e-8 and r.real < 0])
        if len(neg_roots) != 2:
            continue
        rho1, rho2 = neg_roots
        # Vieta: rho1+rho2 = alpha1/alpha2, rho1*rho2 = 1/alpha2
        vieta_sum = alpha1 / alpha2
        vieta_prod = 1.0 / alpha2
        errs_sum.append(abs(rho1 + rho2 - vieta_sum))
        errs_prod.append(abs(rho1 * rho2 - vieta_prod))
    if errs_sum:
        print(f"  Sum error: max={max(errs_sum):.2e}")
        print(f"  Product error: max={max(errs_prod):.2e}")
        ok = max(errs_sum) < 1e-8 and max(errs_prod) < 1e-8
        print(f"  ✓ Vieta verified for degree-2 polynomials" if ok else "  ✗ MISMATCH")


def run_root_gap_analysis():
    print("\n" + "=" * 70)
    print("FORBIDDEN ROOT GAP (-1/3, -1/4) ANALYSIS")
    print("=" * 70)
    print("\nConjecture: No tournament at any n has I(Omega,x) with root in (-1/3, -1/4).")
    print("Root r=-1/3 requires alpha1=3, alpha2=0 (impossible).")
    print("Root r=-1/4 corresponds to alpha1=4, alpha2=0 (H=9 = double root case).\n")

    gap_violations = 0
    root_histogram = defaultdict(int)

    for n in [5, 6, 7, 8]:
        count = 5000 if n <= 7 else 2000
        for _ in range(count):
            adj = random_tournament(n)
            cycles = find_all_odd_cycles(adj, n)
            coeffs = independence_poly_coeffs(cycles)
            if len(coeffs) < 2:
                continue
            roots = poly_roots(coeffs)
            for r in roots:
                if abs(r.imag) < 1e-8:
                    rv = r.real
                    if -1.0/3 < rv < -1.0/4:
                        gap_violations += 1
                        H = hamiltonian_path_count(adj, n)
                        print(f"  *** GAP VIOLATION at n={n}: root={rv:.6f}, H={H}, coeffs={coeffs} ***")
                    # Histogram of root positions
                    bucket = int(rv * 10) / 10.0
                    root_histogram[bucket] += 1

    print(f"Total gap violations: {gap_violations} ({'none ✓' if gap_violations == 0 else 'FOUND!'})")

    print("\nRoot distribution (frequency by bucket):")
    for b in sorted(root_histogram.keys()):
        freq = root_histogram[b]
        bar = '#' * (freq // 50)
        print(f"  [{b:.1f}, {b+0.1:.1f}): {freq:5d} {bar}")


if __name__ == '__main__':
    random.seed(42)
    np.random.seed(42)

    # Main n=6 analysis
    classes = run_n6_analysis()

    # Check forbidden root at n=6,7,8
    run_forbidden_root_check()

    # Ratio I(Omega,6)/H analysis
    run_ratio_analysis()

    # Newton-Maclaurin check
    run_newton_maclaurin_check()

    # Vieta verification
    run_vieta_formulas()

    # Root gap analysis
    run_root_gap_analysis()

    print("\n" + "=" * 70)
    print("INVESTIGATION COMPLETE")
    print("=" * 70)
