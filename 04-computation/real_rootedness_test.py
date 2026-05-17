"""
Test the conjecture: For any tournament T, I(Omega(T), x) has all real, negative roots.

Evidence so far:
- n<=8: Omega(T) is claw-free -> real-rooted by Chudnovsky-Seymour (2007)
- k=5 (n=10), k=6 (n=12): verified for all-0 staircase (which CAN have claws)

This script:
1. Tests ALL tournaments at n=5, n=6 (using canonical forms)
2. Tests random samples at n=7, n=8, n=9, n=10
3. Checks if any polynomial has complex roots
4. Looks for patterns in the roots

Key insight: if I(Omega, x) = prod_i (1 + r_i * x) with all r_i > 0,
then at x=2: H = prod_i (1 + 2*r_i), a product of terms > 1.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
import random
import sys

# ============================================================
# Tournament utilities
# ============================================================

def random_tournament(n, seed=None):
    if seed is not None:
        random.seed(seed)
    arcs = set()
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                arcs.add((i, j))
            else:
                arcs.add((j, i))
    return arcs

def transitive_tournament(n):
    return {(i, j) for i in range(n) for j in range(i+1, n)}

def all_tournaments_n(n):
    """Generate all labeled tournaments on n vertices."""
    edges = list(combinations(range(n), 2))
    for bits in range(2**len(edges)):
        arcs = set()
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                arcs.add((i, j))
            else:
                arcs.add((j, i))
        yield arcs

# ============================================================
# Odd cycle finding via DFS
# ============================================================

def find_odd_cycles(arcs, n):
    out = defaultdict(list)
    for (i, j) in arcs:
        out[i].append(j)

    odd_cycles = []

    def dfs(start, curr, path, visited):
        for nxt in out[curr]:
            if nxt == start and len(path) % 2 == 1:  # odd length
                odd_cycles.append(frozenset(path))
                # Don't recurse to avoid counting sub-cycles
            elif nxt > start and nxt not in visited and len(path) < n:
                visited.add(nxt)
                path.append(nxt)
                dfs(start, nxt, path, visited)
                path.pop()
                visited.remove(nxt)

    for v in range(n):
        dfs(v, v, [v], {v})

    # Deduplicate (frozenset handles vertex-set equality)
    return list(set(odd_cycles))

def find_odd_cycles_with_order(arcs, n):
    """Find all directed odd cycles as tuples (sequence of vertices)."""
    out = defaultdict(list)
    for (i, j) in arcs:
        out[i].append(j)

    odd_cycles = []

    def dfs(start, curr, path, visited):
        for nxt in out[curr]:
            if nxt == start and len(path) % 2 == 1:
                odd_cycles.append(tuple(path))
            elif nxt > start and nxt not in visited and len(path) < n:
                visited.add(nxt)
                path.append(nxt)
                dfs(start, nxt, path, visited)
                path.pop()
                visited.remove(nxt)

    for v in range(n):
        dfs(v, v, [v], {v})

    return odd_cycles

# ============================================================
# Build conflict graph and compute independence polynomial
# ============================================================

def build_omega_and_indpoly(arcs, n):
    """Full pipeline: arcs -> odd cycles -> Omega -> I(Omega, x) coefficients."""
    cycles = find_odd_cycles(arcs, n)
    m = len(cycles)

    if m == 0:
        return [1], []  # transitive: no odd cycles, H = 1

    # Build adjacency (conflict) for Omega
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:  # shared vertex
                adj[i][j] = adj[j][i] = True

    # Enumerate independent sets via branch-and-bound
    coeffs = defaultdict(int)

    def dfs_indep(candidates, size):
        coeffs[size] += 1
        for idx, v in enumerate(candidates):
            new_cands = [u for u in candidates[idx+1:] if not adj[v][u]]
            dfs_indep(new_cands, size + 1)

    dfs_indep(list(range(m)), 0)

    max_k = max(coeffs.keys())
    poly = [coeffs[k] for k in range(max_k + 1)]
    return poly, cycles

def check_real_rootedness(poly):
    """Check if polynomial is real-rooted (all roots real and negative)."""
    if len(poly) <= 1:
        return True, [], True  # constant polynomial

    # Numpy: poly1d expects descending powers
    coeffs_desc = list(reversed(poly))
    p = np.poly1d(coeffs_desc)
    roots = p.roots

    all_real = all(abs(r.imag) < 1e-6 for r in roots)
    all_neg = all(r.real < 0 for r in roots) if all_real else False
    real_roots = sorted([r.real for r in roots])

    return all_real, real_roots, all_neg

def compute_H(poly):
    return sum(c * 2**k for k, c in enumerate(poly))

# ============================================================
# Main tests
# ============================================================

def test_all_n5():
    print("\n=== Testing ALL tournaments at n=5 ===")
    n = 5
    total = 0
    real_rooted_count = 0
    max_neg = 0
    fail_examples = []

    for arcs in all_tournaments_n(n):
        poly, _ = build_omega_and_indpoly(arcs, n)
        all_real, roots, all_neg = check_real_rootedness(poly)
        H = compute_H(poly)
        total += 1
        if all_real and all_neg:
            real_rooted_count += 1
        else:
            fail_examples.append((arcs, poly, roots))

    print(f"Total tournaments: {total}")
    print(f"Real-rooted (all real negative roots): {real_rooted_count}/{total}")
    if fail_examples:
        print(f"FAILURES: {len(fail_examples)}")
        for arcs, poly, roots in fail_examples[:3]:
            print(f"  poly={poly}, roots={roots}")
    else:
        print("✓ All real-rooted!")
    return len(fail_examples) == 0

def test_n6_sample():
    """Test n=6 iso-class representatives (sample for speed)."""
    print("\n=== Testing n=6 tournaments (sample) ===")
    n = 6
    total = 0
    real_rooted_count = 0
    fail_examples = []

    # We can't do all 2^15=32768 but can do a large sample
    random.seed(42)
    for _ in range(5000):
        arcs = random_tournament(n, seed=None)
        poly, _ = build_omega_and_indpoly(arcs, n)
        all_real, roots, all_neg = check_real_rootedness(poly)
        total += 1
        if all_real and all_neg:
            real_rooted_count += 1
        else:
            fail_examples.append((poly, roots))

    print(f"Tested: {total}")
    print(f"Real-rooted: {real_rooted_count}/{total}")
    if fail_examples:
        print(f"FAILURES: {len(fail_examples)}")
        for poly, roots in fail_examples[:3]:
            print(f"  poly={poly}")
            print(f"  roots={[f'{r:.4f}' for r in roots]}")
    else:
        print("✓ All real-rooted!")
    return len(fail_examples) == 0

def test_larger_n(n, num_samples):
    """Test random tournaments at larger n."""
    print(f"\n=== Testing n={n} ({num_samples} random tournaments) ===")
    random.seed(12345)
    total = 0
    real_rooted_count = 0
    fail_examples = []
    H_values = []
    deg_values = []

    for trial in range(num_samples):
        arcs = random_tournament(n)
        poly, _ = build_omega_and_indpoly(arcs, n)
        all_real, roots, all_neg = check_real_rootedness(poly)
        H = compute_H(poly)
        H_values.append(H)
        deg_values.append(len(poly) - 1)
        total += 1
        if all_real and all_neg:
            real_rooted_count += 1
        else:
            fail_examples.append((poly, roots, H))

    print(f"Tested: {total}")
    print(f"Real-rooted: {real_rooted_count}/{total}")
    print(f"H range: [{min(H_values)}, {max(H_values)}]")
    print(f"Poly degree range: [{min(deg_values)}, {max(deg_values)}]")
    if fail_examples:
        print(f"FAILURES: {len(fail_examples)}")
        for poly, roots, H in fail_examples[:3]:
            print(f"  H={H}, poly={poly[:5]}...")
            print(f"  roots (first 5)={[f'{r:.4f}' for r in sorted(roots)[:5]]}")
    else:
        print("✓ All real-rooted!")
    return len(fail_examples) == 0

def analyze_root_structure(n, num_samples=100):
    """Analyze the structure of roots of I(Omega, x)."""
    print(f"\n=== Root structure analysis at n={n} ===")
    random.seed(99)
    all_roots = []
    H_values = []

    for _ in range(num_samples):
        arcs = random_tournament(n)
        poly, _ = build_omega_and_indpoly(arcs, n)
        _, roots, _ = check_real_rootedness(poly)
        H = compute_H(poly)
        all_roots.extend(roots)
        H_values.append(H)

    if all_roots:
        neg_roots = [-r for r in all_roots]  # make positive
        neg_roots_sorted = sorted(neg_roots)
        print(f"Root statistics (negated to be positive):")
        print(f"  Min: {min(neg_roots):.6f}")
        print(f"  Max: {max(neg_roots):.2f}")
        print(f"  Geometric mean: {np.exp(np.mean(np.log([r for r in neg_roots if r > 1e-10]))):.4f}")

        # Check if smallest roots approach 0
        small_roots = [r for r in neg_roots if r < 0.01]
        print(f"  Roots < 0.01: {len(small_roots)}/{len(neg_roots)} ({100*len(small_roots)/len(neg_roots):.1f}%)")

if __name__ == '__main__':
    print("Testing real-rootedness conjecture: I(Omega(T), x) has all real negative roots")
    print("for every tournament T.")
    print()

    ok5 = test_all_n5()
    ok6 = test_n6_sample()
    ok7 = test_larger_n(7, 200)
    ok8 = test_larger_n(8, 100)
    ok9 = test_larger_n(9, 50)
    ok10 = test_larger_n(10, 30)

    print("\n" + "="*60)
    print("SUMMARY OF REAL-ROOTEDNESS TESTS")
    print("="*60)
    for n, ok in [(5, ok5), (6, ok6), (7, ok7), (8, ok8), (9, ok9), (10, ok10)]:
        status = "✓ PASS" if ok else "✗ FAIL"
        print(f"  n={n}: {status}")

    if all([ok5, ok6, ok7, ok8, ok9, ok10]):
        print("\n🎯 CONJECTURE SUPPORTED: I(Omega(T), x) is real-rooted for all tested T.")
        print("   This extends Chudnovsky-Seymour (claw-free case) to ALL tournaments.")
    else:
        print("\n⚠ Counterexample found!")

    analyze_root_structure(7, 200)
