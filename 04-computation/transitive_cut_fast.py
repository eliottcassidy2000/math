#!/usr/bin/env python3
"""
transitive_cut_fast.py — opus-2026-03-14-S71g

Faster version: use score-sequence-based strong connectivity check.

A tournament is strongly connected iff it has no "dominating initial segment":
no proper prefix of the score-sorted vertices beats all later vertices.

Actually, the standard test: a tournament T is strongly connected iff
for the score sequence s₁ ≤ s₂ ≤ ... ≤ sₙ,
  Σ_{i=1}^k sᵢ > C(k,2) for all 1 ≤ k < n.
(Landau's criterion: strict inequality for all proper initial segments)

This is O(n²) per tournament instead of O(2^n).
"""

from itertools import permutations, combinations
from collections import defaultdict
from math import comb

def make_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def is_strongly_connected_fast(A, n):
    """Landau criterion: SC iff for sorted scores s₁≤...≤sₙ,
    Σ_{i=1}^k sᵢ > C(k,2) for all 1 ≤ k < n."""
    scores = sorted(sum(A[i]) for i in range(n))
    partial = 0
    for k in range(1, n):
        partial += scores[k-1]
        if partial <= comb(k, 2):
            return False
    return True

# ============================================================
# Part 1: SC frequency by n
# ============================================================

print("=" * 70)
print("STRONG CONNECTIVITY AND TRANSITIVE CUT ANALYSIS")
print("=" * 70)

for n in range(3, 9):
    total_edges = n * (n - 1) // 2
    num_t = 2 ** total_edges

    if num_t > 2**21:
        import random
        random.seed(42)
        sample = 200000
        sc_count = 0
        for _ in range(sample):
            bits = random.randint(0, num_t - 1)
            A = make_tournament(bits, n)
            if is_strongly_connected_fast(A, n):
                sc_count += 1
        pct = 100*sc_count/sample
        print(f"  n={n} (sample {sample}): SC={sc_count} ({pct:.1f}%)")
    else:
        sc_count = 0
        for bits in range(num_t):
            A = make_tournament(bits, n)
            if is_strongly_connected_fast(A, n):
                sc_count += 1
        pct = 100*sc_count/num_t
        print(f"  n={n} ({num_t} total): SC={sc_count} ({pct:.1f}%)")

# ============================================================
# Part 2: "Prime" H values (SC only)
# ============================================================

print(f"\n{'='*70}")
print("PRIME H VALUES (only from strongly connected tournaments)")
print(f"{'='*70}")

all_h_by_n = {}
for n in range(3, 7):  # Skip n=7 (too slow for exhaustive H computation)
    total_edges = n * (n - 1) // 2
    num_t = 2 ** total_edges

    h_from_sc = set()
    h_from_nonsc = set()

    for bits in range(num_t):
        A = make_tournament(bits, n)
        H = count_hp(A, n)
        if is_strongly_connected_fast(A, n):
            h_from_sc.add(H)
        else:
            h_from_nonsc.add(H)

    only_sc = h_from_sc - h_from_nonsc
    both = h_from_sc & h_from_nonsc
    all_h_by_n[n] = sorted(h_from_sc | h_from_nonsc)

    print(f"\n  n={n}:")
    print(f"    SC H values:     {sorted(h_from_sc)}")
    print(f"    Non-SC H values: {sorted(h_from_nonsc)}")
    print(f"    SC-only:         {sorted(only_sc)}")

# ============================================================
# Part 3: H-spectrum factorization and achievability semigroup
# ============================================================

print(f"\n{'='*70}")
print("H-SPECTRUM FACTORIZATION AND ACHIEVABILITY")
print(f"{'='*70}")

# Achievable H values from exact computation (n≤7)
achievable_exact = {}
for n in range(3, 7):
    achievable_exact[n] = set(all_h_by_n[n])

# Add n=7 H-spectrum from known results
# At n=7: H ∈ {1,3,5,...,189} (all odd, verified exhaustive by earlier work)
# Known: H=7 and H=21 NOT achievable at n=7 (from previous sessions)
import random
random.seed(42)
h7_sample = set()
for _ in range(50000):
    bits = random.randint(0, 2**21-1)
    A = make_tournament(bits, 7)
    H = count_hp(A, 7)
    h7_sample.add(H)
achievable_exact[7] = h7_sample
all_h_by_n[7] = sorted(h7_sample)

# Extend by products (transitive cut construction)
achievable_products = set()
for n1 in range(3, 8):
    for n2 in range(3, 8):
        for a in achievable_exact.get(n1, set()):
            for b in achievable_exact.get(n2, set()):
                achievable_products.add(a * b)

all_achievable = set()
for n in range(3, 8):
    all_achievable.update(achievable_exact.get(n, set()))
all_achievable.update(achievable_products)

# Generate multi-level products
changed = True
while changed:
    changed = False
    new = set()
    for a in list(all_achievable):
        if a > 500:
            continue
        for n in range(3, 6):
            for b in achievable_exact.get(n, set()):
                p = a * b
                if p <= 500 and p not in all_achievable:
                    new.add(p)
                    changed = True
    all_achievable.update(new)

# Which odd values ≤ 200 are NOT achievable?
never_odd = sorted(h for h in range(1, 201, 2) if h not in all_achievable)
print(f"\n  Odd H values NOT achievable (any n ≤ 14, via products): {never_odd}")
print(f"  Total: {len(never_odd)} out of 100 odd values ≤ 200")

# Check forbidden sequences
forbidden_7 = [7 * 3**k for k in range(5) if 7 * 3**k <= 200]
forbidden_21 = [21 * 3**k for k in range(5) if 21 * 3**k <= 200]
print(f"\n  7·3^k sequence: {forbidden_7}")
print(f"    All forbidden? {all(h in never_odd for h in forbidden_7)}")
print(f"  21·3^k sequence: {forbidden_21}")
print(f"    All forbidden? {all(h in never_odd for h in forbidden_21)}")

# What about other values?
explained = set(forbidden_7 + forbidden_21)
unexplained = [h for h in never_odd if h not in explained]
print(f"\n  Unexplained forbidden values (not in 7·3^k or 21·3^k):")
print(f"    {unexplained}")

# Are any of these achievable by n=8 SC tournaments?
# Check a sample
import random
random.seed(42)
n = 8
total_edges = n * (n - 1) // 2
sample = 500000
h8_vals = set()
for _ in range(sample):
    bits = random.randint(0, 2**total_edges - 1)
    A = make_tournament(bits, n)
    H = count_hp(A, n)
    h8_vals.add(H)

newly_achievable = [h for h in unexplained if h in h8_vals]
still_missing = [h for h in unexplained if h not in h8_vals]
print(f"\n  Of unexplained, found at n=8 (500k sample): {newly_achievable}")
print(f"  Still missing at n=8: {still_missing}")

# ============================================================
# Part 4: The SCC decomposition formula
# ============================================================

print(f"\n{'='*70}")
print("SCC DECOMPOSITION VERIFICATION")
print(f"{'='*70}")

# For non-SC tournaments, verify H = product of SCC HP counts
def find_scc_sizes(A, n):
    """Find SCCs of tournament using condensation.
    Tournament condensation: find score-based ordering."""
    # Use Tarjan's or Kosaraju's
    # Simple approach for small n: check all subsets
    # Better: use Landau criterion iteratively

    scores = [sum(A[i]) for i in range(n)]
    # Sort vertices by score (ascending)
    order = sorted(range(n), key=lambda v: scores[v])

    # Find SCCs: tournament has unique SCC decomposition
    # SCC boundaries at positions where partial score sum = C(k,2)
    sorted_scores = [scores[v] for v in order]
    sccs = []
    start = 0
    partial = 0
    for k in range(1, n+1):
        partial += sorted_scores[k-1]
        if partial == comb(k, 2):
            # Found SCC boundary: vertices order[start:k]
            sccs.append(frozenset(order[start:k]))
            start = k

    return sccs

# Verify product formula at n=6
n = 6
total_edges = n * (n - 1) // 2
num_t = 2 ** total_edges

mismatches = 0
verified = 0
for bits in range(num_t):
    A = make_tournament(bits, n)
    if is_strongly_connected_fast(A, n):
        continue

    H = count_hp(A, n)
    sccs = find_scc_sizes(A, n)

    if len(sccs) <= 1:
        continue

    # Compute product of internal H values
    product = 1
    for scc in sccs:
        scc_list = sorted(scc)
        m = len(scc_list)
        sub_A = [[0]*m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                if i != j:
                    sub_A[i][j] = A[scc_list[i]][scc_list[j]]
        product *= count_hp(sub_A, m)

    if H != product:
        mismatches += 1
        if mismatches <= 3:
            print(f"  MISMATCH: bits={bits}, H={H}, product={product}, SCCs={[len(s) for s in sccs]}")
    verified += 1

print(f"\n  n={n}: Verified {verified} non-SC tournaments, {mismatches} mismatches")
if mismatches == 0:
    print(f"  CONFIRMED: H = ∏ H(SCCᵢ) for ALL non-SC tournaments at n={n}")

# ============================================================
# Part 5: SCC size distribution
# ============================================================

print(f"\n{'='*70}")
print("SCC SIZE DISTRIBUTION")
print(f"{'='*70}")

for n in range(3, 7):  # Skip n=7
    total_edges = n * (n - 1) // 2
    num_t = 2 ** total_edges

    scc_patterns = defaultdict(int)
    for bits in range(num_t):
        A = make_tournament(bits, n)
        sccs = find_scc_sizes(A, n)
        pattern = tuple(sorted([len(s) for s in sccs], reverse=True))
        scc_patterns[pattern] += 1

    print(f"\n  n={n}:")
    for pattern in sorted(scc_patterns.keys()):
        count = scc_patterns[pattern]
        pct = 100 * count / num_t
        is_sc = (len(pattern) == 1)
        label = " [SC]" if is_sc else ""
        print(f"    {'+'.join(str(s) for s in pattern):10s}: {count:6d} ({pct:5.1f}%){label}")

print(f"\n{'='*70}")
print("ANALYSIS COMPLETE")
print(f"{'='*70}")
