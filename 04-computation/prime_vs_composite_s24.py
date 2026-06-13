#!/usr/bin/env python3
"""
prime_vs_composite_s24.py — opus-2026-04-05-S24

SYSTEMATIC comparison of tournament invariants at prime vs composite n.

Tang-Yau's stability theorem works for large primes. THM-125 (constant
symbol matrix) is proved for circulant tournaments. The δ parameter
(β_top = (n-1) + δ) is 0 for prime n=7 and 2 for composite n=9.

Key question: WHERE EXACTLY does the prime/composite distinction show up?

We compare across n=3..8 (exhaustive) and n=9 (sampled):
1. H-spectrum: distribution shape, mean, variance, skewness
2. H-maximizer: is it unique? circulant? automorphism group?
3. Iso class count and metagraph: density, spectral gap
4. Cycle structure: c₃ distribution, c₅ distribution
5. Score sequence diversity
6. Self-complementary fraction

For each, we ask: is there a systematic difference between prime and composite n?

n values: 3 (prime), 4 (2²), 5 (prime), 6 (2·3), 7 (prime), 8 (2³), 9 (3²)
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
from math import factorial, gcd
import random
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
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

def count_ham_paths(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))

def score_c3(A, n):
    total = n * (n-1) * (n-2) // 6
    for v in range(n):
        s = sum(A[v])
        total -= s * (s-1) // 2
    return total

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def canonical_form(A, n):
    best = None
    for perm in permutations(range(n)):
        bits = 0
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] == 1:
                    bits |= (1 << idx)
                idx += 1
        if best is None or bits < best:
            best = bits
    return best

def is_circulant(A, n):
    """Check if tournament is circulant (shift-invariant under cyclic permutation)."""
    # Check: A[i][j] = A[(i+1)%n][(j+1)%n] for all i,j
    for i in range(n):
        for j in range(n):
            if A[i][j] != A[(i+1) % n][(j+1) % n]:
                return False
    return True

def automorphism_count(A, n):
    """Count automorphisms (permutations preserving adjacency)."""
    count = 0
    for perm in permutations(range(n)):
        is_auto = True
        for i in range(n):
            for j in range(i+1, n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    is_auto = False
                    break
            if not is_auto:
                break
        if is_auto:
            count += 1
    return count

def is_self_complementary(A, n):
    """Check if T ≅ T^op."""
    Ac = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                Ac[i][j] = 1 - A[i][j]
    return canonical_form(A, n) == canonical_form(Ac, n)

def is_strongly_connected(A, n):
    """Check strong connectivity via DFS."""
    def reachable(adj, start):
        visited = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            for u in range(n):
                if adj[v][u] and u not in visited:
                    stack.append(u)
        return visited

    if len(reachable(A, 0)) < n:
        return False
    # Check reverse
    AT = [[A[j][i] for j in range(n)] for i in range(n)]
    return len(reachable(AT, 0)) == n

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0: return False
    return True

def analyze_n(n, exhaustive=True, num_samples=50000):
    m = n * (n - 1) // 2
    total = 1 << m

    # Collect data
    H_values = []
    c3_values = []
    score_seqs = Counter()
    iso_classes = {}
    sc_count = 0
    circulant_count = 0
    strongly_connected_count = 0
    H_max = 0
    H_max_tournaments = []

    t0 = time.time()

    if exhaustive:
        source = range(total)
        count_target = total
    else:
        source = (random.randint(0, total-1) for _ in range(num_samples))
        count_target = num_samples

    processed = 0
    for bits in source:
        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        c3 = score_c3(A, n)
        ss = score_seq(A, n)

        H_values.append(H)
        c3_values.append(c3)
        score_seqs[ss] += 1

        if H > H_max:
            H_max = H
            H_max_tournaments = [bits]
        elif H == H_max:
            H_max_tournaments.append(bits)

        # For exhaustive small n, compute more
        if exhaustive and n <= 7:
            cf = canonical_form(A, n)
            if cf not in iso_classes:
                iso_classes[cf] = {
                    'H': H, 'c3': c3, 'score': ss,
                    'count': 0,
                    'sc': is_self_complementary(A, n),
                    'circ': is_circulant(A, n),
                    'sc_connected': is_strongly_connected(A, n),
                }
            iso_classes[cf]['count'] += 1

        processed += 1
        if processed % 10000 == 0:
            print(f"  {processed}/{count_target}...    ", end='\r')

    elapsed = time.time() - t0

    H_arr = np.array(H_values, dtype=float)
    c3_arr = np.array(c3_values, dtype=float)

    # H-spectrum statistics
    H_counter = Counter(H_values)
    distinct_H = sorted(H_counter.keys())
    mean_H = np.mean(H_arr)
    var_H = np.var(H_arr)
    std_H = np.std(H_arr)
    skew_H = float(np.mean((H_arr - mean_H)**3) / std_H**3) if std_H > 0 else 0
    kurt_H = float(np.mean((H_arr - mean_H)**4) / std_H**4 - 3) if std_H > 0 else 0

    # Theoretical mean for random tournament
    mean_H_theory = factorial(n) / 2**(n-1)

    # H-maximizer analysis
    H_max_A = bits_to_adj(H_max_tournaments[0], n)
    H_max_circ = is_circulant(H_max_A, n) if len(H_max_tournaments) > 0 else None
    H_max_score = score_seq(H_max_A, n)
    H_max_c3 = score_c3(H_max_A, n)

    # Aut group of H-maximizer (only for small n)
    H_max_aut = None
    if n <= 7:
        H_max_aut = automorphism_count(H_max_A, n)

    # Iso class analysis (only for exhaustive small n)
    num_iso = len(iso_classes) if iso_classes else "N/A"
    num_sc = sum(1 for v in iso_classes.values() if v.get('sc', False)) if iso_classes else "N/A"
    num_circ = sum(1 for v in iso_classes.values() if v.get('circ', False)) if iso_classes else "N/A"

    # Score diversity
    num_score_seqs = len(score_seqs)
    # Regular tournament exists only at odd n
    has_regular = any(s == tuple([n//2]*n) for s in score_seqs) if n % 2 == 1 else False

    # c3 statistics
    mean_c3 = np.mean(c3_arr)
    max_c3 = int(np.max(c3_arr))
    # Theoretical: E[c3] = C(n,3)/4 for random tournament
    mean_c3_theory = n*(n-1)*(n-2)/24

    return {
        'n': n,
        'is_prime': is_prime(n),
        'prime_type': 'PRIME' if is_prime(n) else f'COMPOSITE ({factorize(n)})',
        'total': count_target,
        'elapsed': elapsed,
        # H-spectrum
        'H_max': H_max,
        'H_min': min(H_values),
        'mean_H': mean_H,
        'mean_H_theory': mean_H_theory,
        'std_H': std_H,
        'skew_H': skew_H,
        'kurt_H': kurt_H,
        'distinct_H': len(distinct_H),
        'H_values': distinct_H,
        # H-maximizer
        'H_max_count': len(H_max_tournaments),
        'H_max_circ': H_max_circ,
        'H_max_score': H_max_score,
        'H_max_c3': H_max_c3,
        'H_max_aut': H_max_aut,
        # Iso classes
        'num_iso': num_iso,
        'num_sc': num_sc,
        'num_circ': num_circ,
        # Scores
        'num_score_seqs': num_score_seqs,
        'has_regular': has_regular,
        # Cycles
        'mean_c3': mean_c3,
        'mean_c3_theory': mean_c3_theory,
        'max_c3': max_c3,
        # Correlations
        'corr_c3_H': float(np.corrcoef(c3_arr, H_arr)[0,1]) if std_H > 0 else 0,
    }

def factorize(n):
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return '·'.join(str(f) for f in factors)

# ================================================================
# Run analysis
# ================================================================

print("=" * 70)
print("PRIME vs COMPOSITE: SYSTEMATIC COMPARISON")
print("=" * 70)

results = {}
for n in range(3, 9):
    m = n * (n-1) // 2
    total = 1 << m
    exhaustive = total <= 200000  # exhaustive for n <= 7
    print(f"\nn = {n} ({'PRIME' if is_prime(n) else 'COMPOSITE ' + factorize(n)}), "
          f"m = {m}, {'exhaustive' if exhaustive else 'sampled'}")
    results[n] = analyze_n(n, exhaustive=exhaustive,
                           num_samples=50000 if n == 8 else 10000)
    r = results[n]
    print(f"  H: [{r['H_min']}, {r['H_max']}], mean={r['mean_H']:.2f} (theory {r['mean_H_theory']:.2f}), "
          f"std={r['std_H']:.2f}, skew={r['skew_H']:.2f}, kurt={r['kurt_H']:.2f}")
    print(f"  H-spectrum: {r['distinct_H']} distinct values")
    print(f"  H-max: {r['H_max']}, count={r['H_max_count']}, "
          f"circ={r['H_max_circ']}, score={r['H_max_score']}, c3={r['H_max_c3']}")
    if r['H_max_aut'] is not None:
        print(f"  H-max |Aut| = {r['H_max_aut']}")
    print(f"  Iso classes: {r['num_iso']}, SC: {r['num_sc']}, Circulant: {r['num_circ']}")
    print(f"  c₃: mean={r['mean_c3']:.2f} (theory {r['mean_c3_theory']:.2f}), "
          f"max={r['max_c3']}, corr(c₃,H)={r['corr_c3_H']:.4f}")
    print(f"  Score sequences: {r['num_score_seqs']}, regular exists: {r['has_regular']}")

# ================================================================
# Side-by-side comparison
# ================================================================

print("\n\n" + "=" * 70)
print("SIDE-BY-SIDE: PRIME vs COMPOSITE")
print("=" * 70)

header = f"{'n':>3} {'type':>10} | {'H_max':>6} {'|Aut|':>5} {'circ':>5} | "
header += f"{'E[H]':>7} {'σ(H)':>6} {'skew':>6} {'kurt':>6} | "
header += f"{'#iso':>5} {'#SC':>4} {'r(c3,H)':>7} | {'#scores':>7} {'reg':>3}"
print(f"\n  {header}")
print(f"  {'-'*len(header)}")

for n in sorted(results.keys()):
    r = results[n]
    ptype = "PRIME" if r['is_prime'] else "comp"
    circ = "Y" if r['H_max_circ'] else "N" if r['H_max_circ'] is not None else "?"
    aut = str(r['H_max_aut']) if r['H_max_aut'] is not None else "?"
    iso = str(r['num_iso']) if r['num_iso'] != "N/A" else "?"
    sc = str(r['num_sc']) if r['num_sc'] != "N/A" else "?"
    reg = "Y" if r['has_regular'] else "N"
    print(f"  {n:>3} {ptype:>10} | {r['H_max']:>6} {aut:>5} {circ:>5} | "
          f"{r['mean_H']:>7.1f} {r['std_H']:>6.1f} {r['skew_H']:>6.2f} {r['kurt_H']:>6.2f} | "
          f"{iso:>5} {sc:>4} {r['corr_c3_H']:>7.4f} | {r['num_score_seqs']:>7} {reg:>3}")

# Now print prime-specific vs composite-specific statistics
print("\n\n  PRIME n averages:")
prime_data = {n: r for n, r in results.items() if r['is_prime']}
comp_data = {n: r for n, r in results.items() if not r['is_prime']}

for label, data in [("PRIME", prime_data), ("COMPOSITE", comp_data)]:
    if not data:
        continue
    ns = sorted(data.keys())
    print(f"\n  {label} n = {ns}:")
    # H_max / n! trend
    for n in ns:
        r = data[n]
        ratio_max = r['H_max'] / factorial(n)
        ratio_mean = r['mean_H'] / r['mean_H_theory']
        aut = r['H_max_aut'] if r['H_max_aut'] else "?"
        print(f"    n={n}: H_max/n! = {ratio_max:.6f}, |Aut(H_max)| = {aut}, "
              f"H_max = {r['H_max']}")

# Key ratios
print("\n\n  KEY DIAGNOSTIC RATIOS:")
print(f"  {'n':>3} {'type':>6} | {'H_max/n!':>10} | {'H_max*|Aut|':>12} | {'H_max/(n!/2^(n-1))':>18} | {'skew':>6}")
print(f"  {'-'*3}-{'-'*6}-|-{'-'*10}-|-{'-'*12}-|-{'-'*18}-|-{'-'*6}")
for n in sorted(results.keys()):
    r = results[n]
    ptype = "PRIME" if r['is_prime'] else "comp"
    ratio1 = r['H_max'] / factorial(n)
    aut = r['H_max_aut'] if r['H_max_aut'] else 1
    prod = r['H_max'] * aut
    ratio2 = r['H_max'] / r['mean_H_theory']
    print(f"  {n:>3} {ptype:>6} | {ratio1:>10.6f} | {prod:>12} | {ratio2:>18.4f} | {r['skew_H']:>6.2f}")

# H-spectrum shape comparison
print("\n\n  H-SPECTRUM SHAPE at prime vs composite:")
for n in sorted(results.keys()):
    r = results[n]
    ptype = "PRIME" if r['is_prime'] else "comp"
    vals = r['H_values']
    if len(vals) <= 20:
        print(f"  n={n} ({ptype}): H in {vals}")
    else:
        print(f"  n={n} ({ptype}): H in [{vals[0]}..{vals[-1]}], {len(vals)} values, "
              f"gaps at: {[vals[i+1]-vals[i] for i in range(min(15, len(vals)-1))]}")
