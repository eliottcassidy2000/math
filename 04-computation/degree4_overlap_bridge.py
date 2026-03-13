#!/usr/bin/env python3
"""
degree4_overlap_bridge.py -- kind-pasteur-2026-03-13-S61

THE BRIDGE: Degree-4 Fourier terms encode the OVERLAP WEIGHT structure.

THM-163 showed:
  H(sigma) = H_0 + H_2 + H_4 + ...
  H_0 = n!/2^{n-1} (constant)
  H_2 = c_2 * n * (m - 2*Var_s) (depends only on score sequence)
  H_4 = ??? (depends on something BEYOND score sequence)

The key question: what is H_4 measuring?

Hypothesis: H_4 encodes the CYCLE OVERLAP structure, i.e., how directed odd
cycles are distributed relative to each other. This is exactly the alpha_2
(disjoint pairs) term in the OCF.

The Vitali connection: The degree-2 term captures the "measurable" part of
tournament structure (score regularity). The degree-4 term captures the
"non-measurable" part — the fine structure of cycle interactions that is
NOT determined by scores alone.

Among tournaments with the SAME score sequence, H varies by H_4 alone.
H_4 encodes the overlap weight matrix of Omega(T).

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def tournament_to_sigma(A, n):
    vec = []
    for i in range(n):
        for j in range(i+1, n):
            vec.append(1 if A[i][j] else -1)
    return tuple(vec)


def count_ham_paths(A, n):
    if n <= 1:
        return 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_cycles_on_subset(A, verts):
    """Count directed Hamiltonian cycles on a vertex subset."""
    k = len(verts)
    if k < 3:
        return 0

    sub_A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            sub_A[i][j] = A[verts[i]][verts[j]]

    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if sub_A[v][w]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]

    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and sub_A[v][0]:
            total += dp[(full, v)]
    return total


def compute_H2(A, n):
    """Compute degree-2 contribution to H using the Z_v formula."""
    c2 = math.factorial(n - 2) / (2 ** (n - 2))
    scores = [sum(A[v]) for v in range(n)]
    half = (n - 1) / 2
    Z = [-2 * (s - half)**2 + half for s in scores]
    return c2 * sum(Z)


# ========================================================================
# ANALYSIS 1: H_4 ISOLATES THE NON-SCORE PART OF H
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: H_4 ISOLATES THE NON-SCORE PART")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    # Group tournaments by score sequence
    by_score = defaultdict(list)
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))
        H2 = compute_H2(A, n)
        H4 = H - EH - H2
        by_score[scores].append((bits, H, H2, H4))

    print(f"\nn={n}: Score sequence analysis")
    print(f"  E[H] = {EH}")
    print(f"  {'Score seq':>22s} | count | H range       | H_4 range")
    print(f"  {'':->22s}-+-------+---------------+----------")

    for scores in sorted(by_score.keys()):
        entries = by_score[scores]
        Hs = [e[1] for e in entries]
        H4s = [round(e[3], 4) for e in entries]

        H_range = f"{min(Hs)}-{max(Hs)}" if min(Hs) != max(Hs) else str(min(Hs))
        H4_range = f"{min(H4s)}-{max(H4s)}" if min(H4s) != max(H4s) else str(min(H4s))

        H_var = len(set(Hs)) > 1

        print(f"  {str(scores):>22s} | {len(entries):>5d} | {H_range:>13s} | {H4_range:>10s}"
              + (" <-- H varies!" if H_var else ""))


# ========================================================================
# ANALYSIS 2: WHAT DETERMINES H_4 WITHIN A SCORE CLASS?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: WHAT DETERMINES H_4 WITHIN A SCORE CLASS?")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m
EH = math.factorial(n) / (2 ** (n - 1))

# For score classes where H varies, analyze what distinguishes the tournaments
# Focus on score classes with multiple distinct H values

by_score = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    H2 = compute_H2(A, n)
    H4 = H - EH - H2

    # Count directed 3-cycles
    c3_directed = 0
    c3_sets = set()
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            c3_directed += 1
            c3_sets.add(frozenset([a, b, c]))
        if A[a][c] and A[c][b] and A[b][a]:
            c3_directed += 1
            c3_sets.add(frozenset([a, b, c]))

    # Count disjoint 3-cycle pairs
    c3_list = list(c3_sets)
    disj_pairs = 0
    for i in range(len(c3_list)):
        for j in range(i+1, len(c3_list)):
            if not (c3_list[i] & c3_list[j]):
                disj_pairs += 1

    by_score[scores].append({
        'bits': bits, 'H': H, 'H2': round(H2, 4), 'H4': round(H4, 4),
        'c3_directed': c3_directed, 'c3_sets': len(c3_sets),
        'disj_pairs': disj_pairs
    })

# Analyze score classes where H varies
varying_scores = [s for s, entries in by_score.items()
                  if len(set(e['H'] for e in entries)) > 1]

print(f"\nn={n}: Score classes where H varies: {len(varying_scores)}")

for scores in sorted(varying_scores)[:5]:
    entries = by_score[scores]
    print(f"\n  Score {scores}:")
    print(f"    {'H':>5s} {'H_4':>8s} {'c3_dir':>7s} {'c3_sets':>8s} {'disj':>5s} | count")
    print(f"    {'':->5s}-{'':->8s}-{'':->7s}-{'':->8s}-{'':->5s}-+------")

    # Group by H
    by_H = defaultdict(list)
    for e in entries:
        by_H[e['H']].append(e)

    for H in sorted(by_H.keys()):
        group = by_H[H]
        c3d = set(e['c3_directed'] for e in group)
        c3s = set(e['c3_sets'] for e in group)
        dp = set(e['disj_pairs'] for e in group)
        h4 = set(e['H4'] for e in group)

        print(f"    {H:>5d} {str(sorted(h4)):>10s} {str(sorted(c3d)):>10s} {str(sorted(c3s)):>10s} {str(sorted(dp)):>8s} | {len(group):>5d}")


# ========================================================================
# ANALYSIS 3: H_4 vs DISJOINT CYCLE PAIRS
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: CORRELATION OF H_4 WITH DISJOINT PAIRS")
print("=" * 70)

# For each tournament, compute H_4 and count disjoint 3-cycle pairs
n = 6
m = n * (n - 1) // 2
total = 1 << m
EH = math.factorial(n) / (2 ** (n - 1))

h4_vs_disj = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    H2 = compute_H2(A, n)
    H4 = H - EH - H2

    c3_sets = set()
    for a, b, c in combinations(range(n), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.add(frozenset([a, b, c]))

    c3_list = list(c3_sets)
    disj = 0
    for i in range(len(c3_list)):
        for j in range(i+1, len(c3_list)):
            if not (c3_list[i] & c3_list[j]):
                disj += 1

    # Also count total directed odd cycles (for Omega)
    dir_cycles = 0
    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            dir_cycles += count_directed_cycles_on_subset(A, list(subset))

    h4_vs_disj.append((H4, disj, dir_cycles, H))

# Correlation
H4s = [x[0] for x in h4_vs_disj]
disjs = [x[1] for x in h4_vs_disj]
dcycs = [x[2] for x in h4_vs_disj]
Hs = [x[3] for x in h4_vs_disj]

def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
    sx = (sum((x - mx)**2 for x in xs) / n) ** 0.5
    sy = (sum((y - my)**2 for y in ys) / n) ** 0.5
    return cov / (sx * sy) if sx > 0 and sy > 0 else 0

print(f"\nn=6 correlations:")
print(f"  corr(H_4, disj_3_pairs) = {pearson(H4s, disjs):.6f}")
print(f"  corr(H_4, total_dir_cycles) = {pearson(H4s, dcycs):.6f}")
print(f"  corr(H, disj_3_pairs) = {pearson(Hs, disjs):.6f}")
print(f"  corr(H, total_dir_cycles) = {pearson(Hs, dcycs):.6f}")

# Group by (H_4, disj) to see if they're related
h4_disj_groups = defaultdict(int)
for h4, d, _, _ in h4_vs_disj:
    h4_disj_groups[(round(h4, 4), d)] += 1

print(f"\n  H_4 vs disj_3_pairs joint distribution:")
print(f"  {'H_4':>8s} | {'disj':>5s} | {'count':>6s}")
print(f"  {'':->8s}-+-{'':->5s}-+-{'':->6s}")
for (h4, d), cnt in sorted(h4_disj_groups.items()):
    print(f"  {h4:>8.2f} | {d:>5d} | {cnt:>6d}")


# ========================================================================
# ANALYSIS 4: WITHIN-SCORE VARIATION AT n=7 (sampling)
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: WITHIN-SCORE VARIATION AT n=7 (regular tournaments)")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
EH = math.factorial(n) / (2 ** (n - 1))

# At n=7, regular tournaments have scores (3,3,3,3,3,3,3)
# H_2 is the same for all of them
# H_4 distinguishes them
# We know from earlier: exactly 3 H values for regular n=7 (189, 175, 171)

import random
random.seed(42)

# Sample regular tournaments
regular_data = []
sample_attempts = 0
while len(regular_data) < 5000 and sample_attempts < 200000:
    # Generate random tournament
    bits = random.randint(0, (1 << m) - 1)
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if all(s == 3 for s in scores):
        H = count_ham_paths(A, n)
        H2 = compute_H2(A, n)
        H4 = H - EH - H2

        # Count c3 and c5 directed
        c3_dir = 0
        c3_sets = set()
        for a, b, c in combinations(range(n), 3):
            if A[a][b] and A[b][c] and A[c][a]:
                c3_dir += 1
                c3_sets.add(frozenset([a, b, c]))
            if A[a][c] and A[c][b] and A[b][a]:
                c3_dir += 1
                c3_sets.add(frozenset([a, b, c]))

        # c5 directed (via 5-vertex subsets)
        c5_dir = 0
        for subset in combinations(range(n), 5):
            c5_dir += count_directed_cycles_on_subset(A, list(subset))

        c3_list = list(c3_sets)
        disj_33 = 0
        for i in range(len(c3_list)):
            for j in range(i+1, len(c3_list)):
                if not (c3_list[i] & c3_list[j]):
                    disj_33 += 1

        regular_data.append({
            'H': H, 'H4': round(H4, 4),
            'c3_dir': c3_dir, 'c3_sets': len(c3_sets),
            'c5_dir': c5_dir, 'disj_33': disj_33
        })
    sample_attempts += 1

print(f"\nn=7: Found {len(regular_data)} regular tournaments from {sample_attempts} samples")
print(f"  E[H] = {EH:.2f}")

# Group by H
by_H = defaultdict(list)
for d in regular_data:
    by_H[d['H']].append(d)

print(f"\n  H    | H_4   | c3_dir | c3_sets | c5_dir | disj_33 | count")
print(f"  -----+-------+--------+---------+--------+---------+------")
for H in sorted(by_H.keys()):
    group = by_H[H]
    c3d = set(d['c3_dir'] for d in group)
    c3s = set(d['c3_sets'] for d in group)
    c5d = set(d['c5_dir'] for d in group)
    dp = set(d['disj_33'] for d in group)
    h4 = set(d['H4'] for d in group)
    print(f"  {H:>5d} | {str(sorted(h4)):>10s} | {str(sorted(c3d)):>6s} | {str(sorted(c3s)):>6s} | {str(sorted(c5d)):>8s} | {str(sorted(dp)):>6s} | {len(group)}")


# ========================================================================
# ANALYSIS 5: THE VITALI DIMENSION FORMULA
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: THE VITALI DIMENSION FORMULA")
print("=" * 70)

print("""
The "Vitali dimension" of H is the effective dimensionality of H as a
function on tournament space:

  Total dimension: m = C(n,2)
  Odd-degree (killed by symmetry): sum C(m,k) for k odd
  Even-degree non-zero: sum C(m,k) for k even, up to 2*floor((n-1)/2)

But the EFFECTIVE dimension is much smaller because:
1. Score regularity (captured by H_2) is essentially 1 DOF
2. Cycle structure (captured by H_4) adds a few more DOF
3. Higher terms contribute negligible energy

The "information dimension" of H:
  d_info = -sum p_i log p_i / log(max_i)
where p_i = energy(degree i) / total energy
""")

for n in range(3, 7):
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    H_vals = {}
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        sigma = tournament_to_sigma(A, n)
        H = count_ham_paths(A, n)
        H_vals[sigma] = H

    # Compute total energy by degree
    energy_by_deg = defaultdict(float)
    total_nonzero = 0
    for k in range(0, m + 1, 2):
        for S in combinations(range(m), k):
            coeff = 0
            for sigma, H in H_vals.items():
                chi = 1
                for idx in S:
                    chi *= sigma[idx]
                coeff += H * chi
            coeff /= total
            if abs(coeff) > 1e-10:
                energy_by_deg[k] += coeff ** 2
                total_nonzero += 1

    total_energy = sum(energy_by_deg.values())

    # Shannon entropy of energy distribution
    import math as m_mod
    entropy = 0
    for deg, e in energy_by_deg.items():
        p = e / total_energy
        if p > 0:
            entropy -= p * m_mod.log2(p)

    print(f"\n  n={n}:")
    print(f"    Total dimension: {m}")
    print(f"    Non-zero Fourier terms: {total_nonzero}")
    print(f"    Distinct H values: {len(set(H_vals.values()))}")
    print(f"    Energy entropy: {entropy:.4f} bits")
    print(f"    Effective dimension: ~{2**entropy:.1f} degrees of freedom")

    print(f"    Degree-energy: ", end="")
    for deg in sorted(energy_by_deg):
        pct = 100 * energy_by_deg[deg] / total_energy
        print(f"d{deg}={pct:.1f}% ", end="")
    print()


print("\n" + "=" * 70)
print("DONE.")
