#!/usr/bin/env python3
"""
hidden_invariant_search.py -- kind-pasteur-2026-03-13-S61

Search for the hidden invariant that distinguishes tournaments with
the same lambda histogram and c5 but different H.

The specific case: score (0,2,3,3,4,4,5) at n=7, lambda=specific pattern,
c5=6 for both H=25 and H=29.

What distinguishes them? Candidates:
1. 5-cycle pair-coverage (lambda_5)
2. c7 count (Hamiltonian cycles)
3. The specific ORIENTATION of cycles (not just counts)
4. Higher-order overlap invariants (triples of 3-cycles)
5. Graph-theoretic properties (degree sequence of Omega)
6. The arc structure connecting cycle vertex sets

This is the DEEPEST hidden structure — the part that even the pair-coverage
hierarchy cannot detect.

Author: kind-pasteur-2026-03-13-S61
"""

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


def count_ham_paths(A, n):
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


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def get_lambda_histogram(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))

    pair_cov = {}
    for u, v in combinations(range(n), 2):
        lam = sum(1 for cs in c3_sets if u in cs and v in cs)
        pair_cov[(u, v)] = lam

    return tuple(sorted(pair_cov.values())), c3_sets


# ========================================================================
# ANALYSIS 1: FIND THE AMBIGUOUS CASES AT n=7
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: FINDING ALL AMBIGUOUS TOURNAMENTS AT n=7")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

# Collect tournaments in the first ambiguous class
target_score = (0, 2, 3, 3, 4, 4, 5)
ambig_data = []

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != target_score:
        continue

    H = count_ham_paths(A, n)
    lam_hist, c3_sets = get_lambda_histogram(A, n)

    target_lam = (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3)
    if lam_hist != target_lam:
        continue

    # c5 count
    c5_dir = 0
    for subset in combinations(range(n), 5):
        c5_dir += count_directed_ham_cycles_on_subset(A, list(subset))

    # c7 count
    c7_dir = count_directed_ham_cycles_on_subset(A, list(range(n)))

    # alpha_1, alpha_2
    # Enumerate ALL directed cycles
    all_cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            nc = count_directed_ham_cycles_on_subset(A, list(subset))
            for _ in range(nc):
                all_cycles.append((frozenset(subset), k))

    nc_total = len(all_cycles)

    # alpha_2: count disjoint pairs
    alpha_2 = 0
    for i in range(nc_total):
        for j in range(i+1, nc_total):
            if not (all_cycles[i][0] & all_cycles[j][0]):
                alpha_2 += 1

    ambig_data.append({
        'bits': bits, 'H': H, 'c3': len(c3_sets),
        'c5': c5_dir, 'c7': c7_dir, 'nc': nc_total, 'alpha2': alpha_2,
        'lam': lam_hist
    })

print(f"Found {len(ambig_data)} tournaments with target score and lambda")

# Group by H
by_H = defaultdict(list)
for d in ambig_data:
    by_H[d['H']].append(d)

for H in sorted(by_H.keys()):
    group = by_H[H]
    print(f"\n  H={H}: {len(group)} tournaments")
    c5_set = set(d['c5'] for d in group)
    c7_set = set(d['c7'] for d in group)
    nc_set = set(d['nc'] for d in group)
    a2_set = set(d['alpha2'] for d in group)
    print(f"    c5_dir = {sorted(c5_set)}")
    print(f"    c7_dir = {sorted(c7_set)}")
    print(f"    total cycles = {sorted(nc_set)}")
    print(f"    alpha_2 = {sorted(a2_set)}")

    # Show first few tournaments
    for d in group[:3]:
        print(f"    bits={d['bits']}: c5={d['c5']}, c7={d['c7']}, "
              f"nc={d['nc']}, alpha2={d['alpha2']}")


# ========================================================================
# ANALYSIS 2: WHAT DISTINGUISHES H=25 FROM H=29?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: WHAT DISTINGUISHES H=25 FROM H=29?")
print("=" * 70)

if 25 in by_H and 29 in by_H:
    # Pick one from each
    t25 = by_H[25][0]
    t29 = by_H[29][0]

    A25 = binary_to_tournament(t25['bits'], n)
    A29 = binary_to_tournament(t29['bits'], n)

    print(f"\n  H=25 tournament (bits={t25['bits']}):")
    print(f"    Out-degrees: {[sum(A25[v]) for v in range(n)]}")
    print(f"    c5={t25['c5']}, c7={t25['c7']}, nc={t25['nc']}, alpha2={t25['alpha2']}")

    print(f"\n  H=29 tournament (bits={t29['bits']}):")
    print(f"    Out-degrees: {[sum(A29[v]) for v in range(n)]}")
    print(f"    c5={t29['c5']}, c7={t29['c7']}, nc={t29['nc']}, alpha2={t29['alpha2']}")

    # Check: is the difference in c7 (Hamiltonian cycles)?
    # From OCF: H = 1 + 2*alpha_1 + 4*alpha_2
    # If alpha_1 differs but alpha_2 is same (or vice versa), that explains H difference
    # H=29-25=4. This could be:
    #   2*delta_alpha1 + 4*delta_alpha2 = 4
    # If delta_alpha2=0: delta_alpha1=2 (two more cycles)
    # If delta_alpha2=1: delta_alpha1=0 (one more disjoint pair, same cycles)

    delta_nc = t29['nc'] - t25['nc']
    delta_a2 = t29['alpha2'] - t25['alpha2']
    print(f"\n  Delta: nc(29)-nc(25) = {delta_nc}, alpha2(29)-alpha2(25) = {delta_a2}")
    print(f"  Check: 2*delta_nc + 4*delta_a2 = {2*delta_nc + 4*delta_a2}")
    print(f"  Actual delta H = {29 - 25} = 4")

    # If delta = from alpha_1 (total cycles), then it's the 5-cycle or 7-cycle count
    # that differs despite same lambda and "same c5_dir"
    # Wait: c5 was reported as 6 for both. Let me recheck.

    # Actually, the lambda_completeness script reported c5=6 for both.
    # But if nc (total cycles) differs, then either c5 or c7 must differ.
    print(f"\n  c5: t25={t25['c5']}, t29={t29['c5']}")
    print(f"  c7: t25={t25['c7']}, t29={t29['c7']}")
    print(f"  c3: t25={t25['c3']}, t29={t29['c3']}")

    # If c5 and c3 are the same but nc differs, it must be c7
    if t25['c5'] == t29['c5'] and t25['c3'] == t29['c3']:
        print(f"  => c3 and c5 are SAME. Difference must be in c7!")
        print(f"     c7 difference = {t29['c7'] - t25['c7']}")
        print(f"     This gives alpha_1 difference = {t29['c7'] - t25['c7']}")
        print(f"     2*delta_alpha1 = {2*(t29['c7'] - t25['c7'])} = delta H? "
              f"{'YES' if 2*(t29['c7'] - t25['c7']) + 4*delta_a2 == 4 else 'NO'}")

    # The key insight: SAME lambda histogram, SAME c5_dir, but DIFFERENT c7_dir!
    # This means c7 is NOT determined by (lambda, c5) for non-regular tournaments.
    # The Hamiltonian cycle count contains independent information.


# ========================================================================
# ANALYSIS 3: ALL AMBIGUOUS CASES — WHAT'S THE COMMON PATTERN?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: SYSTEMATIC — WHAT DISTINGUISHES ALL AMBIGUOUS CASES?")
print("=" * 70)

# Check multiple ambiguous score classes
ambig_scores = [
    (0, 2, 3, 3, 4, 4, 5),
    (1, 1, 2, 3, 4, 5, 5),
    (1, 2, 2, 3, 3, 4, 6),
]

for target_sc in ambig_scores:
    print(f"\n  Score {target_sc}:")

    data = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))
        if scores != target_sc:
            continue

        H = count_ham_paths(A, n)
        lam_hist, c3_sets = get_lambda_histogram(A, n)

        c5_dir = 0
        for subset in combinations(range(n), 5):
            c5_dir += count_directed_ham_cycles_on_subset(A, list(subset))

        c7_dir = count_directed_ham_cycles_on_subset(A, list(range(n)))

        data.append({
            'H': H, 'lam': lam_hist, 'c3': len(c3_sets),
            'c5': c5_dir, 'c7': c7_dir
        })

    # Group by lambda
    by_lam = defaultdict(list)
    for d in data:
        by_lam[d['lam']].append(d)

    for lam in sorted(by_lam.keys()):
        group = by_lam[lam]
        H_set = sorted(set(d['H'] for d in group))
        if len(H_set) <= 1:
            continue

        c5_set = sorted(set(d['c5'] for d in group))
        c7_set = sorted(set(d['c7'] for d in group))

        print(f"    lambda=... ({len(lam)} values), H={H_set}")
        print(f"    c5={c5_set}, c7={c7_set}")

        # Group by (c5, c7)
        by_c5c7 = defaultdict(set)
        for d in group:
            by_c5c7[(d['c5'], d['c7'])].add(d['H'])

        still_ambig = False
        for (c5, c7), Hs in sorted(by_c5c7.items()):
            if len(Hs) > 1:
                still_ambig = True
                print(f"      (c5={c5}, c7={c7}) -> H = {sorted(Hs)} [STILL AMBIGUOUS]")
            else:
                print(f"      (c5={c5}, c7={c7}) -> H = {sorted(Hs)[0]}")

        if not still_ambig:
            print(f"    => (c5, c7) RESOLVES the ambiguity!")
        else:
            print(f"    => (c5, c7) does NOT resolve! Deeper invariant needed.")


# ========================================================================
# ANALYSIS 4: THE n=6 AMBIGUOUS CASES
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: n=6 AMBIGUOUS CASES — WHAT RESOLVES?")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

n6_ambig_scores = [
    (1, 2, 2, 3, 3, 4),
    (1, 2, 3, 3, 3, 3),
    (2, 2, 2, 2, 3, 4),
]

for target_sc in n6_ambig_scores:
    print(f"\n  Score {target_sc}:")

    data = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))
        if scores != target_sc:
            continue

        H = count_ham_paths(A, n)
        lam_hist, c3_sets = get_lambda_histogram(A, n)

        c5_dir = 0
        for subset in combinations(range(n), 5):
            c5_dir += count_directed_ham_cycles_on_subset(A, list(subset))

        data.append({
            'H': H, 'lam': lam_hist, 'c3': len(c3_sets), 'c5': c5_dir
        })

    by_lam = defaultdict(list)
    for d in data:
        by_lam[d['lam']].append(d)

    for lam in sorted(by_lam.keys()):
        group = by_lam[lam]
        H_set = sorted(set(d['H'] for d in group))
        if len(H_set) <= 1:
            continue

        c5_set = sorted(set(d['c5'] for d in group))
        print(f"    lambda=..., H={H_set}")
        print(f"    c5={c5_set}")

        # Check if c5 resolves
        by_c5 = defaultdict(set)
        for d in group:
            by_c5[d['c5']].add(d['H'])

        for c5, Hs in sorted(by_c5.items()):
            if len(Hs) > 1:
                print(f"      c5={c5} -> H = {sorted(Hs)} [STILL AMBIGUOUS]")
            else:
                print(f"      c5={c5} -> H = {sorted(Hs)[0]}")


# ========================================================================
# ANALYSIS 5: THE DEEPEST — TRIPLE OVERLAP INVARIANT
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: TRIPLE OVERLAP INVARIANT")
print("=" * 70)

# For cases where (lambda, c5, c7) doesn't resolve, try triple overlap:
# For each triple of 3-cycles (C_i, C_j, C_k), compute the triple overlap
# |V(C_i) cap V(C_j) cap V(C_k)|.

# This is a 3rd-order invariant of the cycle structure.

# Test on the n=7 ambiguous case (0,2,3,3,4,4,5)
n = 7
m = n * (n - 1) // 2
total = 1 << m

target_sc = (0, 2, 3, 3, 4, 4, 5)
target_lam = (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3)

print(f"\n  Score {target_sc}, target lambda:")

data = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != target_sc:
        continue

    lam_hist, c3_sets = get_lambda_histogram(A, n)
    if lam_hist != target_lam:
        continue

    H = count_ham_paths(A, n)
    nc3 = len(c3_sets)

    # Triple overlap histogram
    triple_hist = defaultdict(int)
    for i in range(nc3):
        for j in range(i+1, nc3):
            for k_idx in range(j+1, nc3):
                t_ov = len(c3_sets[i] & c3_sets[j] & c3_sets[k_idx])
                triple_hist[t_ov] += 1

    triple_key = tuple(sorted(triple_hist.items()))

    # Per-vertex cycle count
    per_vertex = tuple(sorted(sum(1 for cs in c3_sets if v in cs) for v in range(n)))

    # Cycle degree: for each 3-cycle, how many other 3-cycles share >= 1 vertex
    cycle_degrees = []
    for i in range(nc3):
        deg = sum(1 for j in range(nc3) if j != i and c3_sets[i] & c3_sets[j])
        cycle_degrees.append(deg)
    cycle_deg_hist = tuple(sorted(cycle_degrees))

    data.append({
        'bits': bits, 'H': H, 'triple': triple_key,
        'per_vertex': per_vertex, 'cycle_deg': cycle_deg_hist,
        'c3': nc3
    })

# Group by H
by_H = defaultdict(list)
for d in data:
    by_H[d['H']].append(d)

for H in sorted(by_H.keys()):
    group = by_H[H]
    triples = set(d['triple'] for d in group)
    per_verts = set(d['per_vertex'] for d in group)
    cycle_degs = set(d['cycle_deg'] for d in group)

    print(f"\n  H={H}: {len(group)} tournaments")
    print(f"    Triple overlap histograms: {len(triples)}")
    for t in sorted(triples):
        cnt = sum(1 for d in group if d['triple'] == t)
        print(f"      {dict(t)}: {cnt}")
    print(f"    Per-vertex cycle distributions: {sorted(per_verts)}")
    print(f"    Cycle degree distributions: {sorted(cycle_degs)}")

# Check: does triple overlap resolve the ambiguity?
by_triple = defaultdict(set)
for d in data:
    by_triple[d['triple']].add(d['H'])

triple_resolves = all(len(v) == 1 for v in by_triple.values())
print(f"\n  Triple overlap resolves ambiguity: {triple_resolves}")

# Check: does cycle degree resolve?
by_cdeg = defaultdict(set)
for d in data:
    by_cdeg[d['cycle_deg']].add(d['H'])

cdeg_resolves = all(len(v) == 1 for v in by_cdeg.values())
print(f"  Cycle degree histogram resolves ambiguity: {cdeg_resolves}")

# Check: does per-vertex cycle count resolve?
by_pv = defaultdict(set)
for d in data:
    by_pv[d['per_vertex']].add(d['H'])

pv_resolves = all(len(v) == 1 for v in by_pv.values())
print(f"  Per-vertex cycle count resolves ambiguity: {pv_resolves}")


print(f"\n\n{'='*70}")
print("DONE.")
print("=" * 70)
