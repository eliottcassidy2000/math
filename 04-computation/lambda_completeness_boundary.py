#!/usr/bin/env python3
"""
lambda_completeness_boundary.py -- kind-pasteur-2026-03-13-S61

Determine the exact boundary where lambda pair-coverage stops being
a complete invariant for H.

Key findings so far:
- n=7 regular: lambda DETERMINES H (3 classes)
- n=6 score (2,2,2,3,3,3): lambda DETERMINES H (4 classes)
- n=7 non-regular: lambda does NOT determine H

Question: for ALL score classes, does lambda determine H?
And if not, what additional invariant is needed?

Also explore the connection to Vitali sets more deeply:
- The "measurability" hierarchy: score -> lambda -> full tournament
- When lambda suffices, the Vitali partition collapses
- When it doesn't, the "non-measurable" c5 information matters

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
    """Compute the sorted pair-coverage histogram of 3-cycles."""
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

    return tuple(sorted(pair_cov.values())), len(c3_sets)


# ========================================================================
# ANALYSIS 1: EXHAUSTIVE n=5 — LAMBDA COMPLETENESS
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: EXHAUSTIVE n=5 — DOES LAMBDA DETERMINE H?")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m

by_score_lambda = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    lam_hist, nc3 = get_lambda_histogram(A, n)

    by_score_lambda[(scores, lam_hist)].append(H)

print(f"Total tournaments: {total}")
print(f"Distinct (score, lambda) classes: {len(by_score_lambda)}")

ambig_count = 0
for (sc, lam), Hs in sorted(by_score_lambda.items()):
    H_set = sorted(set(Hs))
    if len(H_set) > 1:
        ambig_count += 1
        print(f"  AMBIGUOUS: score={sc}, lambda={lam}: H = {H_set} ({len(Hs)} tournaments)")
    else:
        pass  # No ambiguity

if ambig_count == 0:
    print("  Lambda DETERMINES H within each score class at n=5!")
else:
    print(f"  {ambig_count} ambiguous (score, lambda) classes at n=5")

# Also check: does lambda alone (without score) determine H?
by_lambda_only = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    lam_hist, _ = get_lambda_histogram(A, n)
    by_lambda_only[lam_hist].append(H)

lambda_ambig = 0
for lam, Hs in sorted(by_lambda_only.items()):
    H_set = sorted(set(Hs))
    if len(H_set) > 1:
        lambda_ambig += 1
        print(f"  Lambda-only AMBIGUOUS: lambda={lam}: H = {H_set}")

if lambda_ambig == 0:
    print("  Lambda ALONE (without score) determines H at n=5!")
else:
    print(f"  {lambda_ambig} ambiguous lambda classes at n=5")


# ========================================================================
# ANALYSIS 2: EXHAUSTIVE n=6 — LAMBDA COMPLETENESS
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: EXHAUSTIVE n=6 — DOES LAMBDA DETERMINE H?")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

by_score_lambda = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    lam_hist, nc3 = get_lambda_histogram(A, n)

    by_score_lambda[(scores, lam_hist)].append(H)

print(f"Total tournaments: {total}")
print(f"Distinct (score, lambda) classes: {len(by_score_lambda)}")

ambig_count = 0
total_ambig_tours = 0
for (sc, lam), Hs in sorted(by_score_lambda.items()):
    H_set = sorted(set(Hs))
    if len(H_set) > 1:
        ambig_count += 1
        total_ambig_tours += len(Hs)
        print(f"  AMBIGUOUS: score={sc}, lambda={lam}: H = {H_set} ({len(Hs)} tours)")

if ambig_count == 0:
    print("  Lambda DETERMINES H within each score class at n=6!")
else:
    print(f"\n  {ambig_count} ambiguous (score, lambda) classes at n=6")
    print(f"  {total_ambig_tours}/{total} tournaments in ambiguous classes ({100*total_ambig_tours/total:.1f}%)")

# Lambda alone
by_lambda_only = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    lam_hist, _ = get_lambda_histogram(A, n)
    by_lambda_only[lam_hist].append(H)

lambda_ambig = 0
for lam, Hs in sorted(by_lambda_only.items()):
    H_set = sorted(set(Hs))
    if len(H_set) > 1:
        lambda_ambig += 1

if lambda_ambig == 0:
    print("  Lambda ALONE (without score) determines H at n=6!")
else:
    print(f"  {lambda_ambig} ambiguous lambda-only classes at n=6")


# ========================================================================
# ANALYSIS 3: n=7 — WHICH SCORE CLASSES HAVE LAMBDA COMPLETENESS?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: n=7 — LAMBDA COMPLETENESS BY SCORE CLASS")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

# Full enumeration of n=7 is 2^21 = 2M tournaments - feasible but slow
# Let's do it
by_score_lambda = defaultdict(list)
count = 0
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    lam_hist, nc3 = get_lambda_histogram(A, n)

    by_score_lambda[(scores, lam_hist)].append(H)
    count += 1
    if count % 500000 == 0:
        print(f"  Processed {count}/{total}...")

print(f"\nTotal tournaments: {total}")
print(f"Distinct (score, lambda) classes: {len(by_score_lambda)}")

# Count ambiguous classes
ambig_by_score = defaultdict(list)
total_ambig = 0
total_complete = 0
for (sc, lam), Hs in by_score_lambda.items():
    H_set = sorted(set(Hs))
    if len(H_set) > 1:
        ambig_by_score[sc].append((lam, H_set, len(Hs)))
        total_ambig += 1
    else:
        total_complete += 1

print(f"\nComplete (lambda determines H): {total_complete}")
print(f"Ambiguous (lambda does NOT determine H): {total_ambig}")

if total_ambig == 0:
    print("\nLAMBDA DETERMINES H FOR ALL SCORE CLASSES AT n=7!")
else:
    print(f"\nAmbiguous score classes:")
    for sc in sorted(ambig_by_score.keys()):
        cases = ambig_by_score[sc]
        print(f"\n  Score {sc}:")
        for lam, H_set, cnt in sorted(cases, key=lambda x: -len(x[1])):
            if len(H_set) <= 5:
                print(f"    lambda={lam}: H = {H_set} ({cnt} tours)")
            else:
                print(f"    lambda={lam}: H = {H_set[:5]}... ({len(H_set)} values, {cnt} tours)")

    # Summary statistics
    ambig_scores = set(ambig_by_score.keys())
    all_scores = set(sc for (sc, _) in by_score_lambda.keys())
    complete_scores = all_scores - ambig_scores
    print(f"\n  Score classes where lambda is complete: {len(complete_scores)}/{len(all_scores)}")
    print(f"  Score classes where lambda is incomplete: {len(ambig_scores)}/{len(all_scores)}")

    # Which score classes are complete?
    print(f"\n  Complete score classes:")
    for sc in sorted(complete_scores):
        H_vals = set()
        for (s, l), Hs in by_score_lambda.items():
            if s == sc:
                H_vals.update(Hs)
        print(f"    {sc}: H in {sorted(H_vals)}")


# ========================================================================
# ANALYSIS 4: EXTENDED LAMBDA — 5-CYCLE PAIR-COVERAGE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: EXTENDED LAMBDA WITH 5-CYCLE COVERAGE (n=6)")
print("=" * 70)

# If 3-cycle lambda doesn't determine H, maybe adding 5-cycle information helps.
# Define lambda_5(u,v) = number of 5-cycle vertex sets containing both u and v.

n = 6
m = n * (n - 1) // 2
total = 1 << m

# Check ambiguous cases from n=6
by_extended = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    # 3-cycle pair coverage
    lam3_hist, nc3 = get_lambda_histogram(A, n)

    # 5-cycle pair coverage
    c5_vsets = set()
    for subset in combinations(range(n), 5):
        nc5 = count_directed_ham_cycles_on_subset(A, list(subset))
        if nc5 > 0:
            c5_vsets.add(frozenset(subset))

    lam5 = {}
    for u, v in combinations(range(n), 2):
        lam5[(u, v)] = sum(1 for fs in c5_vsets if u in fs and v in fs)
    lam5_hist = tuple(sorted(lam5.values()))

    by_extended[(scores, lam3_hist, lam5_hist)].append(H)

ext_ambig = 0
for key, Hs in by_extended.items():
    H_set = sorted(set(Hs))
    if len(H_set) > 1:
        ext_ambig += 1
        sc, l3, l5 = key
        print(f"  STILL AMBIGUOUS: score={sc}, lam3={l3}, lam5={l5}: H={H_set}")

if ext_ambig == 0:
    print("  Extended lambda (3-cycle + 5-cycle coverage) DETERMINES H at n=6!")
else:
    print(f"  {ext_ambig} still ambiguous with extended lambda at n=6")


# ========================================================================
# ANALYSIS 5: THE HIERARCHY AT n=7 — WHAT ELSE IS NEEDED?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: WHAT ADDITIONAL INVARIANT IS NEEDED AT n=7?")
print("=" * 70)

# For ambiguous cases at n=7, check if adding c5 or c7 information resolves them
n = 7
m = n * (n - 1) // 2
total = 1 << m

if total_ambig > 0:
    # Pick one ambiguous score class and study it
    ambig_sc = sorted(ambig_by_score.keys())[0]
    print(f"\nStudying ambiguous score class {ambig_sc}:")

    # Collect tournaments in this class
    class_data = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))
        if scores != ambig_sc:
            continue

        H = count_ham_paths(A, n)
        lam_hist, nc3 = get_lambda_histogram(A, n)

        # Count c5 (directed)
        c5_dir = 0
        for subset in combinations(range(n), 5):
            c5_dir += count_directed_ham_cycles_on_subset(A, list(subset))

        class_data.append({
            'H': H, 'lam': lam_hist, 'c5': c5_dir, 'nc3': nc3
        })
        if len(class_data) >= 5000:
            break

    # Group by lambda and check if c5 resolves
    by_lam = defaultdict(list)
    for d in class_data:
        by_lam[d['lam']].append(d)

    for lam in sorted(by_lam.keys()):
        group = by_lam[lam]
        H_set = sorted(set(d['H'] for d in group))
        if len(H_set) > 1:
            # Within this lambda class, does c5 resolve?
            by_c5 = defaultdict(set)
            for d in group:
                by_c5[d['c5']].add(d['H'])

            c5_resolves = all(len(v) == 1 for v in by_c5.values())
            if c5_resolves:
                print(f"  lambda={lam}: H={H_set}, c5 RESOLVES: "
                      f"{sorted((c5, sorted(h)) for c5, h in by_c5.items())}")
            else:
                print(f"  lambda={lam}: H={H_set}, c5 does NOT resolve")
                for c5_val, h_vals in sorted(by_c5.items()):
                    if len(h_vals) > 1:
                        print(f"    c5={c5_val}: H={sorted(h_vals)}")

print(f"\n\n{'='*70}")
print("DONE.")
print("=" * 70)
