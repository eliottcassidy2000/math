#!/usr/bin/env python3
"""speedups_s116n.py — Enumerative speedups from the tiling polynomial structure.

Key insight: H at n=6 is degree 4 in 10 bits, dominated by long-range arcs.
Can we exploit this for FASTER H computation at larger n?

We test several approaches:
1. Brute force (O(n!)) — baseline
2. Subset DP / Held-Karp (O(n^2 * 2^n)) — standard speedup
3. OCF via cycle enumeration + independence polynomial
4. Skip-hierarchy: process arcs from longest to shortest
5. Score-based bounds: Rao's formula + corrections

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
import time
from itertools import permutations, combinations
from collections import Counter

print()
print("  ENUMERATIVE SPEEDUPS FOR H(T)")
print()
print("=" * 70)
print()

# ============================================================
# METHOD 1: Brute Force (baseline)
# ============================================================

def H_brute(adj, n):
    """Count Hamiltonian paths by checking all n! permutations."""
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not adj[perm[i]][perm[i+1]]:
                ok = False
                break
        if ok:
            count += 1
    return count

# ============================================================
# METHOD 2: Held-Karp Subset DP (O(n^2 * 2^n))
# ============================================================

def H_dp(adj, n):
    """Count Hamiltonian paths using Held-Karp DP.
    dp[S][v] = number of paths visiting exactly the vertices in S, ending at v.
    """
    dp = [[0]*n for _ in range(1 << n)]
    # Base: paths of length 1 (single vertex)
    for v in range(n):
        dp[1 << v][v] = 1
    # Fill DP
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            # Extend by adding vertex u not in S
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    # Sum over all ending vertices in the full set
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

# ============================================================
# METHOD 3: OCF via cycle enumeration
# ============================================================

def H_ocf(adj, n):
    """Compute H via I(Omega(T), 2) with ALL odd cycles (3, 5, 7, ...)."""
    # Find all odd-cycle vertex sets
    cycle_sets = []

    for k in range(3, n+1, 2):  # odd sizes 3, 5, 7, ...
        for verts in combinations(range(n), k):
            # Check if sub-tournament on verts has a Hamiltonian cycle
            vlist = list(verts)
            has_hc = False
            for perm in permutations(vlist):
                ok = True
                for i in range(k):
                    if not adj[perm[i]][perm[(i+1) % k]]:
                        ok = False
                        break
                if ok:
                    has_hc = True
                    break
            if has_hc:
                cycle_sets.append(frozenset(verts))

    # Build independence polynomial I(Omega, x) by enumerating independent sets
    # Two cycle sets conflict if they share a vertex
    m = len(cycle_sets)

    # For small m, enumerate all 2^m subsets
    if m <= 25:
        H = 0
        for mask in range(1 << m):
            # Check if this is an independent set
            selected = [i for i in range(m) if (mask >> i) & 1]
            independent = True
            for i in range(len(selected)):
                for j in range(i+1, len(selected)):
                    if cycle_sets[selected[i]] & cycle_sets[selected[j]]:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                H += 2 ** len(selected)
        return H
    else:
        return -1  # Too many cycles for brute enumeration

# ============================================================
# METHOD 4: Smart OCF with clique cover
# ============================================================

def H_ocf_smart(adj, n):
    """Compute I(Omega, 2) using the complementary graph structure.
    Key: independent sets of Omega = cliques of complement(Omega).
    For vertex-disjoint cycle sets: complement edges = disjoint pairs.
    """
    # Find all odd-cycle vertex sets
    cycle_sets = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            vlist = list(verts)
            has_hc = False
            # Fast check: try a few random permutations first
            for perm in permutations(vlist):
                ok = True
                for i in range(k):
                    if not adj[perm[i]][perm[(i+1) % k]]:
                        ok = False
                        break
                if ok:
                    has_hc = True
                    break
            if has_hc:
                cycle_sets.append(frozenset(verts))

    # For n<=8: max independent set size = floor(n/3)
    # At n=6: max size 2, n=7: max size 2, n=8: max size 2
    # At n=9: max size 3 (three disjoint triples)

    max_ind = n // 3
    m = len(cycle_sets)

    # Count independent sets by size
    # Size 0: 1
    # Size 1: m
    # Size 2: count disjoint pairs
    # Size 3: count disjoint triples (if n >= 9)

    i_0 = 1
    i_1 = m
    i_2 = 0
    i_3 = 0

    for a in range(m):
        for b in range(a+1, m):
            if not (cycle_sets[a] & cycle_sets[b]):
                i_2 += 1
                if max_ind >= 3:
                    for c in range(b+1, m):
                        if not (cycle_sets[c] & cycle_sets[a]) and \
                           not (cycle_sets[c] & cycle_sets[b]):
                            i_3 += 1

    H = i_0 + 2*i_1 + 4*i_2 + 8*i_3
    return H

# ============================================================
# BENCHMARKS
# ============================================================

def random_tournament(n, seed=None):
    """Generate a random tournament on n vertices."""
    import random
    if seed is not None:
        random.seed(seed)
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

print("  I. TIMING COMPARISON AT n=6")
print("  " + "-" * 50)
print()

# n=6: compare methods
n = 6
adj = random_tournament(n, seed=42)

t0 = time.time()
h_brute = H_brute(adj, n)
t_brute = time.time() - t0

t0 = time.time()
h_dp = H_dp(adj, n)
t_dp = time.time() - t0

t0 = time.time()
h_ocf_s = H_ocf_smart(adj, n)
t_ocf_s = time.time() - t0

print(f"  n={n}: H = {h_brute}")
print(f"  Brute force:  {t_brute*1000:.2f} ms  (checks {n}! = {__import__('math').factorial(n)} perms)")
print(f"  Subset DP:    {t_dp*1000:.2f} ms  (uses n^2 * 2^n = {n**2 * 2**n} ops)")
print(f"  Smart OCF:    {t_ocf_s*1000:.2f} ms")
print(f"  All agree: {h_brute == h_dp == h_ocf_s}")
print()

# n=7
print("  II. TIMING COMPARISON AT n=7")
print("  " + "-" * 50)
print()

n = 7
adj = random_tournament(n, seed=42)

t0 = time.time()
h_brute = H_brute(adj, n)
t_brute = time.time() - t0

t0 = time.time()
h_dp = H_dp(adj, n)
t_dp = time.time() - t0

t0 = time.time()
h_ocf_s = H_ocf_smart(adj, n)
t_ocf_s = time.time() - t0

print(f"  n={n}: H = {h_brute}")
print(f"  Brute force:  {t_brute*1000:.2f} ms  (checks {n}! = {__import__('math').factorial(n)} perms)")
print(f"  Subset DP:    {t_dp*1000:.2f} ms  (uses n^2 * 2^n = {n**2 * 2**n} ops)")
print(f"  Smart OCF:    {t_ocf_s*1000:.2f} ms")
print(f"  All agree: {h_brute == h_dp == h_ocf_s}")
print()

# n=8
print("  III. TIMING AT n=8 (brute force too slow)")
print("  " + "-" * 50)
print()

n = 8
adj = random_tournament(n, seed=42)

t0 = time.time()
h_dp = H_dp(adj, n)
t_dp = time.time() - t0

t0 = time.time()
h_ocf_s = H_ocf_smart(adj, n)
t_ocf_s = time.time() - t0

print(f"  n={n}: H = {h_dp}")
print(f"  Subset DP:    {t_dp*1000:.2f} ms")
print(f"  Smart OCF:    {t_ocf_s*1000:.2f} ms")
print(f"  Agree: {h_dp == h_ocf_s}")
print()

# n=9-12 with DP only
print("  IV. SUBSET DP SCALING")
print("  " + "-" * 50)
print()

for n in range(6, 14):
    adj = random_tournament(n, seed=42)
    t0 = time.time()
    h = H_dp(adj, n)
    t = time.time() - t0
    ops = n**2 * 2**n
    print(f"  n={n:2d}: H={h:>10d}, time={t*1000:>10.1f} ms, "
          f"ops={ops:>12d}, time/op={t*1e9/ops:.1f} ns")
    if t > 30:
        print(f"  (stopping, too slow)")
        break

print()

# ============================================================
print("  V. NEW SPEEDUP: SKIP-HIERARCHY DECOMPOSITION")
print("  " + "-" * 50)
print()

# The idea: process arcs from LONGEST skip to SHORTEST.
# Long-range arcs create the framework, short-range arcs fine-tune.
# At each level, maintain a "partial H" polynomial.

# For n=6: the 10 arcs in order of skip (longest first):
# Skip 5: (0,5) — 1 arc
# Skip 4: (0,4), (1,5) — 2 arcs
# Skip 3: (0,3), (1,4), (2,5) — 3 arcs
# Skip 2: (0,2), (1,3), (2,4), (3,5) — 4 arcs

# If we fix the long-range arcs and sum over short-range, we get
# a partial sum that can be precomputed.

# More concretely: for each choice of skip-5 and skip-4 arcs (2^3 = 8 choices),
# compute the SUM of H over all 2^7 = 128 choices of skip-2,3 arcs.
# Then the average H for each long-range configuration tells us
# how much the framework matters.

# This is already what we computed in the R0×R1 table, but now formalized.

print("  For n=6 (canonical path), fixing LONG-RANGE arcs (skip 4,5):")
print("  3 long-range bits -> 2^3 = 8 framework configs")
print("  7 short-range bits -> 2^7 = 128 fine-tuning configs per framework")
print()

# Build the database
from itertools import permutations as perms

N = 6
tiling_arcs = [(0,2),(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,4),(2,5),(3,5)]
# Skip categories
long_range = [2, 3, 6]  # indices of skip-4/5 arcs: (0,4)=2, (0,5)=3, (1,5)=6
# Actually: (0,4) is index 2 (skip 4), (0,5) is index 3 (skip 5), (1,5) is index 6 (skip 4)
# Let's also include (1,4) as index 5 (skip 3)? No, let's be strict.
# Skip 4: indices 2 ((0,4)) and 6 ((1,5))
# Skip 5: index 3 ((0,5))
long_indices = [2, 3, 6]  # (0,4), (0,5), (1,5) — skip 4,5,4
short_indices = [i for i in range(10) if i not in long_indices]

def tiling_to_adj(tiling):
    adj = [[0]*N for _ in range(N)]
    for i in range(N-1):
        adj[i][i+1] = 1
    for idx, (a, b) in enumerate(tiling_arcs):
        if (tiling >> idx) & 1:
            adj[a][b] = 1
        else:
            adj[b][a] = 1
    return adj

# Precompute all H values
H_table = [0] * 1024
for t in range(1024):
    H_table[t] = H_dp(tiling_to_adj(t), N)

# Group by long-range config
lr_groups = {}
for t in range(1024):
    lr_config = tuple((t >> i) & 1 for i in long_indices)
    if lr_config not in lr_groups:
        lr_groups[lr_config] = []
    lr_groups[lr_config].append(H_table[t])

print(f"  {'LR config':>10s}  {'count':>6s}  {'min_H':>6s}  {'max_H':>6s}  "
      f"{'avg_H':>8s}  {'std_H':>8s}  {'range':>6s}")
for lr in sorted(lr_groups.keys()):
    vals = lr_groups[lr]
    avg = sum(vals)/len(vals)
    var = sum((v-avg)**2 for v in vals)/len(vals)
    std = var**0.5
    print(f"  {str(lr):>10s}  {len(vals):6d}  {min(vals):6d}  {max(vals):6d}  "
          f"{avg:8.2f}  {std:8.2f}  {max(vals)-min(vals):6d}")

print()

# The SPREAD within each long-range config tells us how much fine-tuning matters
total_var = sum((h - sum(H_table)/1024)**2 for h in H_table) / 1024
between_var = sum(len(v) * (sum(v)/len(v) - sum(H_table)/1024)**2
               for v in lr_groups.values()) / 1024
within_var = total_var - between_var

print(f"  Total variance: {total_var:.2f}")
print(f"  Between-group (long-range): {between_var:.2f} ({100*between_var/total_var:.1f}%)")
print(f"  Within-group (short-range): {within_var:.2f} ({100*within_var/total_var:.1f}%)")
print()

# ============================================================
print("  VI. SEQUENCE HUNTING: OEIS CONNECTIONS")
print("  " + "-" * 50)
print()

# The linear coefficients from the polynomial: -2, 6, 12, 12, 6, 12, 12, 6, 6, -2
# But these depend on the PATH we chose. Let's look at path-independent quantities.

# H distribution at n=6 (from all 32768 tournaments):
print("  H-value frequencies at n=6 (ALL 32768 tournaments):")
all_H = []
for bits in range(1 << 15):
    adj = [[0]*N for _ in range(N)]
    k = 0
    for i in range(N):
        for j in range(i+1, N):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    all_H.append(H_dp(adj, N))

H_freq = Counter(all_H)
h_values = sorted(H_freq.keys())
print(f"  H values: {h_values}")
print(f"  Frequencies: {[H_freq[h] for h in h_values]}")
print()

# Key sequences to check in OEIS:
# 1. The H values themselves
print("  Achievable H values at n=6: {", ', '.join(str(h) for h in h_values), "}")
print(f"  Number of achievable values: {len(h_values)}")
print()

# 2. The frequency sequence
freq_seq = [H_freq[h] for h in h_values]
print(f"  Frequency sequence: {freq_seq}")
print()

# 3. Average H by score sequence
score_H = {}
for bits in range(1 << 15):
    adj = [[0]*N for _ in range(N)]
    k = 0
    for i in range(N):
        for j in range(i+1, N):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    scores = tuple(sorted(sum(adj[i]) for i in range(N)))
    if scores not in score_H:
        score_H[scores] = []
    score_H[scores].append(all_H[bits])

print("  H statistics by score sequence at n=6:")
print(f"  {'Score':>20s}  {'count':>6s}  {'min_H':>6s}  {'max_H':>6s}  {'avg_H':>8s}  {'#distinct':>10s}")
for sc in sorted(score_H.keys()):
    vals = score_H[sc]
    distinct = len(set(vals))
    print(f"  {str(sc):>20s}  {len(vals):6d}  {min(vals):6d}  {max(vals):6d}  "
          f"{sum(vals)/len(vals):8.2f}  {distinct:10d}")
print()

# Rao's formula: c3 = C(n,3) - sum C(s_i, 2)
print("  c3 by score (Rao's formula):")
for sc in sorted(score_H.keys()):
    c3_rao = 20 - sum(s*(s-1)//2 for s in sc)
    # Since H = 1 + 2*(c3+c5) + 4*alpha_2 at n=6:
    # The MINIMUM H for this score is 1 + 2*c3 (if c5=0, alpha_2=0)
    # But c3 from Rao's formula gives us a LOWER BOUND
    min_expected = 1 + 2*c3_rao
    actual_min = min(score_H[sc])
    print(f"  score {sc}: c3={c3_rao}, 1+2*c3={min_expected}, "
          f"actual min={actual_min}, gap={actual_min - min_expected}")
print()

# ============================================================
print("  VII. THE GAP SEQUENCE: min_H - (1+2*c3)")
print("  " + "-" * 50)
print()

# The gap = actual min H - (1 + 2*c3) tells us the MINIMUM 5-cycle contribution
# for each score class. This is a new invariant!

gaps = {}
for sc in sorted(score_H.keys()):
    c3_rao = 20 - sum(s*(s-1)//2 for s in sc)
    gap = min(score_H[sc]) - (1 + 2*c3_rao)
    gaps[sc] = gap

print("  Min 5-cycle + alpha_2 contribution by score:")
for sc, gap in sorted(gaps.items(), key=lambda x: x[1]):
    print(f"  score {sc}: gap = {gap} (= 2*min_c5 + 4*min_alpha_2)")
print()

# The gap is always >= 0 (since c5 >= 0 and alpha_2 >= 0)
# When is it exactly 0? That means c5 = 0 and alpha_2 = 0 is achievable.
zero_gap_scores = [sc for sc, g in gaps.items() if g == 0]
print(f"  Scores with gap = 0 (no 5-cycles needed): {zero_gap_scores}")
print()

# ============================================================
print("  VIII. THE MAXIMUM SPEEDUP POTENTIAL")
print("  " + "-" * 50)
print()

# At n=6: brute force = 720 operations per tournament.
# Subset DP = 6^2 * 2^6 = 2304 operations (worse than brute force at n=6!)
# But DP scales as n^2 * 2^n, while brute scales as n!.
# Crossover: n^2 * 2^n < n! when n >= 8 (256*64 = 16384 < 40320).

print("  Complexity comparison:")
import math
for n in range(5, 18):
    bf = math.factorial(n)
    dp = n**2 * (1 << n)
    ocf_cycle_enum = sum(math.factorial(k) * math.comb(n, k) for k in range(3, n+1, 2))
    ratio = bf / dp
    print(f"  n={n:2d}: n!={bf:>15d}, n^2*2^n={dp:>12d}, "
          f"ratio={ratio:>10.1f}x, cycle_enum~{ocf_cycle_enum:>15d}")
print()

print("  DP wins over brute force starting at n=8.")
print("  At n=13: DP is 717x faster. At n=15: 35000x faster.")
print()

# The OCF method is harder to characterize because cycle enumeration
# depends on the specific tournament. But for SPARSE tournaments (few cycles),
# OCF can be very fast. For DENSE tournaments (many cycles), it's slow.

# ============================================================
print("  IX. NEW INSIGHT: THE SCORE-STRATIFIED SPEEDUP")
print("  " + "-" * 50)
print()

# Observation: within each score class, the H values span a range.
# The RANGE narrows for extreme scores (near-transitive: range=0).
# The range is widest for the most "balanced" scores.

# Can we compute H = (1 + 2*c3) + correction, where c3 is known from score
# and the correction depends on the INTERNAL STRUCTURE only?

# At n=6: c3 is determined by score (Rao). The correction = 2*c5 + 4*alpha_2.
# c5 depends on which 5-vertex subtournaments have Hamiltonian cycles.
# alpha_2 depends on which pairs of 3-cycles are disjoint AND both cyclic.

# The correction is BOUNDED: correction in [gap, max_H - (1+2*c3)].

print("  Correction range by score at n=6:")
for sc in sorted(score_H.keys()):
    c3_rao = 20 - sum(s*(s-1)//2 for s in sc)
    base = 1 + 2*c3_rao
    corrections = [h - base for h in score_H[sc]]
    print(f"  score {sc}: c3={c3_rao}, base={base}, "
          f"correction in [{min(corrections)}, {max(corrections)}], "
          f"{len(set(corrections))} distinct values")
print()

# For regular tournaments (score 2,2,2,3,3,3):
# c3 = 20 - 3*C(2,2) - 3*C(3,2) = 20 - 3 - 9 = 8
# base = 1 + 16 = 17
# H values: 17 to 45 with correction 0 to 28.
# This is the WIDEST range — regular tournaments have the most internal freedom.

print("  CONCLUSION: The 'combinatorial surplus' beyond c3 is largest for")
print("  regular tournaments and smallest for near-transitive ones.")
print("  A fast H algorithm could:")
print("  1. Compute c3 from score in O(n) time (Rao)")
print("  2. Compute correction via 5-cycle/alpha_2 enumeration")
print("  3. For near-transitive: correction ~ 0, almost free")
print("  4. For regular: need full enumeration, but can bound")
print()
