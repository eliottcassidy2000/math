#!/usr/bin/env python3
"""
Z_11 Circulant Tournament OCF Structure Analysis

For each of the 32 circulant tournaments on Z_11, compute:
- H (Hamiltonian path count) via Held-Karp DP
- All odd directed cycles (lengths 3,5,7,9,11)
- OCF decomposition: alpha_k via inclusion-exclusion on vertex subsets
  I(Omega,x) at x=2 computed via vertex-subset DP
- Verify H = I(Omega, 2)
- Spectral spread

Key insight: Instead of building the conflict graph explicitly (too many cycles),
compute I(Omega, 2) = sum over vertex subsets V' of (-1)^{n-|V'|} * ... NO.

Better: Use the cycle-space formulation directly.
H(T) = sum_{S subset of [n], S independent in Omega} 2^{|S|}

Alternative efficient approach:
For tournament on n=11 vertices, the conflict graph Omega has vertices = odd cycles.
Two cycles conflict if they share a tournament vertex.

I(Omega, x) = sum_k alpha_k x^k where alpha_k = # independent sets of size k.

EFFICIENT METHOD: Compute I(Omega, 2) via vertex (tournament vertex) DP.
For each tournament vertex subset U, count the number of collections of
vertex-disjoint odd cycles using exactly the vertices in U. Call this c(U).
Then I(Omega, 2) = sum_U c(U) * 2^{number of cycles in collection}... no,
I(Omega, 2) = sum over independent sets S in Omega, of 2^|S|.
Each independent set S is a set of pairwise vertex-disjoint odd cycles.
If S has k cycles covering vertex set U, contribution is 2^k.

So: I(Omega, 2) = sum over vertex subsets U of f(U)
where f(U) = sum over ways to partition U into vertex-disjoint odd directed cycles
             of 2^{number of cycles in the partition}.

This is a DP over subsets of [n]! With n=11 that's 2^11 = 2048 subsets.
For each subset U, we need to know which odd cycles use exactly U.

f(emptyset) = 1 (the empty collection, 2^0 = 1)
f(U) = sum over odd cycles C contained in U with min(C) = min(U):
        2 * f(U \ V(C))

The "min(C) = min(U)" constraint ensures each partition is counted once.
And the factor 2 accounts for the 2^k weighting (each cycle contributes factor 2).

Then I(Omega, 2) = sum_U f(U).
Wait, we want I(Omega,2) = sum over all collections of vertex-disjoint odd cycles,
weighted by 2^{size of collection}. And f(U) counts such collections using vertex set U.
So I(Omega,2) = sum_U f(U), yes.

But actually H = I(Omega,2) = sum over independent sets S, 2^|S|.
And f(U) = sum over independent sets S whose union of vertex sets is U, 2^|S|.
So I(Omega,2) = sum_U f(U). Correct.

For the DP:
- Pre-compute for each subset U of size 3,5,7,9,11: the number of directed
  odd cycles on exactly that vertex set. Call this numcyc(U).
- Then f(emptyset) = 1
- For nonempty U with smallest element v:
  f(U) = sum over subsets W of U with v in W, |W| odd, |W|>=3:
          numcyc(W) * 2 * f(U \ W)

This gives us H = sum_U f(U) AND the alpha_k decomposition if we track cycle count.

With n=11 and 2^11=2048 subsets, this is very fast!
"""

import itertools
import numpy as np
from collections import defaultdict, Counter
import time

n = 11

pairs = [(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)]

def get_connection_set(choices):
    S = set()
    for i, c in enumerate(choices):
        S.add(pairs[i][c])
    return frozenset(S)

def build_adjacency(S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            j = (i + s) % n
            A[i][j] = 1
    return A

def held_karp_count(A):
    """Count all Hamiltonian paths using Held-Karp DP."""
    full = (1 << n) - 1
    # Use array for speed
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask][v]
            if c == 0:
                continue
            for u in range(n):
                if (mask & (1 << u)) == 0 and A[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full][v] for v in range(n))

def count_directed_cycles_on_subset(A, subset_mask, verts):
    """Count the number of directed cycles using exactly the vertices in subset.
    verts: list of vertex indices corresponding to set bits in subset_mask.
    Returns: number of distinct directed cycles on these vertices.
    """
    k = len(verts)
    if k < 3 or k % 2 == 0:
        return 0

    # Fix start vertex as verts[0] (smallest), enumerate Hamiltonian cycles
    # on this subset returning to start
    start = verts[0]
    rest = verts[1:]
    full_local = (1 << k) - 1

    # Held-Karp for Hamiltonian cycles on subset
    # dp[mask][v_idx] = number of paths from start visiting mask, ending at v_idx
    # mask bit i corresponds to verts[i]
    dp_local = [[0]*k for _ in range(1 << k)]
    # Start at index 0
    dp_local[1][0] = 1

    for mask in range(1, 1 << k):
        for vi in range(k):
            if not (mask & (1 << vi)):
                continue
            c = dp_local[mask][vi]
            if c == 0:
                continue
            for ui in range(1, k):  # don't go back to 0 until the end
                if (mask & (1 << ui)) == 0 and A[verts[vi]][verts[ui]]:
                    dp_local[mask | (1 << ui)][ui] += c

    # Count cycles: paths visiting all vertices, with edge back to start
    count = 0
    for vi in range(1, k):
        if dp_local[full_local][vi] > 0 and A[verts[vi]][start]:
            count += dp_local[full_local][vi]

    return count

def spectral_spread(A):
    An = np.array(A)
    eigenvalues = np.linalg.eigvalsh(An + An.T)
    return (max(eigenvalues) - min(eigenvalues)) / 2

def compute_ocf_dp(A):
    """
    Compute OCF decomposition using subset DP.
    Returns: (H_value, alpha_dict) where alpha_dict[k] = number of independent
    sets of size k in the odd-cycle conflict graph.

    Method: DP over subsets of [n].
    g[mask] = polynomial in formal variable (tracking number of cycles)
    g[mask] = list where g[mask][k] = number of ways to partition the vertices
              in mask into k vertex-disjoint odd directed cycles.

    Then alpha_k = sum_mask g[mask][k], and H = sum_k alpha_k * 2^k.
    """
    # Step 1: Precompute numcyc[mask] for all odd-size subsets
    print("  Precomputing cycle counts per subset...", flush=True)
    t0 = time.time()
    numcyc = [0] * (1 << n)

    for size in [3, 5, 7, 9, 11]:
        count_at_size = 0
        for subset in itertools.combinations(range(n), size):
            mask = sum(1 << v for v in subset)
            c = count_directed_cycles_on_subset(A, mask, list(subset))
            numcyc[mask] = c
            count_at_size += c
        print(f"    {size}-cycles: {count_at_size}", flush=True)

    t1 = time.time()
    print(f"  Cycle precomputation: {t1 - t0:.1f}s", flush=True)

    # Step 2: DP over subsets
    # g[mask][k] = number of ways to cover exactly mask using k vertex-disjoint odd cycles
    # Transition: for nonempty mask, let v = lowest set bit.
    # g[mask][k] = sum over subsets W of mask containing v, |W| odd, |W|>=3:
    #              numcyc[W] * g[mask ^ W][k-1]
    # Base: g[0][0] = 1

    max_cycles = n // 3  # at most floor(11/3) = 3 cycles (with 3-cycles)
    # Actually up to 3 disjoint 3-cycles = 9 vertices, or 1 three + 1 five = 8, etc.
    # Max independent set size could be 3 (three 3-cycles covering 9 vertices)
    # or could have more with larger cycles? No, 3 disjoint cycles of size 3 = 9 verts.
    # 11 vertices: can't fit 4 disjoint 3-cycles. Could fit 3+3+5=11.
    # So max k = 3.
    max_k = 4  # generous bound

    # Use dict for sparse storage
    g = [None] * (1 << n)
    g[0] = [0] * (max_k + 1)
    g[0][0] = 1

    for mask in range(1, 1 << n):
        # Find lowest set bit
        v = 0
        while not (mask & (1 << v)):
            v += 1

        gm = [0] * (max_k + 1)

        # Enumerate all odd-sized subsets of mask containing v
        # Iterate over subsets of mask that include v
        # mask without v:
        rest = mask ^ (1 << v)
        # Iterate over subsets of rest, add v back
        sub = rest
        while True:
            W = sub | (1 << v)
            sz = bin(W).count('1')
            if sz >= 3 and sz % 2 == 1 and numcyc[W] > 0:
                complement = mask ^ W
                gc = g[complement]
                if gc is not None:
                    nc = numcyc[W]
                    for k in range(max_k):
                        if gc[k] > 0:
                            gm[k + 1] += nc * gc[k]
            if sub == 0:
                break
            sub = (sub - 1) & rest

        if any(x > 0 for x in gm):
            g[mask] = gm

    # Step 3: Aggregate
    alpha = defaultdict(int)
    alpha[0] = 1  # empty set of cycles

    for mask in range(1, 1 << n):
        if g[mask] is not None:
            for k in range(1, max_k + 1):
                alpha[k] += g[mask][k]

    H_check = sum(alpha[k] * (2 ** k) for k in alpha)
    return H_check, dict(alpha)


# ===== MAIN =====
print("=" * 80)
print("Z_11 CIRCULANT TOURNAMENT OCF STRUCTURE ANALYSIS")
print("=" * 80)
print(flush=True)

QR_11 = frozenset({1, 3, 4, 5, 9})
print(f"Paley QR_11 = {sorted(QR_11)}")
print(f"Pairs: {pairs}\n")

# Enumerate all 32 circulant tournaments
results = []
t_global = time.time()

print("--- Held-Karp DP for all 32 circulants ---", flush=True)
for bits in range(32):
    choices = tuple((bits >> i) & 1 for i in range(5))
    S = get_connection_set(choices)
    A = build_adjacency(S)
    H = held_karp_count(A)
    spread = spectral_spread(A)
    is_paley = (S == QR_11)
    results.append({
        'S': S, 'H': H, 'spread': spread, 'is_paley': is_paley, 'A': A,
    })

print(f"Done in {time.time() - t_global:.1f}s\n")

# Sort by H descending
results.sort(key=lambda r: -r['H'])

print(f"{'Rank':>4} {'S':>22} {'H':>10} {'spread':>8} {'Paley':>6}")
print("-" * 60)
for i, r in enumerate(results):
    tag = " ***" if r['is_paley'] else ""
    print(f"{i+1:>4} {str(sorted(r['S'])):>22} {r['H']:>10} {r['spread']:>8.3f}{tag}")

max_H = results[0]['H']
maximizers = [r for r in results if r['H'] == max_H]
print(f"\nH-maximizer: H = {max_H}")
for r in maximizers:
    print(f"  S = {sorted(r['S'])}, Paley = {r['is_paley']}")

paley_result = [r for r in results if r['is_paley']][0]
paley_rank = [i for i, r in enumerate(results) if r['is_paley']][0] + 1
print(f"Paley rank by H: {paley_rank}/32, H = {paley_result['H']}")

# Distinct H values
H_values = sorted(set(r['H'] for r in results), reverse=True)
print(f"\nDistinct H values: {len(H_values)}")
H_counts = Counter(r['H'] for r in results)
for h in H_values:
    print(f"  H = {h}: {H_counts[h]} tournaments")

# Complement check
print(f"\nComplement check (S vs 11-S):")
for r in results[:5]:
    S = r['S']
    S_comp = frozenset(n - s for s in S)
    H_comp = [r2['H'] for r2 in results if r2['S'] == S_comp]
    print(f"  S={sorted(S)}, comp={sorted(S_comp)}, H={r['H']}, H_comp={H_comp[0] if H_comp else '?'}")

# OCF decomposition for one representative per distinct H
print("\n" + "=" * 80)
print("OCF DECOMPOSITION (one representative per distinct H)")
print("=" * 80)

seen_H = set()
representatives = []
for r in results:
    if r['H'] not in seen_H:
        seen_H.add(r['H'])
        representatives.append(r)

# Always include Paley explicitly
paley_in = any(r['is_paley'] for r in representatives)

for r in representatives:
    tag = " [PALEY]" if r['is_paley'] else ""
    print(f"\n{'='*60}")
    print(f"S = {sorted(r['S'])}{tag}, H = {r['H']}")
    print(f"{'='*60}", flush=True)

    H_check, alpha = compute_ocf_dp(r['A'])
    match = "VERIFIED" if H_check == r['H'] else f"MISMATCH ({H_check} vs {r['H']})"
    print(f"\n  OCF: H = {H_check} ({match})")

    max_k = max(alpha.keys())
    for k in range(max_k + 1):
        ak = alpha.get(k, 0)
        print(f"    alpha_{k} = {ak:>10}  (2^{k} * alpha_{k} = {ak * 2**k:>10})")

    r['alpha'] = alpha

# Summary table
print("\n" + "=" * 80)
print("SUMMARY TABLE")
print("=" * 80)

max_deg = max(max(r.get('alpha', {0:1}).keys()) for r in representatives)
header = f"{'S':>22} {'H':>8}"
for k in range(max_deg + 1):
    header += f" {'a'+str(k):>10}"
header += "  Paley"
print(header)
print("-" * len(header))

for r in representatives:
    a = r.get('alpha', {})
    tag = "  ***" if r['is_paley'] else ""
    line = f"{str(sorted(r['S'])):>22} {r['H']:>8}"
    for k in range(max_deg + 1):
        line += f" {a.get(k, 0):>10}"
    line += tag
    print(line)

# Contribution breakdown
print(f"\n--- Contribution breakdown ---")
for r in representatives:
    a = r.get('alpha', {})
    tag = " [PALEY]" if r['is_paley'] else ""
    total = r['H']
    parts = []
    for k in range(max_deg + 1):
        c = a.get(k, 0) * 2**k
        pct = 100.0 * c / total
        parts.append(f"2^{k}*{a.get(k,0)}={c} ({pct:.1f}%)")
    print(f"  {sorted(r['S'])}{tag}:")
    print(f"    {' + '.join(parts)} = {total}")

# Alpha_1 dominance
print(f"\n--- Alpha_1 (linear term) comparison ---")
p_a1 = paley_result.get('alpha', {}).get(1, 0)
print(f"Paley: alpha_1 = {p_a1}")
for r in representatives:
    if not r['is_paley']:
        a1 = r.get('alpha', {}).get(1, 0)
        print(f"  S={sorted(r['S'])}: alpha_1 = {a1}, diff = {p_a1 - a1:+d}")

# Check: does Paley maximize alpha_1?
all_a1 = [(r.get('alpha', {}).get(1, 0), sorted(r['S']), r['is_paley']) for r in representatives]
all_a1.sort(reverse=True)
print(f"\nAlpha_1 ranking:")
for a1, s, p in all_a1:
    tag = " ***" if p else ""
    print(f"  alpha_1 = {a1}: S = {s}{tag}")

# Does alpha_1 alone determine H ordering?
print(f"\nDoes alpha_1 ordering match H ordering?")
h_order = [sorted(r['S']) for r in representatives]  # already sorted by H desc
a1_order = [s for _, s, _ in all_a1]
print(f"  H ordering:      {h_order}")
print(f"  alpha_1 ordering: {a1_order}")
print(f"  Match: {h_order == a1_order}")

print(f"\nTotal runtime: {time.time() - t_global:.1f}s")
