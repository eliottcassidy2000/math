#!/usr/bin/env python3
"""
Test palindromic N(a,b,j) for circulant tournaments at n=9.

Uses DP for Hamiltonian path counting (much faster than brute-force permutations).

kind-pasteur-2026-03-06-S25d
"""

from functools import lru_cache

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def compute_N_dp(T, n):
    """Compute N(a,b,j) using DP over bitmask Hamiltonian paths."""
    # dp[mask][last] = number of Ham paths visiting exactly 'mask' ending at 'last'
    # Also track consecutive pairs by position

    # N[a][b][j] = #{paths with {a,b} at positions {j,j+1}}
    N = [[[0]*(n-1) for _ in range(n)] for _ in range(n)]
    H = 0

    full = (1 << n) - 1

    # We need to enumerate paths, not just count them.
    # For N(a,b,j) we need to know the position of consecutive pairs.
    # Unfortunately DP alone doesn't easily give per-position consecutive info.
    # We need a modified DP that tracks the path length (= position of last vertex).

    # dp[mask][last] = list of counts? No - we need:
    # dp[mask][last] = number of paths of length popcount(mask)-1 visiting mask, ending at last.
    # The position of 'last' in such a path is popcount(mask)-1.

    # For consecutive pair (a,b) at position (j, j+1):
    # This means a is at position j, b is at position j+1.
    # So: there's a partial path of length j+1 visiting some mask, with second-to-last=a, last=b.
    # The mask has j+2 = popcount(mask) vertices.
    # Then this extends to a full Ham path.

    # So: N_directed(a,b,j) = sum over masks with popcount=j+2, a and b in mask, last=b, second-to-last=a:
    #   dp_prefix[mask][b] (with a as predecessor) * dp_suffix_count[full ^ mask | (1<<b)][b]
    # This is complex. Let me use a different approach.

    # dp_prefix[mask][last] = # paths visiting exactly mask, ending at last
    # dp_suffix[mask][first] = # paths visiting exactly mask, starting at first

    # Build dp_prefix
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            if (mask >> last) & 1 == 0:
                continue
            key = (mask, last)
            if key not in dp:
                continue
            cnt = dp[key]
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if T.get((last, nxt), 0) == 0:
                    continue
                nkey = (mask | (1 << nxt), nxt)
                dp[nkey] = dp.get(nkey, 0) + cnt

    # Build dp_suffix: # of paths visiting mask, starting at first
    # dp_suf[mask][first] = # Ham paths on vertex set=mask starting at first
    dp_suf = {}
    for v in range(n):
        dp_suf[(1 << v, v)] = 1

    for popcount in range(2, n + 1):
        for mask in range(1, 1 << n):
            if bin(mask).count('1') != popcount:
                continue
            for first in range(n):
                if (mask >> first) & 1 == 0:
                    continue
                total = 0
                for nxt in range(n):
                    if nxt == first:
                        continue
                    if (mask >> nxt) & 1 == 0:
                        continue
                    if T.get((first, nxt), 0) == 0:
                        continue
                    rest_mask = mask & ~(1 << first)
                    total += dp_suf.get((rest_mask, nxt), 0)
                if total > 0:
                    dp_suf[(mask, first)] = total

    # Now compute N(a,b,j):
    # A path has a at position j, b at position j+1 means:
    # prefix = path of length j+1 visiting some S, ending at a
    # edge a->b (or b->a for the other direction)
    # suffix = path of length n-j-2 visiting V\S\{b}, starting at b... wait
    # Actually: path positions 0..n-1. a at pos j, b at pos j+1.
    # prefix visits positions 0..j (j+1 vertices), ending at a.
    # suffix visits positions j+1..n-1 (n-j-1 vertices), starting at b.
    # prefix mask has j+1 = popcount vertices including a.
    # suffix mask = full ^ prefix_mask, but wait: suffix includes b.
    # prefix includes positions 0..j with vertices in prefix_mask (includes a, not b).
    # suffix includes positions j+1..n-1 with vertices starting at b.
    # suffix_mask = (full ^ prefix_mask) must include b.

    # N_directed(a,b,j) where a->b:
    # = sum over prefix_masks with popcount=j+1, a in mask, b not in mask:
    #     dp_prefix[(mask, a)] * dp_suffix[(full ^ mask, b)]
    #   provided T[(a,b)] = 1

    # N_directed(b,a,j) where b->a would need b at pos j, a at pos j+1.
    # But we want N(a,b,j) = #{paths with {a,b} at {j,j+1}} regardless of direction.
    # = N_directed(a,b,j) + N_directed(b,a,j)
    # where N_directed(x,y,j) requires T[(x,y)]=1 or edge x->y exists.

    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            if T.get((a,b), 0) == 0:
                continue
            # Count paths with a at position j, b at position j+1 (edge a->b used)
            for j in range(n-1):
                pc = j + 1  # popcount of prefix mask
                # Enumerate prefix masks with exactly pc bits, containing a, not containing b
                for prefix_mask in range(1, 1 << n):
                    if bin(prefix_mask).count('1') != pc:
                        continue
                    if (prefix_mask >> a) & 1 == 0:
                        continue
                    if (prefix_mask >> b) & 1:
                        continue
                    prefix_count = dp.get((prefix_mask, a), 0)
                    if prefix_count == 0:
                        continue
                    suffix_mask = full ^ prefix_mask
                    suffix_count = dp_suf.get((suffix_mask, b), 0)
                    if suffix_count == 0:
                        continue
                    # This contributes to N[a][b][j] AND N[b][a][j] (symmetrized)
                    N[a][b][j] += prefix_count * suffix_count
                    N[b][a][j] += prefix_count * suffix_count

    # Count H
    H = sum(dp.get((full, v), 0) for v in range(n))

    return N, H


# ============================================================
# n=9: All circulant tournaments
# ============================================================
print("=" * 70)
print("n=9: Palindromic N test for ALL circulant tournaments")
print("=" * 70)

n = 9
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

print(f"  {len(gen_sets)} circulant tournaments")

pal_count = 0
non_pal_count = 0

for gs_idx, gs in enumerate(sorted(gen_sets)):
    T = circulant_tournament(n, gs)
    N, H = compute_N_dp(T, n)

    # By circulant symmetry, check only pairs involving vertex 0
    is_pal = True
    for b in range(1, n):
        for j in range(n-1):
            if N[0][b][j] != N[0][b][n-2-j]:
                is_pal = False
                break
        if not is_pal:
            break

    if is_pal:
        pal_count += 1
        print(f"  [{gs_idx+1}/{len(gen_sets)}] gs={sorted(gs)}: H={H}, palindromic=YES")
    else:
        non_pal_count += 1
        print(f"  [{gs_idx+1}/{len(gen_sets)}] gs={sorted(gs)}: H={H}, palindromic=NO")
        for b in range(1, min(n, 3)):
            nab = [N[0][b][j] for j in range(n-1)]
            alt = sum((-1)**j * nab[j] for j in range(n-1))
            print(f"    N(0,{b}) = {nab}, alt_sum = {alt}")

print(f"\n  Palindromic: {pal_count}, Non-palindromic: {non_pal_count}")
if non_pal_count == 0:
    print("  ==> ALL n=9 circulants have palindromic N!")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
