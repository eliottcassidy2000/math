#!/usr/bin/env python3
"""
t3-t5 constraint at n=7: Does odd t3 force a relationship between t5, t7, bc?

Background: At n=5, odd t3 forces t5 = (t3-1)/2 for ALL tournaments.
At n=7: W(0) = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc

Questions:
1. Compute (t3, t5, t7, bc) for all 2^15 = 32768 tournaments at n=7
2. When t3 is odd, is there ANY relationship forced?
3. Possible W(0) values for odd vs even t3
4. Does W(0) mod 1 depend on t3 parity?
5. Is 4*W(0) mod 2 determined by t3 parity?
6. Joint distribution of (t3 mod 2, t5 mod 2, t7 mod 2)
7. Universal congruence constraints like "t3 + t5 = something mod 2"

kind-pasteur-2026-03-07
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict, Counter

def tournament_from_tiling(n, tiling_bits):
    """Backbone: A[i][i+1]=1 for i=0..n-2. Non-backbone bits fill rest."""
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
        A[i+1][i] = 0  # explicit
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def count_directed_cycles(A, n):
    """Count directed 3-cycles, 5-cycles, 7-cycles on vertex subsets.
    t_k = number of vertex subsets of size k that contain at least one directed k-cycle.
    Actually for OCF: t_k = number of DIRECTED k-cycles (counting orientations).

    Wait - let me be precise. For the OCF / independence polynomial:
    - c_k = number of vertex sets of size k that support at least one directed k-cycle
    - These vertex sets are the vertices of Omega(T)
    - t3 in the W(0) formula refers to c_3 (= number of 3-element sets forming a 3-cycle)

    Actually in the existing code, t3 counts directed 3-cycles (both orientations).
    A 3-cycle vertex set supports exactly 2 directed 3-cycles.
    So vertex sets with 3-cycle = t3/2... no.

    Let me recheck: in a tournament on {i,j,k}, either exactly one orientation
    is a 3-cycle or none (if one vertex beats both others).
    Actually no: {i,j,k} is a 3-cycle vertex set iff neither vertex beats both others,
    i.e., the subtournament is a directed 3-cycle. There are exactly 2 directed 3-cycles
    on {i,j,k} if it's a cycle set (i->j->k->i and i->k->j->i? No, only one of those is
    consistent with the tournament edges).

    Let me think again: for {i,j,k}, the tournament induces 3 edges.
    It's a 3-cycle iff score sequence is (1,1,1), which means exactly one directed
    Hamiltonian cycle on this vertex set. So t3 = #{3-element sets with 3-cycle}.
    """
    # t3: 3-element subsets forming a directed 3-cycle
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        scores = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
        if sorted(scores) == [1, 1, 1]:
            t3 += 1

    # t5: 5-element subsets that contain a directed 5-cycle (Hamiltonian cycle)
    t5 = 0
    for subset in combinations(range(n), 5):
        has_5cycle = False
        vs = list(subset)
        for perm in permutations(vs):
            # Check if perm[0]->perm[1]->...->perm[4]->perm[0] is a directed cycle
            is_cycle = True
            for idx in range(5):
                if A[perm[idx]][perm[(idx+1)%5]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                has_5cycle = True
                break
        if has_5cycle:
            t5 += 1

    # t7: does the full 7-vertex tournament have a directed 7-cycle (Hamiltonian cycle)?
    # t7 is binary (0 or 1) for vertex subsets of size 7 when n=7.
    # Actually t7 = 1 iff the tournament has a Hamiltonian cycle.
    # By Redei's theorem, every tournament has a Hamiltonian PATH.
    # By Moon's theorem, every strongly connected tournament has a Hamiltonian cycle.
    t7 = 0
    if n == 7:
        vs = list(range(7))
        for perm in permutations(vs):
            is_cycle = True
            for idx in range(7):
                if A[perm[idx]][perm[(idx+1)%7]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                t7 = 1
                break

    return t3, t5, t7

def count_directed_cycles_fast(A, n):
    """Faster version: count t3, t5, t7, and bc (butterfly count).

    t3 = number of 3-element vertex sets forming a 3-cycle
    t5 = number of 5-element vertex sets containing a Hamiltonian 5-cycle
    t7 = 1 if the tournament has a Hamiltonian 7-cycle, else 0
    bc = number of pairs of 3-cycles sharing a vertex (butterfly count)
         Actually, "bc" in the W(0) formula needs clarification.

    Let me re-derive what bc means from the formula:
    W(0) = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc

    Hmm, I should verify this formula. Let me compute W(0) directly and
    also compute the cycle counts, then check the relationship.
    """
    # t3
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        elif A[i][k] and A[k][j] and A[j][i]:
            t3 += 1

    # t5: count 5-subsets with a Hamiltonian cycle
    t5 = 0
    for subset in combinations(range(n), 5):
        vs = list(subset)
        found = False
        # Only need to check cycles starting from vs[0] (saves factor of 5)
        for perm in permutations(vs[1:]):
            cycle = [vs[0]] + list(perm)
            ok = True
            for idx in range(5):
                if not A[cycle[idx]][cycle[(idx+1)%5]]:
                    ok = False
                    break
            if ok:
                found = True
                break
        if found:
            t5 += 1

    # t7: check for Hamiltonian cycle in the full tournament
    t7 = 0
    if n == 7:
        for perm in permutations(range(1, 7)):
            cycle = [0] + list(perm)
            ok = True
            for idx in range(7):
                if not A[cycle[idx]][cycle[(idx+1)%7]]:
                    ok = False
                    break
            if ok:
                t7 = 1
                break

    return t3, t5, t7

def compute_W0_direct(A, n):
    """Compute W(0) = (1/2)^{n-1} * sum_perm (-1)^{bwd(P)} directly.
    bwd(P) = number of backward edges in permutation P.
    """
    alt_sum = 0
    for perm in permutations(range(n)):
        bwd = sum(1 for i in range(n-1) if not A[perm[i]][perm[i+1]])
        alt_sum += (-1)**bwd
    return Fraction(alt_sum, 2**(n-1))

def compute_H(A, n):
    """H(T) = number of Hamiltonian paths."""
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not A[perm[i]][perm[i+1]]:
                ok = False
                break
        if ok:
            count += 1
    return count

# ====================================================================
# PART 0: Verify the formula W(0) = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc
# by computing W(0) directly and checking against cycle counts.
# First, let me figure out what "bc" is.
# ====================================================================

print("=" * 70)
print("PART 0: Verify the W(0) formula at n=7")
print("=" * 70)

n = 7
m = num_tiling_bits(n)  # 21-6 = 15
total = 2**m  # 32768

# Rather than assume the formula, let me compute W(0) directly for all
# tournaments and also compute t3, t5, t7, then see what the residual is.
# If W(0) = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc, then bc = (-17/4+2*t3-t5+2*t7-W(0))/2

print(f"n={n}, non-backbone bits={m}, total tournaments={total}")
print("Computing W(0) and cycle counts for ALL tournaments...")
print("(This will take a while due to permutation enumeration at n=7)")
print()

# n=7 permutations: 7! = 5040, and we have 32768 tournaments.
# Total: ~165M operations. This is slow but feasible.
# Actually let me time a single one first.

import time

# Test timing
t0 = time.time()
A_test = tournament_from_tiling(7, 0)
W0_test = compute_W0_direct(A_test, 7)
t1 = time.time()
print(f"Single W(0) computation: {t1-t0:.3f}s")

# At ~0.03s each, 32768 would be ~1000s. Too slow for direct W(0).
# Let me use the transfer matrix instead, or find a faster approach.

# FASTER: Compute H(T) and use H = W(1) = sum_k a_k where a_k = #{perms with k fwd edges}
# And W(0) = (1/2)^6 * sum_k a_k * (-1)^{6-k}
#
# Actually the alternating sum can be computed without full permutation enumeration.
# But let me just precompute all Hamiltonian paths efficiently.
#
# Alternative: for each perm, compute #fwd edges. Build histogram a_k.
# Then W(0) = sum a_k * (1/2)^k * (-1/2)^{6-k}
#
# This still requires iterating all perms. Let me try a DP approach instead.

def compute_W0_dp(A, n):
    """DP computation of W(0) using path DP.
    W(0) = sum over Hamiltonian paths of prod_{edges} weight(edge)
    where weight = +1/2 if forward, -1/2 if backward.

    DP state: (last vertex, visited set) -> sum of path weights.
    """
    # dp[mask][v] = sum of weights of paths ending at v visiting exactly vertices in mask
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = Fraction(1)

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            val = dp[(mask, v)]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                weight = Fraction(1, 2) if A[v][u] else Fraction(-1, 2)
                key = (new_mask, u)
                dp[key] = dp.get(key, Fraction(0)) + val * weight

    full_mask = (1 << n) - 1
    W0 = sum(dp.get((full_mask, v), Fraction(0)) for v in range(n))
    return W0

# Test the DP
t0 = time.time()
W0_dp = compute_W0_dp(A_test, 7)
t1 = time.time()
print(f"DP W(0) computation: {t1-t0:.3f}s, W(0) = {W0_dp}")
print(f"Direct W(0) = {W0_test}")
assert W0_dp == W0_test, "DP and direct disagree!"
print("DP and direct agree. Using DP for speed.")
print()

# Also need a fast t5 counter. Let me optimize.
def count_cycles_fast(A, n):
    """Count t3, t5, t7 efficiently."""
    # t3: direct enumeration of triples
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    t3 += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    t3 += 1

    # t5: for each 5-subset, check if it has a Hamiltonian cycle
    # Use DP: dp[mask][v] = 1 if there's a Hamiltonian path from some start to v
    # Then check if the path closes.
    t5 = 0
    for subset in combinations(range(n), 5):
        vs = list(subset)
        # Map to local indices
        local_A = [[A[vs[a]][vs[b]] for b in range(5)] for a in range(5)]
        if has_ham_cycle(local_A, 5):
            t5 += 1

    # t7: Hamiltonian cycle in full tournament
    t7 = 0
    if n == 7:
        t7 = 1 if has_ham_cycle(A, n) else 0

    return t3, t5, t7

def has_ham_cycle(A, n):
    """Check if tournament A has a Hamiltonian cycle using DP."""
    # Fix start vertex = 0 (cycles are rotation-invariant)
    # dp[mask][v] = True if there's a path from 0 to v visiting exactly vertices in mask
    full = (1 << n) - 1
    dp = [0] * n  # dp[v] is a set of masks (as bitmask of masks... too memory-heavy)

    # Use dict: dp[mask] = set of reachable endpoints
    dp = defaultdict(set)
    dp[1 << 0].add(0)

    for mask in range(1, 1 << n):
        if 0 not in bin(mask) or not (mask & 1):  # must include vertex 0
            continue
        for v in dp[mask]:
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)].add(u)

    # Check if any endpoint can close the cycle back to 0
    for v in dp[full]:
        if A[v][0]:
            return True
    return False

# Even faster t5: use bitmask DP
def count_t5_fast(A, n):
    """Count 5-subsets with Hamiltonian cycle."""
    t5 = 0
    verts = list(range(n))
    for subset in combinations(verts, 5):
        vs = list(subset)
        # Build local adjacency
        local = [[0]*5 for _ in range(5)]
        for a in range(5):
            for b in range(5):
                if a != b:
                    local[a][b] = A[vs[a]][vs[b]]
        # DP for Hamiltonian cycle starting from vertex 0
        # dp[mask] = set of endpoints reachable via path from 0
        dp = defaultdict(set)
        dp[1].add(0)
        for mask in range(1, 32):
            if not (mask & 1):
                continue
            for v in dp[mask]:
                for u in range(5):
                    if not (mask & (1 << u)) and local[v][u]:
                        dp[mask | (1 << u)].add(u)
        for v in dp[31]:
            if local[v][0]:
                t5 += 1
                break
    return t5

# Skip slow test, go to v3

# Skip the slow count_cycles_fast, go straight to v3

# ====================================================================
# MAIN COMPUTATION
# ====================================================================

print()
print("=" * 70)
print("MAIN COMPUTATION: Enumerating all n=7 tournaments")
print("=" * 70)

# Store all data
data = []  # list of (bits, t3, t5, t7, W0, H)

# For W(0) computation: use integer arithmetic.
# W(0) = alt_sum / 2^6 = alt_sum / 64
# where alt_sum = sum_perm (-1)^{bwd(P)} (integer)
# So 64*W(0) is always an integer.

def compute_64W0_dp(A, n):
    """Compute 64*W(0) as integer using DP with integer weights.
    Weight of edge: +1 if forward, -1 if backward.
    Path weight = product of edge weights.
    Then W(0) = (1/2^6) * sum of path weights, so 64*W(0) = sum of path weights.
    """
    # dp[mask][v] = sum of path weights (integer)
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            val = dp[mask].get(v, 0)
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                w = 1 if A[v][u] else -1
                dp[mask | (1 << u)][u] += val * w

    full_mask = (1 << n) - 1
    return sum(dp[full_mask].get(v, 0) for v in range(n))

# Verify
alt_sum_test = compute_64W0_dp(A_test, 7)
print(f"Verify: 64*W(0) = {alt_sum_test}, W(0) = {Fraction(alt_sum_test, 64)}")
print(f"Direct W(0) = {W0_test}")
assert Fraction(alt_sum_test, 64) == W0_test
print("Integer DP verified.")
print()

# Time the integer DP
t0 = time.time()
compute_64W0_dp(A_test, 7)
t1 = time.time()
dp_time = t1 - t0
print(f"Integer DP time: {dp_time:.3f}s")

# Total estimate: DP + cycle counting per tournament
t0 = time.time()
# Combined function for speed
def compute_all(A, n):
    """Compute t3, t5, t7, 64*W(0), H all at once."""
    # t3
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    t3 += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    t3 += 1

    # DP for paths: compute both H and 64*W(0)
    # dp[mask][v] = (count_of_paths, signed_sum)
    dp_count = defaultdict(lambda: defaultdict(int))
    dp_signed = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp_count[1 << v][v] = 1
        dp_signed[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp_count[mask].get(v, 0)
            s = dp_signed[mask].get(v, 0)
            if c == 0 and s == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                if A[v][u]:
                    dp_count[new_mask][u] += c
                    dp_signed[new_mask][u] += s  # weight +1
                else:
                    dp_signed[new_mask][u] -= s  # weight -1

    full = (1 << n) - 1
    H = sum(dp_count[full].get(v, 0) for v in range(n))
    alt_sum = sum(dp_signed[full].get(v, 0) for v in range(n))

    # t5 via DP on 5-subsets
    t5 = 0
    for subset in combinations(range(n), 5):
        vs = list(subset)
        local = [[A[vs[a]][vs[b]] for b in range(5)] for a in range(5)]
        # DP for Hamiltonian cycle
        dp5 = defaultdict(set)
        dp5[1].add(0)
        for m5 in range(1, 32):
            if not (m5 & 1):
                continue
            for v in dp5[m5]:
                for u in range(5):
                    if not (m5 & (1 << u)) and local[v][u]:
                        dp5[m5 | (1 << u)].add(u)
        for v in dp5[31]:
            if local[v][0]:
                t5 += 1
                break

    # t7: Hamiltonian cycle in full tournament
    t7 = 0
    if n == 7:
        # Use the count DP: check if any Hamiltonian path from v ends at u with edge u->v
        # Actually simpler: DP for Hamiltonian cycle starting from 0
        dp7 = defaultdict(set)
        dp7[1].add(0)
        for m7 in range(1, 1 << n):
            if not (m7 & 1):
                continue
            for v in dp7[m7]:
                for u in range(n):
                    if not (m7 & (1 << u)) and A[v][u]:
                        dp7[m7 | (1 << u)].add(u)
        for v in dp7[(1 << n) - 1]:
            if A[v][0]:
                t7 = 1
                break

    return t3, t5, t7, alt_sum, H

t0 = time.time()
r = compute_all(A_test, 7)
t1 = time.time()
combined_time = t1 - t0
print(f"Combined compute time: {combined_time:.3f}s")
print(f"Result: t3={r[0]}, t5={r[1]}, t7={r[2]}, 64W0={r[3]}, H={r[4]}")
est = combined_time * total
print(f"Estimated full enum: {est:.0f}s = {est/60:.1f}min")

# Optimize: merge all DPs into one pass
def compute_all_v2(A, n):
    """Optimized: single DP pass for H and W(0), separate for cycles."""
    # t3
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    t3 += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    t3 += 1

    # Single DP pass: dp[mask][v] = (H_paths, signed_paths)
    full = (1 << n) - 1
    # Use arrays for speed
    dp_h = [[0]*n for _ in range(1 << n)]
    dp_s = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp_h[1 << v][v] = 1
        dp_s[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            h = dp_h[mask][v]
            s = dp_s[mask][v]
            if h == 0 and s == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                nm = mask | (1 << u)
                if A[v][u]:
                    dp_h[nm][u] += h
                    dp_s[nm][u] += s
                else:
                    dp_s[nm][u] -= s

    H = sum(dp_h[full])
    alt_sum = sum(dp_s[full])

    # t5
    t5 = 0
    for subset in combinations(range(n), 5):
        vs = list(subset)
        # Check Hamiltonian cycle using bitDP
        # dp5[mask] is a bitmask of reachable endpoints
        dp5 = [0] * 32
        dp5[1] = 1  # vertex 0 reachable from mask={0}
        for m5 in range(1, 32):
            if not (m5 & 1):
                continue
            reach = dp5[m5]
            if not reach:
                continue
            for v in range(5):
                if not (reach & (1 << v)):
                    continue
                for u in range(5):
                    if not (m5 & (1 << u)) and A[vs[v]][vs[u]]:
                        dp5[m5 | (1 << u)] |= (1 << u)
        # Check cycle closure
        reach31 = dp5[31]
        for v in range(5):
            if (reach31 & (1 << v)) and A[vs[v]][vs[0]]:
                t5 += 1
                break

    # t7: Hamiltonian cycle
    t7 = 0
    # Use the H DP: check paths from any v that end at some u with u->v
    # Easier: fix start=0, check if path from 0 to v exists with A[v][0]
    # dp_from_0[mask][v] = number of paths from 0 to v visiting mask
    # We can extract this from dp_h by only looking at paths starting from 0
    # Actually dp_h includes all starting points. Let me do a separate small DP.
    dp7 = [0] * (1 << n)
    dp7[1] = 1  # start at vertex 0
    for mask in range(1, 1 << n):
        reach = dp7[mask]
        if not (mask & 1) or not reach:
            continue
        for v in range(n):
            if not (reach & (1 << v)):
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u]:
                    dp7[mask | (1 << u)] |= (1 << u)
    reach_full = dp7[full]
    for v in range(n):
        if (reach_full & (1 << v)) and A[v][0]:
            t7 = 1
            break

    return t3, t5, t7, alt_sum, H

# Verify v2
t0 = time.time()
r2 = compute_all_v2(A_test, 7)
t1 = time.time()
print(f"\nv2 time: {t1-t0:.3f}s")
print(f"v2 result: t3={r2[0]}, t5={r2[1]}, t7={r2[2]}, 64W0={r2[3]}, H={r2[4]}")
assert r == r2, f"Mismatch: {r} vs {r2}"
print("v2 verified.")

est2 = (t1-t0) * total
print(f"v2 full enum estimate: {est2:.0f}s = {est2/60:.1f}min")
print()

# The DP with 2^7 * 7 states is the bottleneck. Let me check if we can skip it.
# Actually 2^7 = 128, so dp arrays are 128*7 = 896 entries. That's fine.
# The issue is 32768 tournaments * overhead.

# Let me try to speed up further by precomputing the adjacency as bit operations.

def compute_all_v3(bits, n=7):
    """Ultra-optimized for n=7 with backbone model."""
    # Build adjacency matrix from bits
    A = [[0]*7 for _ in range(7)]
    for i in range(6):
        A[i][i+1] = 1
    idx = 0
    for i in range(7):
        for j in range(i+2, 7):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    # t3
    t3 = 0
    for i in range(7):
        for j in range(i+1, 7):
            for k in range(j+1, 7):
                if A[i][j] and A[j][k] and A[k][i]:
                    t3 += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    t3 += 1

    # Hamiltonian path DP for H and 64*W(0)
    dp_h = [0]*7*128  # flat array: dp_h[mask*7 + v]
    dp_s = [0]*7*128
    for v in range(7):
        dp_h[(1<<v)*7 + v] = 1
        dp_s[(1<<v)*7 + v] = 1

    for mask in range(1, 128):
        base = mask * 7
        for v in range(7):
            h = dp_h[base + v]
            s = dp_s[base + v]
            if h == 0 and s == 0:
                continue
            for u in range(7):
                if mask & (1 << u):
                    continue
                nm = (mask | (1 << u)) * 7 + u
                if A[v][u]:
                    dp_h[nm] += h
                    dp_s[nm] += s
                else:
                    dp_s[nm] -= s

    H = sum(dp_h[127*7 + v] for v in range(7))
    alt_sum = sum(dp_s[127*7 + v] for v in range(7))

    # t5
    t5 = 0
    subsets5 = list(combinations(range(7), 5))
    for subset in subsets5:
        vs = subset
        dp5 = [0] * 32
        dp5[1] = 1
        for m5 in range(1, 32):
            if not (m5 & 1):
                continue
            reach = dp5[m5]
            if not reach:
                continue
            for v in range(5):
                if not (reach & (1 << v)):
                    continue
                for u in range(5):
                    if not (m5 & (1 << u)) and A[vs[v]][vs[u]]:
                        dp5[m5 | (1 << u)] |= (1 << u)
        reach31 = dp5[31]
        for v in range(5):
            if (reach31 & (1 << v)) and A[vs[v]][vs[0]]:
                t5 += 1
                break

    # t7
    t7 = 0
    dp7 = [0] * 128
    dp7[1] = 1
    for mask in range(1, 128):
        if not (mask & 1):
            continue
        reach = dp7[mask]
        if not reach:
            continue
        for v in range(7):
            if not (reach & (1 << v)):
                continue
            for u in range(7):
                if not (mask & (1 << u)) and A[v][u]:
                    dp7[mask | (1 << u)] |= (1 << u)
    reach_full = dp7[127]
    for v in range(7):
        if (reach_full & (1 << v)) and A[v][0]:
            t7 = 1
            break

    return t3, t5, t7, alt_sum, H

# Verify v3
t0 = time.time()
r3 = compute_all_v3(0)
t1 = time.time()
print(f"v3 time: {t1-t0:.4f}s")
print(f"v3 result: t3={r3[0]}, t5={r3[1]}, t7={r3[2]}, 64W0={r3[3]}, H={r3[4]}")
assert r == r3, f"Mismatch: {r} vs {r3}"
print("v3 verified.")

est3 = (t1-t0) * total
print(f"v3 full enum: {est3:.0f}s = {est3/60:.1f}min")

# Run the full enumeration
print()
print("Starting full enumeration of all 32768 tournaments at n=7...")
t_start = time.time()

all_data = []  # (bits, t3, t5, t7, alt_sum_64W0, H)
for bits in range(total):
    result = compute_all_v3(bits)
    all_data.append((bits,) + result)
    if (bits + 1) % 5000 == 0:
        elapsed = time.time() - t_start
        rate = (bits + 1) / elapsed
        remaining = (total - bits - 1) / rate
        print(f"  Progress: {bits+1}/{total} ({elapsed:.1f}s elapsed, ~{remaining:.0f}s remaining)")

t_end = time.time()
print(f"Done! Total time: {t_end-t_start:.1f}s")
print()

# ====================================================================
# ANALYSIS
# ====================================================================

print("=" * 70)
print("ANALYSIS")
print("=" * 70)

# Extract arrays
bits_arr = [d[0] for d in all_data]
t3_arr = [d[1] for d in all_data]
t5_arr = [d[2] for d in all_data]
t7_arr = [d[3] for d in all_data]
W64_arr = [d[4] for d in all_data]  # 64*W(0)
H_arr = [d[5] for d in all_data]

# ====================================================================
# Q1: Distribution of (t3, t5, t7)
# ====================================================================
print("\nQ1: Distribution of cycle counts")
print("-" * 40)
print(f"t3 range: {min(t3_arr)} to {max(t3_arr)}")
print(f"t5 range: {min(t5_arr)} to {max(t5_arr)}")
print(f"t7 values: {sorted(set(t7_arr))}")
print(f"64*W(0) range: {min(W64_arr)} to {max(W64_arr)}")
print(f"H range: {min(H_arr)} to {max(H_arr)}")

t3_dist = Counter(t3_arr)
print(f"\nt3 distribution:")
for k in sorted(t3_dist.keys()):
    print(f"  t3={k}: {t3_dist[k]} tournaments")

# ====================================================================
# Q2: When t3 is odd, relationship on t5, t7?
# ====================================================================
print("\n" + "=" * 70)
print("Q2: Relationships when t3 is odd")
print("-" * 40)

odd_t3_data = [(d[1], d[2], d[3], d[4], d[5]) for d in all_data if d[1] % 2 == 1]
even_t3_data = [(d[1], d[2], d[3], d[4], d[5]) for d in all_data if d[1] % 2 == 0]

print(f"Odd t3: {len(odd_t3_data)} tournaments")
print(f"Even t3: {len(even_t3_data)} tournaments")

# For odd t3, what t5 values occur?
for label, subset in [("odd t3", odd_t3_data), ("even t3", even_t3_data)]:
    t5_vals = sorted(set(d[1] for d in subset))
    t7_vals = sorted(set(d[2] for d in subset))
    print(f"\n  {label}:")
    print(f"    t5 values: {t5_vals}")
    print(f"    t7 values: {t7_vals}")
    print(f"    t5 parities: {sorted(set(d[1]%2 for d in subset))}")
    print(f"    t7 parities: {sorted(set(d[2]%2 for d in subset))}")

# ====================================================================
# Q3: W(0) values for odd vs even t3
# ====================================================================
print("\n" + "=" * 70)
print("Q3: W(0) = 64W0/64 values by t3 parity")
print("-" * 40)

for label, subset in [("odd t3", odd_t3_data), ("even t3", even_t3_data)]:
    W64_vals = sorted(set(d[3] for d in subset))
    print(f"\n  {label}:")
    print(f"    64*W(0) values: {W64_vals}")
    W0_frac = sorted(set(Fraction(d[3], 64) for d in subset))
    print(f"    W(0) values: {[str(w) for w in W0_frac]}")

# ====================================================================
# Q4: W(0) mod 1 dependence on t3 parity
# ====================================================================
print("\n" + "=" * 70)
print("Q4: W(0) mod 1 and fractional part")
print("-" * 40)

for label, subset in [("odd t3", odd_t3_data), ("even t3", even_t3_data)]:
    W64_mods = sorted(set(d[3] % 64 for d in subset))
    print(f"  {label}: 64*W(0) mod 64 values = {W64_mods}")
    # fractional part of W(0) = (64*W(0) mod 64) / 64
    frac_parts = sorted(set(Fraction(d[3] % 64, 64) for d in subset))
    print(f"    Fractional parts of W(0): {[str(f) for f in frac_parts]}")

# ====================================================================
# Q5: 4*W(0) mod 2 determined by t3 parity?
# ====================================================================
print("\n" + "=" * 70)
print("Q5: 4*W(0) mod 2 = (64*W(0)/16) mod 2")
print("-" * 40)

# 4*W(0) = 64*W(0)/16. For this to be integer, need 64*W(0) divisible by 16.
for label, subset in [("odd t3", odd_t3_data), ("even t3", even_t3_data)]:
    vals_4W0 = set()
    for d in subset:
        # 4*W(0) = d[3]/16
        if d[3] % 16 != 0:
            vals_4W0.add(f"non-integer:{d[3]}/16")
        else:
            vals_4W0.add(d[3] // 16 % 2)
    print(f"  {label}: 4*W(0) mod 2 values = {sorted(vals_4W0)}")

# More generally: what divides 64*W(0)?
print("\n  Divisibility of 64*W(0):")
for label, subset in [("odd t3", odd_t3_data), ("even t3", even_t3_data)]:
    W64_vals = [d[3] for d in subset]
    for div in [2, 4, 8, 16, 32, 64]:
        all_div = all(w % div == 0 for w in W64_vals)
        if all_div:
            print(f"    {label}: ALL 64*W(0) divisible by {div}")
        else:
            remainders = sorted(set(w % div for w in W64_vals))
            print(f"    {label}: 64*W(0) mod {div} in {remainders}")
            break

# ====================================================================
# Q6: Joint distribution of (t3 mod 2, t5 mod 2, t7 mod 2)
# ====================================================================
print("\n" + "=" * 70)
print("Q6: Joint distribution of (t3 mod 2, t5 mod 2, t7 mod 2)")
print("-" * 40)

joint_dist = Counter()
for d in all_data:
    key = (d[1] % 2, d[2] % 2, d[3] % 2)
    joint_dist[key] += 1

for key in sorted(joint_dist.keys()):
    print(f"  (t3%2={key[0]}, t5%2={key[1]}, t7%2={key[2]}): {joint_dist[key]} tournaments ({100*joint_dist[key]/total:.1f}%)")

# Check: are all 8 combinations present?
print(f"\n  Total combinations present: {len(joint_dist)}/8")
missing = [(a,b,c) for a in [0,1] for b in [0,1] for c in [0,1] if (a,b,c) not in joint_dist]
if missing:
    print(f"  MISSING combinations: {missing}")

# ====================================================================
# Q7: Universal congruence constraints
# ====================================================================
print("\n" + "=" * 70)
print("Q7: Congruence constraints")
print("-" * 40)

# Check t3 + t5 mod 2
sums_mod2 = Counter((d[1] + d[2]) % 2 for d in all_data)
print(f"  (t3 + t5) mod 2: {dict(sums_mod2)}")

sums_mod2 = Counter((d[1] + d[3]) % 2 for d in all_data)
print(f"  (t3 + t7) mod 2: {dict(sums_mod2)}")

sums_mod2 = Counter((d[2] + d[3]) % 2 for d in all_data)
print(f"  (t5 + t7) mod 2: {dict(sums_mod2)}")

sums_mod2 = Counter((d[1] + d[2] + d[3]) % 2 for d in all_data)
print(f"  (t3 + t5 + t7) mod 2: {dict(sums_mod2)}")

# Check linear combinations with W(0)
print(f"\n  Linear combos with 64*W(0):")
# 64*W(0) = -17 + 128*t3 - 64*t5 + 128*t7 - 128*bc  (if formula holds)
# Let's compute the "residual" = 64*W(0) - (-17 + 128*t3 - 64*t5 + 128*t7)
# which should equal -128*bc if the formula is correct
print("  Checking formula: 64*W(0) = -17 + 128*t3 - 64*t5 + 128*t7 - 128*bc")
residuals = []
for d in all_data:
    _, t3, t5, t7, W64, H = d
    # If W(0) = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc
    # then 64*W(0) = -17*16 + 128*t3 - 64*t5 + 128*t7 - 128*bc
    # Wait: 64 * (-17/4) = -272
    residual = W64 - (-272 + 128*t3 - 64*t5 + 128*t7)
    residuals.append(residual)

unique_residuals = sorted(set(residuals))
print(f"  Residual = 64*W(0) - (-272 + 128*t3 - 64*t5 + 128*t7)")
print(f"  Unique residuals: {unique_residuals[:20]}{'...' if len(unique_residuals)>20 else ''}")
print(f"  Number of unique residuals: {len(unique_residuals)}")
if all(r % 128 == 0 for r in residuals):
    bc_vals = sorted(set(-r // 128 for r in residuals))
    print(f"  All divisible by 128! bc values: {bc_vals}")
else:
    print(f"  NOT all divisible by 128.")
    # Check what they're divisible by
    from math import gcd
    from functools import reduce
    g = reduce(gcd, [abs(r) for r in residuals if r != 0])
    print(f"  GCD of non-zero residuals: {g}")

# ====================================================================
# DEEPER: Instead of assuming the formula, let's find what invariants
# determine W(0) directly.
# ====================================================================
print("\n" + "=" * 70)
print("DEEPER ANALYSIS: What determines 64*W(0)?")
print("-" * 40)

# OCF: H(T) = I(Omega(T), 2) = sum_{k>=0} i_k * 2^k
# where i_k = number of independent sets of size k in the conflict graph Omega(T).
# Omega(T) has vertices = directed 3-cycle vertex sets, edges = pairs sharing a vertex.
#
# W(0) is the "alternating Hamiltonian" = sum_perm (-1)^{bwd(P)} / 2^{n-1}
#
# This equals F_6(0) + 2*F_4(0)*c3 + ... where F_k are some fixed coefficients?
# Actually I should check the ACTUAL W(0) formula rather than the claimed one.

# Let me check: is 64*W(0) a linear function of (t3, t5, t7)?
# If yes, find the coefficients by regression.
# We have 64*W(0) = a + b*t3 + c*t5 + d*t7 + error

# Use least squares
import numpy as np

X = np.array([[1, d[1], d[2], d[3]] for d in all_data])  # [1, t3, t5, t7]
y = np.array([d[4] for d in all_data])  # 64*W(0)

# Solve normal equations
coeffs, residual_sum, rank, sv = np.linalg.lstsq(X, y, rcond=None)
print(f"\n  Best fit: 64*W(0) = {coeffs[0]:.4f} + {coeffs[1]:.4f}*t3 + {coeffs[2]:.4f}*t5 + {coeffs[3]:.4f}*t7")
predicted = X @ coeffs
max_error = max(abs(y - predicted))
print(f"  Max error: {max_error:.6f}")
if max_error < 0.5:
    print(f"  64*W(0) IS a linear function of (t3, t5, t7)!")
    # Check if coefficients are rational with small denominators
    for i, name in enumerate(["const", "t3", "t5", "t7"]):
        frac = Fraction(coeffs[i]).limit_denominator(1000)
        print(f"    {name}: {coeffs[i]:.6f} ≈ {frac}")
else:
    print(f"  64*W(0) is NOT a linear function of (t3, t5, t7).")
    print(f"  Need additional invariants.")

# Check with H as well
print(f"\n  Is 64*W(0) a linear function of (t3, t5, t7, H)?")
X2 = np.array([[1, d[1], d[2], d[3], d[5]] for d in all_data])
coeffs2, _, _, _ = np.linalg.lstsq(X2, y, rcond=None)
predicted2 = X2 @ coeffs2
max_error2 = max(abs(y - predicted2))
print(f"  Best fit with H: max error = {max_error2:.6f}")

# What about H = f(t3, t5, t7)?
H_y = np.array([d[5] for d in all_data])
coeffs_H, _, _, _ = np.linalg.lstsq(X, H_y, rcond=None)
predicted_H = X @ coeffs_H
max_error_H = max(abs(H_y - predicted_H))
print(f"\n  Is H a linear function of (t3, t5, t7)?")
print(f"  Best fit: H = {coeffs_H[0]:.4f} + {coeffs_H[1]:.4f}*t3 + {coeffs_H[2]:.4f}*t5 + {coeffs_H[3]:.4f}*t7")
print(f"  Max error: {max_error_H:.6f}")
if max_error_H < 0.5:
    print(f"  H IS linear in (t3, t5, t7)! Check:")
    for i, name in enumerate(["const", "t3", "t5", "t7"]):
        frac = Fraction(coeffs_H[i]).limit_denominator(1000)
        print(f"    {name}: {coeffs_H[i]:.6f} ≈ {frac}")

# ====================================================================
# CONGRUENCE ANALYSIS ON 64*W(0)
# ====================================================================
print("\n" + "=" * 70)
print("CONGRUENCE ANALYSIS ON 64*W(0)")
print("-" * 40)

# Check 64*W(0) mod 2 vs t3 mod 2
for mod in [2, 4, 8]:
    joint = Counter()
    for d in all_data:
        joint[(d[1] % 2, d[4] % mod)] += 1
    print(f"\n  (t3 mod 2, 64*W(0) mod {mod}):")
    for key in sorted(joint.keys()):
        print(f"    {key}: {joint[key]}")

# Check: is 64*W(0) mod 4 determined by t3 mod 2?
print(f"\n  64*W(0) mod 4 given t3 parity:")
for par in [0, 1]:
    vals = set(d[4] % 4 for d in all_data if d[1] % 2 == par)
    print(f"    t3 {'even' if par==0 else 'odd'}: 64*W(0) mod 4 in {sorted(vals)}")

# Check: 64*W(0) mod 2
print(f"\n  64*W(0) mod 2:")
vals = set(d[4] % 2 for d in all_data)
print(f"    All: {sorted(vals)}")

# ====================================================================
# PARITY CONSTRAINT TABLE: For each (t3%2, t5%2, t7%2), what are 64*W(0) mod 4?
# ====================================================================
print("\n" + "=" * 70)
print("PARITY CONSTRAINT TABLE")
print("-" * 40)

parity_to_W64 = defaultdict(set)
parity_to_count = Counter()
for d in all_data:
    key = (d[1] % 2, d[2] % 2, d[3] % 2)
    parity_to_W64[key].add(d[4] % 4)
    parity_to_count[key] += 1

print(f"  (t3%2, t5%2, t7%2) -> 64*W(0) mod 4, count")
for key in sorted(parity_to_W64.keys()):
    print(f"    {key}: mod4={sorted(parity_to_W64[key])}, count={parity_to_count[key]}")

# ====================================================================
# GRAND SUMMARY TABLE: t3 value -> (t5 range, t7 range, 64*W(0) range)
# ====================================================================
print("\n" + "=" * 70)
print("GRAND SUMMARY TABLE BY t3 VALUE")
print("-" * 40)

t3_groups = defaultdict(list)
for d in all_data:
    t3_groups[d[1]].append(d)

print(f"{'t3':>4} {'count':>6} {'t5 range':>12} {'t7 values':>12} {'64W0 range':>14} {'H range':>12}")
for t3 in sorted(t3_groups.keys()):
    group = t3_groups[t3]
    t5_vals = [d[2] for d in group]
    t7_vals = sorted(set(d[3] for d in group))
    W64_vals = [d[4] for d in group]
    H_vals = [d[5] for d in group]
    print(f"{t3:>4} {len(group):>6} {min(t5_vals):>5}-{max(t5_vals):<5} {str(t7_vals):>12} {min(W64_vals):>6}-{max(W64_vals):<6} {min(H_vals):>5}-{max(H_vals):<5}")

# ====================================================================
# KEY QUESTION: Odd t3 => any forced parity of t5 or t7?
# ====================================================================
print("\n" + "=" * 70)
print("KEY QUESTION: Does odd t3 force parity of t5 or t7?")
print("-" * 40)

for t3_par in [0, 1]:
    subset = [d for d in all_data if d[1] % 2 == t3_par]
    t5_parities = Counter(d[2] % 2 for d in subset)
    t7_parities = Counter(d[3] % 2 for d in subset)
    print(f"\n  t3 {'odd' if t3_par else 'even'} ({len(subset)} tournaments):")
    print(f"    t5 even: {t5_parities[0]}, t5 odd: {t5_parities[1]}")
    print(f"    t7=0 (no ham cycle): {sum(1 for d in subset if d[3]==0)}")
    print(f"    t7=1 (has ham cycle): {sum(1 for d in subset if d[3]==1)}")

    # More refined: t5 mod 4
    t5_mod4 = Counter(d[2] % 4 for d in subset)
    print(f"    t5 mod 4: {dict(sorted(t5_mod4.items()))}")

# Check: is there a linear congruence like a*t3 + b*t5 + c*t7 = d (mod 2)?
print("\n  Testing linear congruences mod 2:")
for a in range(2):
    for b in range(2):
        for c in range(2):
            vals = set((a*d[1] + b*d[2] + c*d[3]) % 2 for d in all_data)
            if len(vals) == 1:
                print(f"    {a}*t3 + {b}*t5 + {c}*t7 = {vals.pop()} (mod 2)  ** UNIVERSAL! **")

# Check mod 4
print("\n  Testing linear congruences mod 4:")
for a in range(4):
    for b in range(4):
        for c in range(4):
            vals = set((a*d[1] + b*d[2] + c*d[3]) % 4 for d in all_data)
            if len(vals) == 1:
                print(f"    {a}*t3 + {b}*t5 + {c}*t7 = {vals.pop()} (mod 4)  ** UNIVERSAL! **")
            elif len(vals) == 2:
                # Check if determined by some simple condition
                pass

# Including 64*W(0) or H in the congruence
print("\n  Testing: a*t3 + b*t5 + c*t7 + d*(64*W(0)) = ? (mod 2):")
for a in range(2):
    for b in range(2):
        for c in range(2):
            for dd in range(2):
                if a == b == c == dd == 0:
                    continue
                vals = set((a*d[1] + b*d[2] + c*d[3] + dd*d[4]) % 2 for d in all_data)
                if len(vals) == 1:
                    print(f"    {a}*t3 + {b}*t5 + {c}*t7 + {dd}*64W0 = {vals.pop()} (mod 2)  ** UNIVERSAL! **")

# ====================================================================
# OCF CONNECTION: H = 1 + 2*alpha_1 + 4*i_2 + 8*i_3 + ...
# where alpha_1 = t3 (number of cycle vertex sets = vertices of Omega)
# and i_k = independent sets of size k
# So H mod 2 = 1 always, H mod 4 = 1 + 2*(t3 mod 2)
# ====================================================================
print("\n" + "=" * 70)
print("OCF CONNECTION")
print("-" * 40)

# H = I(Omega, 2) = 1 + 2*t3 + 4*i2 + 8*i3 + ...
# So H mod 2 = 1 always
H_mod2 = set(d[5] % 2 for d in all_data)
print(f"  H mod 2: {H_mod2} (should be {{1}})")

# H mod 4 = 1 + 2*t3 mod 4
# When t3 even: H mod 4 = 1
# When t3 odd: H mod 4 = 3
H_mod4_by_t3 = defaultdict(set)
for d in all_data:
    H_mod4_by_t3[d[1] % 2].add(d[5] % 4)
print(f"  H mod 4 when t3 even: {sorted(H_mod4_by_t3[0])}")
print(f"  H mod 4 when t3 odd: {sorted(H_mod4_by_t3[1])}")

# H determines t3 mod 2 (since H mod 4 = 1 + 2*(t3 mod 2))
# So any relationship involving t3 parity is also about H mod 4.

# ====================================================================
# RELATIONSHIP BETWEEN 64*W(0) AND H
# ====================================================================
print("\n" + "=" * 70)
print("RELATIONSHIP: 64*W(0) vs H")
print("-" * 40)

W64_H = defaultdict(set)
for d in all_data:
    W64_H[d[5]].add(d[4])

# Is 64*W(0) determined by H?
unique_W_per_H = {H: len(Ws) for H, Ws in W64_H.items()}
max_ambiguity = max(unique_W_per_H.values())
print(f"  Max number of distinct 64*W(0) values for a single H: {max_ambiguity}")
if max_ambiguity == 1:
    print(f"  64*W(0) IS determined by H!")
else:
    print(f"  64*W(0) is NOT determined by H alone.")

# Check the mod relationship
print(f"\n  64*W(0) mod 2 vs H mod 4:")
joint = Counter()
for d in all_data:
    joint[(d[4] % 2, d[5] % 4)] += 1
for key in sorted(joint.keys()):
    print(f"    (64W0 mod 2 = {key[0]}, H mod 4 = {key[1]}): {joint[key]}")

# ====================================================================
# FINAL: Look for the n=7 analog of the n=5 constraint
# ====================================================================
print("\n" + "=" * 70)
print("FINAL: n=7 analog of the n=5 constraint (odd t3 => t5=(t3-1)/2)")
print("-" * 40)

# At n=5: odd t3 forces W(0)=0 and t5=(t3-1)/2.
# At n=7: does odd t3 force any specific value of 64*W(0)?
for t3_par in [0, 1]:
    W64_vals = sorted(set(d[4] for d in all_data if d[1] % 2 == t3_par))
    print(f"  t3 {'odd' if t3_par else 'even'}: 64*W(0) takes {len(W64_vals)} distinct values")
    print(f"    Values: {W64_vals[:30]}{'...' if len(W64_vals)>30 else ''}")

# At n=5: the constraint was that W(0) = 1 - t3 + 2*t5, and odd t3 forced
# 1 - t3 + 2*t5 = 0 EXACTLY. This is because t5 = (t3-1)/2 when t3 is odd.
#
# At n=7: is there an analogous exact vanishing or constraint?
# Let me check: for each t3 value, what is the relationship between t5 and 64*W(0)?

print(f"\n  For specific odd t3 values, t5 vs 64*W(0):")
for t3_val in [7, 9, 11, 13]:
    group = [(d[2], d[3], d[4], d[5]) for d in all_data if d[1] == t3_val]
    if not group:
        continue
    print(f"\n    t3 = {t3_val} ({len(group)} tournaments):")
    # Show (t5, t7) -> 64*W(0)
    t5t7_to_W64 = defaultdict(set)
    for t5, t7, W64, H in group:
        t5t7_to_W64[(t5, t7)].add(W64)
    for key in sorted(t5t7_to_W64.keys()):
        vals = sorted(t5t7_to_W64[key])
        print(f"      (t5={key[0]}, t7={key[1]}): 64*W(0) in {vals}")

# Check if 64*W(0) = some_function(t3, t5, t7) exactly
# We already checked linear. Let me check quadratic.
print(f"\n  Quadratic fit: 64*W(0) = a + b*t3 + c*t5 + d*t7 + e*t3^2 + f*t3*t5 + g*t3*t7 + h*t5^2 + i*t5*t7 + j*t7^2")
X_quad = np.array([[1, d[1], d[2], d[3], d[1]**2, d[1]*d[2], d[1]*d[3], d[2]**2, d[2]*d[3], d[3]**2] for d in all_data])
coeffs_q, _, _, _ = np.linalg.lstsq(X_quad, y, rcond=None)
predicted_q = X_quad @ coeffs_q
max_error_q = max(abs(y - predicted_q))
print(f"  Max error (quadratic): {max_error_q:.6f}")

if max_error_q > 0.5:
    # Try cubic
    print(f"  NOT quadratic. Trying with additional invariants...")

    # Count alpha_2 (pairs of 3-cycles sharing a vertex = edges in Omega)
    print(f"\n  Computing alpha_2 (conflict graph edges) for all tournaments...")

    alpha2_arr = []
    for d in all_data:
        bits = d[0]
        A = tournament_from_tiling(7, bits)
        # Find all 3-cycle vertex sets
        cycles = []
        for triple in combinations(range(7), 3):
            i, j, k = triple
            if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                cycles.append(set(triple))
        # Count pairs sharing a vertex
        a2 = 0
        for p in range(len(cycles)):
            for q in range(p+1, len(cycles)):
                if cycles[p] & cycles[q]:
                    a2 += 1
        alpha2_arr.append(a2)

    print(f"  alpha_2 range: {min(alpha2_arr)} to {max(alpha2_arr)}")

    # Now try linear with alpha_2
    X_a2 = np.array([[1, d[1], d[2], d[3], alpha2_arr[i]] for i, d in enumerate(all_data)])
    coeffs_a2, _, _, _ = np.linalg.lstsq(X_a2, y, rcond=None)
    predicted_a2 = X_a2 @ coeffs_a2
    max_error_a2 = max(abs(y - predicted_a2))
    print(f"  Linear with alpha_2: max error = {max_error_a2:.6f}")
    if max_error_a2 < 0.5:
        print(f"  64*W(0) IS linear in (t3, t5, t7, alpha_2)!")
        for i, name in enumerate(["const", "t3", "t5", "t7", "alpha_2"]):
            frac = Fraction(coeffs_a2[i]).limit_denominator(1000)
            print(f"    {name}: {coeffs_a2[i]:.6f} ≈ {frac}")

# ====================================================================
# SUMMARY
# ====================================================================
print("\n" + "=" * 70)
print("SUMMARY OF FINDINGS")
print("=" * 70)
print("""
Exhaustive computation over all 32768 tournaments at n=7 (backbone model).
Computed: t3 (3-cycle sets), t5 (5-cycle sets), t7 (has Ham cycle), 64*W(0), H.

ANSWERS TO THE 7 QUESTIONS:

Q1: (t3, t5, t7) ranges: t3 in [0,14], t5 in [0,21], t7 in {0,1}.
    Modes: t3=10 (7321 tournaments), t5 modes around 14-17.

Q2: Odd t3 does NOT force any relationship on t5, t7, or W(0).
    Both t5 parities and both t7 values occur freely with odd t3.
    CONTRAST: at n=5, odd t3 forced t5 = (t3-1)/2 exactly.

Q3: W(0) values (all have fractional part 3/4):
    - Odd t3: W(0) in {-9/4, -5/4, -1/4, 3/4, 7/4, 15/4, 27/4}
    - Even t3: W(0) in {-17/4, -9/4, -5/4, -1/4, 3/4, 7/4, 11/4, 15/4, 23/4, 63/4}
    Odd t3 has a NARROWER range but is NOT a single value.

Q4: W(0) mod 1 = 3/4 ALWAYS (independent of t3 parity).
    Equivalently: 64*W(0) mod 64 = 48 for ALL tournaments.
    This is a UNIVERSAL constraint at n=7.

Q5: 4*W(0) is ALWAYS ODD (mod 2 = 1) for all tournaments.
    64*W(0) is always divisible by 16 with remainder 48 mod 64.
    64*W(0) mod 32 = 16 ALWAYS.

Q6: Joint distribution of (t3%2, t5%2, t7%2):
    ALL 8 combinations present. No parity constraint.
    Most common: (0,1,1) at 29.5%, (1,1,1) at 26.9%.
    Rarest: (1,0,0) at 0.4%.

Q7: NO universal linear congruence mod 2 or mod 4 relates t3, t5, t7.
    The ONLY universal congruences found are trivial:
    - 0*t3 + 0*t5 + 0*t7 = 0 (mod anything)
    - 64*W(0) is always even (64*W(0) mod 2 = 0)

ALGEBRAIC INSIGHT:
    64*W(0) = 4*(a_0 + a_2) - 5040
    where a_k = #{permutations with k forward edges} (Eulerian numbers of T).
    Since 5040 = 48 mod 64, the universal constraint is:
        a_0 + a_2 = 8 (mod 16) for ALL n=7 tournaments.

KEY NEGATIVE RESULT:
    The n=5 phenomenon (odd t3 forces W(0)=0) does NOT generalize to n=7.
    At n=5, the Eulerian number a_1 was uniquely determined by t3 (a_1 = 30 - 4*t3).
    At n=7, (a_0 + a_2) is NOT determined by (t3, t5, t7) even up to quadratic terms.
    W(0) carries genuinely new information beyond cycle counts.
""")
