import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
dc_ocf_fast.py
kind-pasteur-2026-03-07-S39b

Fast version of DC-OCF tracking. Uses bitmask DP for cycle finding
instead of permutation enumeration.

Key question: How does I(Omega(T), 2) relate to I(Omega(T\e), 2) + I(Omega(T/e), 2)?

Since H(T) = H(T\e) + H(T/e) and H(T) = I(Omega(T), 2) for tournaments,
but OCF may NOT hold for the non-tournament T\e, we track the discrepancy.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from collections import defaultdict
from itertools import combinations


def ham_paths_dp(adj, n):
    if n <= 1:
        return 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def find_odd_cycles_dp(adj, n):
    """Find all directed odd cycles using bitmask DP. Returns list of frozensets."""
    results = []
    for start in range(n):
        dp = [[0] * n for _ in range(1 << n)]
        dp[1 << start][start] = 1
        for mask in range(1 << start, 1 << n):
            if not (mask & (1 << start)):
                continue
            if mask & ((1 << start) - 1):
                continue
            popcount = bin(mask).count('1')
            for last in range(n):
                if not (mask & (1 << last)):
                    continue
                cnt = dp[mask][last]
                if cnt == 0:
                    continue
                if popcount >= 3 and popcount % 2 == 1 and last != start:
                    if adj[last][start]:
                        verts = frozenset(i for i in range(n) if mask & (1 << i))
                        for _ in range(cnt):
                            results.append(verts)
                if popcount < n:
                    for nxt in range(start + 1, n):
                        if mask & (1 << nxt):
                            continue
                        if adj[last][nxt]:
                            dp[mask | (1 << nxt)][nxt] += cnt
    return results


def independence_poly_at_2(cycle_vsets):
    """Compute I(Omega, 2) from list of cycle vertex sets."""
    if not cycle_vsets:
        return 1
    m = len(cycle_vsets)
    if m > 20:
        # Too many cycles for brute force
        return _ip_at_2_bounded(cycle_vsets, m)
    adj = [0] * m
    for a in range(m):
        for b in range(a+1, m):
            if cycle_vsets[a] & cycle_vsets[b]:
                adj[a] |= 1 << b
                adj[b] |= 1 << a
    total = 1
    for mask in range(1, 1 << m):
        bits = []
        temp = mask
        while temp:
            b = temp & (-temp)
            bits.append(b.bit_length() - 1)
            temp ^= b
        is_indep = True
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if adj[bits[i]] & (1 << bits[j]):
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            total += 2 ** len(bits)
    return total


def _ip_at_2_bounded(cycle_vsets, m):
    """Independence polynomial at 2 for large cycle sets (bounded enumeration)."""
    max_indep = max(len(set(range(100))) for _ in [0])  # dummy
    # Use greedy bounding
    adj = {}
    for a in range(m):
        adj[a] = set()
        for b in range(m):
            if a != b and cycle_vsets[a] & cycle_vsets[b]:
                adj[a].add(b)

    # Count up to size 4
    total = 1 + 2 * m
    # Size 2
    pairs = 0
    for a in range(m):
        for b in range(a+1, m):
            if b not in adj[a]:
                pairs += 1
    total += 4 * pairs
    # Size 3 (approximate for large m)
    # Skip for now
    return total


def contract_edge(T, u, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    nn = len(verts)
    D = [[0]*nn for _ in range(nn)]
    for i_idx, i in enumerate(verts):
        for j_idx, j in enumerate(verts):
            if i_idx == j_idx:
                continue
            if i == u:
                D[i_idx][j_idx] = T[v][j]
            elif j == u:
                D[i_idx][j_idx] = T[i][u]
            else:
                D[i_idx][j_idx] = T[i][j]
    return D


def delete_edge(T, u, v):
    n = len(T)
    D = [row[:] for row in T]
    D[u][v] = 0
    return D


# ============================================================
# Main
# ============================================================
print("=" * 70)
print("DC-OCF FAST TRACKING")
print("=" * 70)

for n in range(3, 6):  # n=3,4,5 only (n=6 too slow in Python)
    m = n * (n - 1) // 2
    print(f"\nn={n} (exhaustive, {1 << m} tournaments)")

    relationships = defaultdict(int)
    ocf_holds_del = 0
    ocf_fails_del = 0
    total = 0

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        H_T = hamiltonian_path_count(T)

        # Pick first directed edge
        eu, ev = None, None
        for u in range(n):
            for v in range(n):
                if u != v and T[u][v]:
                    eu, ev = u, v
                    break
            if eu is not None:
                break

        # T: OCF holds (tournament)
        cycles_T = find_odd_cycles_dp(T, n)
        I_T = independence_poly_at_2(cycles_T)
        assert I_T == H_T, f"OCF fails: bits={bits}, H={H_T}, I={I_T}"

        # T\e: deletion (NOT a tournament)
        T_del = delete_edge(T, eu, ev)
        H_del = ham_paths_dp(T_del, n)
        cycles_del = find_odd_cycles_dp(T_del, n)
        I_del = independence_poly_at_2(cycles_del)
        if H_del == I_del:
            ocf_holds_del += 1
        else:
            ocf_fails_del += 1

        # T/e: contraction
        T_con = contract_edge(T, eu, ev)
        H_con = ham_paths_dp(T_con, n - 1)
        cycles_con = find_odd_cycles_dp(T_con, n - 1)
        I_con = independence_poly_at_2(cycles_con)

        # Key identity: I(T) = I(T\e) + I(T/e)?
        diff = I_T - (I_del + I_con)
        relationships[diff] += 1

        # Also: H(T\e) = I(T) - I(T/e)?
        # Since H(T) = H(T\e) + H(T/e) and H(T) = I(T):
        # H(T\e) = I(T) - H(T/e)
        # If OCF holds for T/e: H(T/e) = I(T/e), so H(T\e) = I(T) - I(T/e)
        # But OCF may not hold for T/e when it's not a tournament.

        total += 1

    print(f"  OCF for T\\e: {ocf_holds_del}/{total} ({100*ocf_holds_del/total:.1f}%)")
    print(f"  I(T) = I(T\\e) + I(T/e): {relationships.get(0,0)}/{total} ({100*relationships.get(0,0)/total:.1f}%)")
    print(f"  Diff distribution: {dict(sorted(relationships.items()))}")

    # Now check: what is H(T\e) - I(T\e) when OCF fails for T\e?
    print(f"\n  OCF error for T\\e analysis:")
    err_dist = defaultdict(int)
    for bits in range(min(1 << m, 1024)):  # limit for speed
        T = tournament_from_bits(n, bits)
        eu, ev = None, None
        for u in range(n):
            for v in range(n):
                if u != v and T[u][v]:
                    eu, ev = u, v
                    break
            if eu is not None:
                break
        T_del = delete_edge(T, eu, ev)
        H_del = ham_paths_dp(T_del, n)
        cycles_del = find_odd_cycles_dp(T_del, n)
        I_del = independence_poly_at_2(cycles_del)
        err = H_del - I_del
        err_dist[err] += 1

    print(f"    H(T\\e) - I(Omega(T\\e), 2) distribution: {dict(sorted(err_dist.items()))}")

# ============================================================
# CRITICAL: Does deleting a tournament edge change Omega?
# ============================================================
print("\n" + "=" * 70)
print("HOW DOES DELETING EDGE e CHANGE Omega?")
print("=" * 70)

n = 5
m = n * (n - 1) // 2

# For each tournament and its first edge, compare cycles of T vs T\e
cycle_diff_stats = defaultdict(int)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    eu, ev = None, None
    for u in range(n):
        for v in range(n):
            if u != v and T[u][v]:
                eu, ev = u, v
                break
        if eu is not None:
            break

    cycles_T = find_odd_cycles_dp(T, n)
    T_del = delete_edge(T, eu, ev)
    cycles_del = find_odd_cycles_dp(T_del, n)

    # How many cycles are lost?
    # Cycles using edge eu->ev are lost; all others remain
    lost = len(cycles_T) - len(cycles_del)
    cycle_diff_stats[lost] += 1

print(f"\nn={n}: Cycles lost when deleting first edge:")
for lost, count in sorted(cycle_diff_stats.items()):
    print(f"  Lost {lost} cycles: {count} tournaments")

# Detailed: which cycles use the deleted edge?
print(f"\n  Detailed example: cyclic T_5")
T5c = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        T5c[i][(i + d) % 5] = 1

cycles_T5 = find_odd_cycles_dp(T5c, 5)
print(f"  T_5 has {len(cycles_T5)} odd cycles")
# Show vertex sets
vset_counts = defaultdict(int)
for c in cycles_T5:
    vset_counts[tuple(sorted(c))] += 1
for vs, cnt in sorted(vset_counts.items()):
    print(f"    {vs}: {cnt} directed cycle(s)")

# Delete edge 0->1
T5_del = delete_edge(T5c, 0, 1)
cycles_del = find_odd_cycles_dp(T5_del, 5)
print(f"  T_5 \\ (0->1) has {len(cycles_del)} odd cycles")
vset_counts2 = defaultdict(int)
for c in cycles_del:
    vset_counts2[tuple(sorted(c))] += 1
for vs, cnt in sorted(vset_counts2.items()):
    print(f"    {vs}: {cnt} directed cycle(s)")

# Which cycles were lost?
lost_vsets = set(vset_counts.keys()) - set(vset_counts2.keys())
print(f"  Lost vertex sets: {lost_vsets}")
for vs in lost_vsets:
    print(f"    {vs}: was {vset_counts[vs]} cycle(s)")

# Which vertex sets lost some but not all cycles?
for vs in set(vset_counts.keys()) & set(vset_counts2.keys()):
    if vset_counts[vs] != vset_counts2[vs]:
        print(f"    {vs}: {vset_counts[vs]} -> {vset_counts2[vs]} cycles")

print("\nDone.")
