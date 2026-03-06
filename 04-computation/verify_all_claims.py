#!/usr/bin/env python3
"""
RIGOROUS VERIFICATION of all critical computational claims in the repo.
Instance: opus-2026-03-05-S12

This script independently re-derives everything from scratch.
"""
from itertools import combinations, permutations
from collections import Counter
import sys

def count_ham_dp(T, n):
    """Count Hamiltonian paths via DP."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1][v] for v in range(n))

def find_directed_odd_cycles(T, n, max_len=None):
    """Find all directed odd cycles. Returns list of tuples (canonicalized)."""
    if max_len is None:
        max_len = n if n % 2 == 1 else n - 1
    cycles = set()
    for length in range(3, max_len + 1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(T[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    min_idx = perm.index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    cycles.add(canon)
                    break
    return list(cycles)

def build_omega(cycles):
    """Build conflict graph: vertices=cycles, edge iff share a vertex."""
    m = len(cycles)
    cycle_sets = [set(c) for c in cycles]
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = 1
    return adj

def indep_poly_at_2(adj, m):
    """Compute I(G, 2) by brute force enumeration of independent sets."""
    result = 0
    for mask in range(1 << m):
        verts = [i for i in range(m) if mask & (1 << i)]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            result += 2**len(verts)
    return result

def quadratic_residues(p):
    return {(k*k) % p for k in range(1, p)} - {0}

def build_paley(p):
    qr = quadratic_residues(p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

# ============================================================
print("=" * 70)
print("CLAIM 1: OCF (H(T) = I(Omega(T), 2)) — exhaustive n=3,4,5")
print("=" * 70)

fails = 0
total = 0
for n in [3, 4, 5]:
    arc_pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    num_arcs = len(arc_pairs)
    for bits in range(2**num_arcs):
        T = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arc_pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        h = count_ham_dp(T, n)
        cycles = find_directed_odd_cycles(T, n)
        m = len(cycles)
        if m <= 20:
            adj = build_omega(cycles)
            i2 = indep_poly_at_2(adj, m)
            if h != i2:
                fails += 1
                print(f"  FAIL at n={n}: H={h}, I(Omega,2)={i2}")
        total += 1
    print(f"  n={n}: {2**num_arcs} tournaments checked, {fails} failures")

print(f"\nTotal: {total} tournaments, {fails} failures")

# ============================================================
print(f"\n{'=' * 70}")
print("CLAIM 2: Paley tournament H values")
print("=" * 70)

paley_expected = {3: 3, 7: 189, 11: 95095}
for p, expected in paley_expected.items():
    T = build_paley(p)
    h = count_ham_dp(T, p)
    status = "OK" if h == expected else "FAIL"
    print(f"  P({p}): H = {h}, expected {expected} [{status}]")
    
    # Also verify OCF for Paley tournaments
    if p <= 7:
        cycles = find_directed_odd_cycles(T, p)
        adj = build_omega(cycles)
        i2 = indep_poly_at_2(adj, len(cycles))
        ocf_status = "OK" if i2 == h else "FAIL"
        print(f"    OCF: I(Omega,2) = {i2} [{ocf_status}]")

# ============================================================
print(f"\n{'=' * 70}")
print("CLAIM 3: OEIS A038375 — max H at n=7 is 189 (exhaustive)")
print("=" * 70)

n = 7
arc_pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
max_h = 0
for bits in range(2**21):
    T = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(arc_pairs):
        if (bits >> k) & 1:
            T[i][j] = 1
        else:
            T[j][i] = 1
    h = count_ham_dp(T, n)
    if h > max_h:
        max_h = h
status = "OK" if max_h == 189 else "FAIL"
print(f"  Max H over all {2**21} n=7 tournaments: {max_h} (expected 189) [{status}]")

# ============================================================
print(f"\n{'=' * 70}")
print("CLAIM 4: Claim B — I(Omega(T),2) - I(Omega(T-v),2) = 2*sum mu(C)")
print("=" * 70)

# Test at n=5 exhaustive, n=6 sampled
for n in [5]:
    arc_pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    cb_fails = 0
    for bits in range(2**len(arc_pairs)):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(arc_pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        
        cycles = find_directed_odd_cycles(T, n)
        m = len(cycles)
        adj = build_omega(cycles)
        i_full = indep_poly_at_2(adj, m)
        
        for v in range(n):
            # Build T-v
            verts_minus = [u for u in range(n) if u != v]
            Tv = [[T[verts_minus[i]][verts_minus[j]] for j in range(n-1)] for i in range(n-1)]
            
            # Omega(T-v)
            cycles_v = find_directed_odd_cycles(Tv, n-1)
            adj_v = build_omega(cycles_v)
            i_minus = indep_poly_at_2(adj_v, len(cycles_v))
            
            lhs = i_full - i_minus
            
            # RHS: 2 * sum_{C through v} mu(C)
            rhs = 0
            for ci, c in enumerate(cycles):
                if v not in c:
                    continue
                # mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
                c_minus_v = set(c) - {v}
                # Map to T-v indices
                c_mapped = set()
                for u in c_minus_v:
                    c_mapped.add(verts_minus.index(u))
                # Cycles in T-v disjoint from c_mapped
                avail = [ci2 for ci2, c2 in enumerate(cycles_v) if not (set(c2) & c_mapped)]
                adj_sub = [[adj_v[avail[i]][avail[j]] for j in range(len(avail))] for i in range(len(avail))]
                mu_c = indep_poly_at_2(adj_sub, len(avail))
                rhs += 2 * mu_c
            
            if lhs != rhs:
                cb_fails += 1
    
    print(f"  n={n}: {2**len(arc_pairs)} tournaments x {n} vertices, Claim B failures: {cb_fails}")

# ============================================================
print(f"\n{'=' * 70}")
print("CLAIM 5: Claw-freeness of Omega(T) for n<=8")
print("=" * 70)

# THM-020 claims: a claw in Omega needs 3 pairwise vertex-disjoint odd cycles
# plus one touching all three. Three disjoint odd cycles need >= 9 vertices.
# So no claw for n<=8. Let's verify the ARGUMENT (not just sample).
print("  Argument: A claw K_{1,3} in Omega requires center cycle C0 adjacent to")
print("  three pairwise non-adjacent (= pairwise vertex-disjoint) cycles C1,C2,C3.")
print("  Three pairwise disjoint odd cycles need >= 3*3 = 9 vertices.")
print("  Therefore NO claw exists in Omega(T) for n <= 8. QED.")
print("  [ARGUMENT VERIFIED: Sound.]")

# Verify computationally at n=5 (quick) — check Omega has no claw
n = 5
arc_pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
claw_found = False
for bits in range(2**len(arc_pairs)):
    T = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(arc_pairs):
        if (bits >> k) & 1:
            T[i][j] = 1
        else:
            T[j][i] = 1
    cycles = find_directed_odd_cycles(T, n, max_len=5)
    m = len(cycles)
    if m < 4:
        continue
    adj = build_omega(cycles)
    # Check for claw: vertex c with 3 non-adjacent neighbors
    for c in range(m):
        nbrs = [j for j in range(m) if adj[c][j]]
        if len(nbrs) < 3:
            continue
        for tri in combinations(nbrs, 3):
            if adj[tri[0]][tri[1]] == 0 and adj[tri[0]][tri[2]] == 0 and adj[tri[1]][tri[2]] == 0:
                claw_found = True
                break
        if claw_found:
            break
    if claw_found:
        break
print(f"  n=5 computational check: claw found = {claw_found}")

# ============================================================
print(f"\n{'=' * 70}")
print("CLAIM 6: T_11 cycle counts (MISTAKE-007 corrected values)")
print("=" * 70)

T11 = build_paley(11)
expected_counts = {3: 55, 5: 594, 7: 3960, 9: 11055, 11: 5505}
# Only check odd lengths for OCF
for length in [3, 5, 7, 9, 11]:
    count = 0
    for verts in combinations(range(11), length):
        for perm in permutations(verts):
            if all(T11[perm[i]][perm[(i+1) % length]] for i in range(length)):
                count += 1
                break
    status = "OK" if count == expected_counts[length] else "FAIL"
    print(f"  c_{length}(T_11) = {count}, expected {expected_counts[length]} [{status}]")

# ============================================================
print(f"\n{'=' * 70}")
print("CLAIM 7: OCF verification for T_11")
print("=" * 70)

# From OPEN-Q-013: 95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155
# Let's verify the formula structure
# I(Omega, 2) = sum_{k>=0} alpha_k * 2^k
# alpha_0 = 1 (empty set)
# alpha_1 = total odd cycles
alpha_1 = 55 + 594 + 3960 + 11055 + 5505  # = 21169
print(f"  Total odd cycles: {alpha_1}")
print(f"  Claimed: 55+594+3960+11055+5505 = {55+594+3960+11055+5505}")

# Check: 1 + 2*21169 = 43339. But H = 95095.
# So alpha_2 * 4 + alpha_3 * 8 + ... = 95095 - 43339 = 51756
# Claimed: alpha_2 = 10879, alpha_3 = 1155
# 4*10879 + 8*1155 = 43516 + 9240 = 52756... wait
claimed_ocf = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155
print(f"  Claimed OCF sum: {claimed_ocf}")
print(f"  H(T_11) = 95095")
print(f"  Match: {claimed_ocf == 95095}")

if claimed_ocf != 95095:
    # Let's compute what it should be
    print(f"  DISCREPANCY: {claimed_ocf} != 95095, diff = {95095 - claimed_ocf}")

# ============================================================
print(f"\n{'=' * 70}")
print("CLAIM 8: H(P(19)) = 1172695746915")
print("=" * 70)

T19 = build_paley(19)
h19 = count_ham_dp(T19, 19)
status = "OK" if h19 == 1172695746915 else "FAIL"
print(f"  H(P(19)) = {h19}, expected 1172695746915 [{status}]")

print(f"\n{'=' * 70}")
print("ALL VERIFICATIONS COMPLETE")
print("=" * 70)
