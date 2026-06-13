#!/usr/bin/env python3
"""
Deep investigation of H=661 tournament at n=8.
- Is it the global maximum (not just among SC)?
- What are its properties (Aut, contains P(7)?, D_v uniformity)?
- How does it compare to T_657?
"""
from itertools import combinations, permutations
from collections import Counter

def count_ham_dp(T, n):
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

def find_odd_cycles_3(T, n):
    """Find all directed 3-cycles."""
    cycles = []
    for a, b, c in combinations(range(n), 3):
        for perm in permutations([a, b, c]):
            if T[perm[0]][perm[1]] and T[perm[1]][perm[2]] and T[perm[2]][perm[0]]:
                canon = min((perm[i:] + perm[:i]) for i in range(3))
                cycles.append(tuple(canon))
                break
    return list(set(cycles))

def is_isomorphic(T1, T2, n):
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if T1[i][j] != T2[perm[i]][perm[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return perm
    return None

# The H=661 tournament from the search
T661 = [
    [0,0,0,0,1,1,1,0],
    [1,0,0,0,0,1,0,1],
    [1,1,0,0,0,0,1,1],
    [1,1,1,0,0,0,0,1],
    [0,1,1,1,0,0,0,0],
    [0,0,1,1,1,0,0,0],
    [0,1,0,1,1,1,0,0],
    [1,0,0,0,1,1,1,0],
]

# T_657 from previous work
T657 = [
    [0,0,1,0,0,1,1,0],
    [1,0,0,0,0,1,0,1],
    [0,1,0,0,0,0,1,1],
    [1,1,1,0,0,0,0,0],
    [1,1,1,1,0,0,0,0],
    [0,0,1,1,1,0,0,1],
    [0,1,0,1,1,1,0,0],
    [1,0,0,1,1,0,1,0],
]

n = 8

# Verify H values
h661 = count_ham_dp(T661, n)
h657 = count_ham_dp(T657, n)
print(f"H(T661) = {h661}")
print(f"H(T657) = {h657}")

# Score sequences
scores661 = [sum(T661[i][j] for j in range(n) if j!=i) for i in range(n)]
scores657 = [sum(T657[i][j] for j in range(n) if j!=i) for i in range(n)]
print(f"\nT661 scores: {scores661} (sorted: {sorted(scores661)})")
print(f"T657 scores: {scores657} (sorted: {sorted(scores657)})")

# Build P(7)
P7 = [[0]*7 for _ in range(7)]
qr = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in qr:
            P7[i][j] = 1

# Check T661 vertex deletions
print(f"\nT661 vertex deletions:")
for v in range(n):
    verts = [u for u in range(n) if u != v]
    subT = [[T661[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
    sub_scores = [sum(subT[i][j] for j in range(7) if j!=i) for i in range(7)]
    sub_h = count_ham_dp(subT, 7)
    is_regular = all(d == 3 for d in sub_scores)
    iso_p7 = is_isomorphic(subT, P7, 7) if is_regular else None
    marker = " <-- P(7)!" if iso_p7 else ""
    print(f"  T661-{v}: scores={sorted(sub_scores)}, H={sub_h}, regular={is_regular}{marker}")

# Check isomorphism between T661 and T657
print(f"\nT661 isomorphic to T657? {is_isomorphic(T661, T657, n) is not None}")

# Automorphism group
print(f"\nT661 automorphisms:")
auts = []
for perm in permutations(range(n)):
    ok = True
    for i in range(n):
        for j in range(i+1, n):
            if T661[i][j] != T661[perm[i]][perm[j]]:
                ok = False
                break
        if not ok:
            break
    if ok:
        auts.append(perm)
print(f"  |Aut(T661)| = {len(auts)}")
for a in auts:
    if a != tuple(range(n)):
        print(f"  {a}")

# 3-cycles and mu analysis
cycles661 = find_odd_cycles_3(T661, n)
cycles657 = find_odd_cycles_3(T657, n)
print(f"\nT661: {len(cycles661)} 3-cycles")
print(f"T657: {len(cycles657)} 3-cycles")

# D_v for T661
print(f"\nD_v (sum of mu over 3-cycles through v) for T661:")
for v in range(n):
    # For each 3-cycle through v, compute mu = I(Omega_complement, 2)
    Dv = 0
    for cyc in cycles661:
        if v not in cyc:
            continue
        complement_verts = [u for u in range(n) if u not in cyc]
        # Find 3-cycles in complement
        comp_cycles = []
        for a, b, c in combinations(complement_verts, 3):
            for perm in permutations([a, b, c]):
                if T661[perm[0]][perm[1]] and T661[perm[1]][perm[2]] and T661[perm[2]][perm[0]]:
                    canon = min((perm[i:] + perm[:i]) for i in range(3))
                    comp_cycles.append(tuple(canon))
                    break
        comp_cycles = list(set(comp_cycles))
        # Build Omega on complement cycles
        m = len(comp_cycles)
        cycle_sets = [set(c) for c in comp_cycles]
        adj = [[0]*m for _ in range(m)]
        for i in range(m):
            for j in range(i+1, m):
                if cycle_sets[i] & cycle_sets[j]:
                    adj[i][j] = adj[j][i] = 1
        # Independence polynomial at x=2 (restricted to 3-cycles)
        mu = 0
        for mask in range(1 << m):
            vv = [i for i in range(m) if mask & (1 << i)]
            ok = True
            for i in range(len(vv)):
                for j in range(i+1, len(vv)):
                    if adj[vv[i]][vv[j]]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                mu += 2**len(vv)
        Dv += mu
    print(f"  D_{v} = {Dv}")

# Check self-converse with alpha = reverse
alpha = [7,6,5,4,3,2,1,0]
T661_rev = [[T661[alpha[j]][alpha[i]] for j in range(n)] for i in range(n)]
print(f"\nT661 self-converse (alpha=reverse)? {T661 == T661_rev}")
