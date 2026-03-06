#!/usr/bin/env python3
"""
Check: what tournament(s) achieve H=661 at n=8?
Also verify H(T_657)=657 and check if 661 appears among SC tournaments.
"""
from itertools import permutations

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

def build_sc8():
    """Build all self-converse tournaments on 8 vertices with alpha=(7,6,5,4,3,2,1,0)."""
    n = 8
    alpha = [7, 6, 5, 4, 3, 2, 1, 0]
    free_arcs = []
    determined = set()
    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in determined:
                continue
            ai, aj = alpha[i], alpha[j]
            u, v = min(aj, ai), max(aj, ai)
            if (u, v) == (i, j):
                free_arcs.append(((i, j), None))
                determined.add((i, j))
            else:
                free_arcs.append(((i, j), (aj, ai)))
                determined.add((i, j))
                determined.add((u, v))
    num_free = len(free_arcs)
    tournaments = []
    for bits in range(1 << num_free):
        T = [[0]*n for _ in range(n)]
        for k, (primary, linked) in enumerate(free_arcs):
            val = (bits >> k) & 1
            i, j = primary
            T[i][j] = val; T[j][i] = 1 - val
            if linked:
                a, b = linked
                T[a][b] = val; T[b][a] = 1 - val
        tournaments.append(T)
    return tournaments

# Build all SC tournaments with this specific anti-automorphism
print("Building SC tournaments with alpha=(7,6,5,4,3,2,1,0)...")
sc_tours = build_sc8()
print(f"Total: {len(sc_tours)}")

# Find max H among SC tournaments
max_h = 0
max_tours = []
h_dist = {}
for T in sc_tours:
    h = count_ham_dp(T, 8)
    h_dist[h] = h_dist.get(h, 0) + 1
    if h > max_h:
        max_h = h
        max_tours = [T]
    elif h == max_h:
        max_tours.append(T)

print(f"\nMax H among SC(alpha=reverse) tournaments: {max_h}")
print(f"Number achieving max: {len(max_tours)}")

# Check if 661 appears
if 661 in h_dist:
    print(f"\nH=661 count: {h_dist[661]}")
else:
    print(f"\nH=661 NOT found among these SC tournaments")

# Show top H values
top = sorted(h_dist.keys(), reverse=True)[:10]
print(f"\nTop 10 H values: {top}")
for h in top[:5]:
    print(f"  H={h}: {h_dist[h]} tournaments")

# Check the max tournament's properties
if max_tours:
    T = max_tours[0]
    scores = sorted([sum(T[i][j] for j in range(8) if j!=i) for i in range(8)])
    print(f"\nMax tournament score sequence: {scores}")
    print("Adjacency matrix:")
    for row in T:
        print(f"  {''.join(str(x) for x in row)}")
