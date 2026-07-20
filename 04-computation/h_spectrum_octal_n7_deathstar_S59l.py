#!/usr/bin/env python3
"""
death-star-2026-07-19-S59l -- HYP-8045(iii): the octal H-spectrum sieve at n=7
(dormant testable from T832/mac-mini-S1, surfaced by the concept sweep).
Gaps in odd[1,45] at n<=6 are {7,21,35,39}; 7 and 21 are proved/strongly-conj
PERMANENT (THM-343, THM-115); 35 and 39 are conjectured to first appear at n=7
but never checked. Test: large-sample n=7 tournaments, exact H by Held-Karp DP,
report which odd values in [1,60] are realized.
"""
import random, time
random.seed(59)
N = 7
def H_count(adj):
    # Held-Karp: paths[(mask, last)] = number of ham paths on mask ending at last
    dp = [[0]*N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1 << N):
        for last in range(N):
            c = dp[mask][last]
            if not c: continue
            for nxt in range(N):
                if mask & (1 << nxt): continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    full = (1 << N) - 1
    return sum(dp[full][v] for v in range(N))

seen = {}
t0 = time.time()
trials = 250000
for it in range(trials):
    adj = [[0]*N for _ in range(N)]
    for i in range(N):
        for j in range(i+1, N):
            if random.random() < 0.5: adj[i][j] = 1
            else: adj[j][i] = 1
    h = H_count(adj)
    seen[h] = seen.get(h, 0) + 1
odd_seen = sorted(v for v in seen if v % 2 == 1 and v <= 60)
missing = [v for v in range(1, 61, 2) if v not in seen]
print(f"n=7 random sample {trials} tournaments ({time.time()-t0:.0f}s)")
print(f"odd H values realized (<=60): {odd_seen}")
print(f"odd H values MISSING (<=60): {missing}")
print(f"35 realized: {35 in seen} ({seen.get(35,0)}x); 39 realized: {39 in seen} ({seen.get(39,0)}x)")
print(f"7 absent: {7 not in seen}; 21 absent: {21 not in seen}")
print(f"mod-8 residues of missing: {sorted(set(v % 8 for v in missing))}")
