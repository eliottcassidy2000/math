#!/usr/bin/env python3
"""probes.py: homogeneous readers as module insertions.  The cheapest readers are a single directed cycle through the
module: l exterior vertices a_1 -> ... -> a_l in a transitive chain with T -> a_1 and a_l -> T (all other exterior arcs
forward).  By the odd-cycle formula the count is R_l(T) = H(T) + 2 * Pi_{par}(T), where Pi_par sums H(T - P) over the
directed paths P of T whose length has the parity making |P| + l odd (the chain is transitive, so no cycle avoids T).
Check this formula against the direct count, then ask for each strong core T (orders 3-7) and l = 2..5 whether the probe
lands exactly on 3H + 2 or 3H + 4 (breaks the divisor 3 of the midpoint without leaving the window), and the overshoot
R_l - 3H.  Also the readout gcd of the probe family versus the divisors of H."""
import subprocess
from math import gcd
exec(open('nearjoin_check.py').read().split('# (1)')[0])   # helpers
def pi_parity(T):
    n = len(T); N = 1 << n; dp, H = ham_all(T)
    # number of directed paths of T covering S (any order) = H[S]; so Pi_odd = sum over S with |S| odd of H[S] * H[V-S]
    odd = sum(H[S] * H[(N-1) & ~S] for S in range(1, N) if bin(S).count('1') % 2 == 1)
    even = sum(H[S] * H[(N-1) & ~S] for S in range(1, N) if bin(S).count('1') % 2 == 0)
    return odd, even, H[N-1]
def probe(T, l):
    n = len(T); m = n + l; adj = [[0]*m for _ in range(m)]
    for i in range(n):
        for j in range(n): adj[i][j] = T[i][j]
    for i in range(n):
        for j in range(n, m): adj[i][j] = 1
    for i in range(n, m):
        for j in range(i+1, m): adj[i][j] = 1
    for i in range(n): adj[m-1][i] = 1; adj[i][m-1] = 0        # a_l -> T
    return ham_all(adj)[1][-1]
bad = 0; checks = 0; hits = {2: [], 4: []}; overs = {}; gcd_stats = []
for n in (3, 4, 5, 6, 7):
    for T in classes(n):
        if not strong(T): continue
        odd, even, H = pi_parity(T)
        g = H
        for l in (2, 3, 4, 5):
            # cycles through the module: tau -> a_j -> ... -> a_l -> tau for every j (the core beats a_1..a_{l-1}, a_l beats the core),
            # lengths 2..l+1, at most one mixed cycle (every one uses a_l): R_l = H + 2 (ceil(l/2) Pi_even + floor(l/2) Pi_odd)
            # the exterior segment of a mixed cycle is an increasing subsequence of a_1..a_l that starts at some a_j, j < l (the core
            # beats a_1..a_{l-1}) and ends at a_l (which beats the core): C(l-1, m-1) segments of size m >= 2, and the run P in the
            # core must have |P| + m odd: R_l = H + 2 sum_{m=2}^{l} C(l-1, m-1) Pi_{|P| = m+1 mod 2}
            from math import comb
            R = H + 2 * sum(comb(l-1, m-1) * (odd if (m % 2 == 0) else even) for m in range(2, l+1))
            if n + l <= 11:
                checks += 1; Rd = probe(T, l)
                if Rd != R: bad += 1
                R = Rd
            g = gcd(g, R)
            d = R - 3*H
            overs.setdefault(l, []).append(d)
            if d == 2: hits[2].append((n, H, l))
            if d == 4: hits[4].append((n, H, l))
        gcd_stats.append((H, g))
print(f"formula R_l = H + 2 sum_m C(l-1,m-1) Pi_parity(m) checked directly: {checks} cases, {bad} failures")
print("probes landing exactly on 3H + 2 (n, H, l):", hits[2][:12], "count", len(hits[2]))
print("probes landing exactly on 3H + 4:", hits[4][:12], "count", len(hits[4]))
for l in (2, 3, 4, 5):
    ds = overs[l]; print(f"  l={l}: min overshoot R - 3H = {min(ds)}, median {sorted(ds)[len(ds)//2]}, share with R < 3H + 2: {sum(1 for d in ds if d < 2)/len(ds):.3f}")
import collections
print("readout gcd of {H, R_2..R_5} equal to H (divisor never broken by these probes):", sum(1 for H, g in gcd_stats if g == H), "of", len(gcd_stats), "; gcd > 1:", sum(1 for H, g in gcd_stats if g > 1))
