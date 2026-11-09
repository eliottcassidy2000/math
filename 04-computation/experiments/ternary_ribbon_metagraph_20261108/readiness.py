#!/usr/bin/env python3
"""readiness.py NMAX HMAX -- near-join readiness for the ternary return forest (P9).
For every tournament T of order <= NMAX (gentourng) and every vertex t: h = H(T), R = R_t(T) = sum_{X subset V - t}
h(T[X]) s_t(T - X) (h(empty) = 1).  The near-join (two new vertices forming a 3-cycle with t, both beaten by T - t)
has count h + 2R (HYP-9453, k = 1) and is strong when T is.  Digit d in {2, 4} of the return forest 3h + d is
near-join-ready at h iff some (T, t) with H(T) = h has R = h + d/2.  Reports, for odd h <= HMAX attained at order
<= NMAX: the R-set (R <= 3h), readiness for d = 2 and d = 4, and whether a STRONG T supplies it."""
import sys, subprocess
from collections import defaultdict
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
HMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 199
def gentourng(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    res = []
    for bits in out:
        adj = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits[k] == '1': adj[i] |= 1 << j
                else: adj[j] |= 1 << i
                k += 1
        res.append(adj)
    return res
def tables(n, adj):
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, full):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            out = adj[v] & ~mask
            while out:
                u = (out & -out).bit_length() - 1; out &= out - 1
                dp[mask | (1 << u)][u] += c
    h = [sum(dp[mask]) for mask in range(full)]; h[0] = 1
    radj = [0] * n
    for p in range(n):
        o = adj[p]
        while o:
            q = (o & -o).bit_length() - 1; o &= o - 1; radj[q] |= 1 << p
    dp2 = [[0] * n for _ in range(full)]
    for v in range(n): dp2[1 << v][v] = 1
    for mask in range(1, full):
        for v in range(n):
            c = dp2[mask][v]
            if not c: continue
            out = radj[v] & ~mask
            while out:
                u = (out & -out).bit_length() - 1; out &= out - 1
                dp2[mask | (1 << u)][u] += c
    return h, dp2   # dp2[mask][v] = paths of T[mask] starting at v
def is_strong(adj, n):
    seen = 1; stack = [0]
    while stack:
        x = stack.pop(); out = adj[x] & ~seen
        while out:
            u = (out & -out).bit_length() - 1; out &= out - 1; seen |= 1 << u; stack.append(u)
    if seen != (1 << n) - 1: return False
    radj = [sum(1 << u for u in range(n) if (adj[u] >> v) & 1) for v in range(n)]
    seen = 1; stack = [0]
    while stack:
        x = stack.pop(); out = radj[x] & ~seen
        while out:
            u = (out & -out).bit_length() - 1; out &= out - 1; seen |= 1 << u; stack.append(u)
    return seen == (1 << n) - 1
Rset = defaultdict(set); Rstrong = defaultdict(set); attained = defaultdict(lambda: 99)
for n in range(1, NMAX + 1):
    for adj in gentourng(n):
        h, s = tables(n, adj)
        full = (1 << n) - 1; H = h[full]
        if H > HMAX: continue
        attained[H] = min(attained[H], n)
        st = is_strong(adj, n)
        for t in range(n):
            rest = full & ~(1 << t); R = 0; sub = rest
            while True:
                R += h[sub] * s[full & ~sub][t]
                if sub == 0: break
                sub = (sub - 1) & rest
            if R <= 3 * H:
                Rset[H].add(R)
                if st: Rstrong[H].add(R)
print(f"# near-join readiness, tournaments of order <= {NMAX}, odd h <= {HMAX}")
print("# h  born  R=h+1 (strong)  R=h+2 (strong)  |R-set <= 3h|  first R-gaps")
fail2 = []; fail4 = []
for H in sorted(Rset):
    r1 = (H + 1) in Rset[H]; r2 = (H + 2) in Rset[H]
    s1 = (H + 1) in Rstrong[H]; s2 = (H + 2) in Rstrong[H]
    gaps = [R for R in range(H, 3 * H + 1) if R not in Rset[H]][:6]
    print(f"{H:4d} {attained[H]:2d}  {int(r1)} ({int(s1)})  {int(r2)} ({int(s2)})  {len(Rset[H]):4d}  {gaps}")
    if not r1: fail2.append(H)
    if not r2: fail4.append(H)
print(f"# h with no R = h+1 (digit 2 not near-join-ready at order <= {NMAX}): {fail2}")
print(f"# h with no R = h+2 (digit 4 not near-join-ready at order <= {NMAX}): {fail4}")
