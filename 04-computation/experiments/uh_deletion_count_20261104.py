"""uh_deletion_count_20261104.py (kind-pasteur-2026-11-04-S10) -- HYP-9464 (B).
For every class in an hcdump file (lines 'T n=.. H=.. hc=.. rows: r0 r1 ...', r_v = out-neighbour bitmask of v),
count the vertices v such that T - v has exactly one Hamiltonian cycle, and print the distribution per (n, hc).
Usage: python3 uh_deletion_count_20261104.py FILE [FILE ...]"""
import re, sys, collections


def hcycles(rows):
    n = len(rows)
    if n < 3:
        return 0
    M = 1 << n
    dp = [[0] * n for _ in range(M)]
    dp[1][0] = 1
    for m in range(1, M, 2):
        for e in range(n):
            c = dp[m][e]
            if not c:
                continue
            nx = rows[e] & ~m & (M - 1)
            while nx:
                w = (nx & -nx).bit_length() - 1
                nx &= nx - 1
                dp[m | 1 << w][w] += c
    return sum(dp[M - 1][e] for e in range(1, n) if rows[e] & 1)


def delete(rows, v):
    keep = [u for u in range(len(rows)) if u != v]
    idx = {u: i for i, u in enumerate(keep)}
    return [sum(1 << idx[w] for w in keep if rows[u] >> w & 1) for u in keep]


for path in sys.argv[1:]:
    dist = collections.defaultdict(collections.Counter)
    for line in open(path):
        if not line.startswith("T n="):
            continue
        n = int(re.search(r"n=(\d+)", line).group(1))
        hc = int(re.search(r"hc=(\d+)", line).group(1))
        rows = [int(x, 16) for x in line.split("rows:")[1].split()]
        k = sum(1 for v in range(n) if hcycles(delete(rows, v)) == 1)
        dist[(n, hc)][k] += 1
    for (n, hc), c in sorted(dist.items()):
        print(f"{path}: n={n} hc={hc} classes={sum(c.values())} #(v: T - v uniquely Hamiltonian) -> {dict(sorted(c.items()))}")
