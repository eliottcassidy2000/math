#!/usr/bin/env python3
"""verify_md_py.py - independent pure-Python module defect, lambda and Hamiltonian-path count for one tournament
(paleyeven lane; used to re-check D8_144).  usage: python3 verify_md_py.py "n code" """
import sys
n, code = sys.argv[1].split()[:2]; n = int(n)
out = [0]*n; k = 0
for i in range(n):
    for j in range(i+1, n):
        if code[k] == '1': out[i] |= 1 << j
        else: out[j] |= 1 << i
        k += 1
inn = [0]*n
for a in range(n):
    for b in range(n):
        if out[a] >> b & 1: inn[b] |= 1 << a
full = (1 << n) - 1
md = lam = 10**9
pc = int.bit_count
for M in range(1, full):
    cut = sum(pc(out[z] & ~M & full) for z in range(n) if M >> z & 1)
    lam = min(lam, cut)
    if pc(M) < 2: continue
    d = 0
    for z in range(n):
        if not (M >> z & 1):
            a = pc(out[z] & M); b = pc(inn[z] & M); d += a if a < b else b
    md = min(md, d)
# Hamiltonian paths by DP over subsets
dp = [[0]*n for _ in range(1 << n)]
for a in range(n): dp[1 << a][a] = 1
for S in range(1, 1 << n):
    row = dp[S]
    for a in range(n):
        c = row[a]
        if c:
            nx = out[a] & ~S
            while nx:
                b = (nx & -nx).bit_length() - 1; nx &= nx - 1
                dp[S | (1 << b)][b] += c
print("n", n, "md", md, "lambda", lam, "H", sum(dp[full]))
