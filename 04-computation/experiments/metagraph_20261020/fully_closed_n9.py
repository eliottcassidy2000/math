#!/usr/bin/env python3
"""Fully closed tournaments (every Hamiltonian path closes to a Hamiltonian cycle, i.e. H = n*HC) among all strong
tournaments of order n (gentourng -c). Known: n = 3 (C3), n = 5 (a -> C3 -> b -> a, (H,|Aut|) = (15,3)); none at 4, 6, 7, 8."""
import subprocess, sys
for n in [int(a) for a in sys.argv[1:]] or [9]:
    out = subprocess.run(['gentourng', '-c', '-q', str(n)], capture_output=True, text=True).stdout.split()
    full = (1 << n) - 1; found = []
    for line in out:
        adjm = [[0] * n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[k] == '1': adjm[i][j] = 1
                else: adjm[j][i] = 1
                k += 1
        adj = [[j for j in range(n) if adjm[i][j]] for i in range(n)]
        dp = [[0] * n for _ in range(1 << n)]
        for v in range(n): dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            row = dp[mask]
            for v in range(n):
                c = row[v]
                if not c: continue
                for w in adj[v]:
                    if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
        H = sum(dp[full])
        # Hamiltonian cycles: paths starting at 0 whose end beats 0 (count paths from 0 = dp restricted)
        dp0 = [[0] * n for _ in range(1 << n)]; dp0[1][0] = 1
        for mask in range(1, 1 << n):
            if not mask & 1: continue
            row = dp0[mask]
            for v in range(n):
                c = row[v]
                if not c: continue
                for w in adj[v]:
                    if not (mask >> w) & 1: dp0[mask | (1 << w)][w] += c
        HC = sum(dp0[full][v] for v in range(n) if adjm[v][0])
        if H == n * HC: found.append((line, H, HC))
    print(f"n = {n}: strong classes {len(out)}, fully closed (H = n*HC): {len(found)} {found[:5]}", flush=True)
