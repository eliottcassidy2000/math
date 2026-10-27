"""Brute-force cross-check of lemmaB_all_cycles.c at orders 7 and 9 (pure Python, permutations instead of DFS):
counts (T, v, Hamilton cycle of T - v) instances, blue rotations closed at v, and cycles with all rotations blue."""
import subprocess, sys
from itertools import permutations
for n in (7, 9):
    m = (n - 1) // 2
    lines = subprocess.run(["gentourng", "-q", f"-d{m}D{m}", str(n)], capture_output=True, text=True).stdout.split()
    inst = blue_tot = allblue = 0
    for s in lines:
        A = [[0] * n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == '1': A[i][j] = 1
                else: A[j][i] = 1
                k += 1
        for v in range(n):
            W = [u for u in range(n) if u != v]; c0 = W[0]
            for rest in permutations(W[1:]):
                C = (c0,) + rest; M = len(C)
                if not all(A[C[i]][C[(i + 1) % M]] for i in range(M)): continue
                inst += 1; rots = blue = 0
                for t in range(M):
                    if not A[C[(t - 1) % M]][v]: continue
                    rots += 1
                    P = [C[(t + i) % M] for i in range(M)] + [v]
                    sig = {P[i]: P[n - 1 - i] for i in range(n)}
                    if all(A[x][y] == A[sig[y]][sig[x]] for x in range(n) for y in range(n) if x != y): blue += 1
                blue_tot += blue; allblue += (blue == rots)
    print(f"n={n}: instances {inst}; blue rotations {blue_tot}; all-blue cycles {allblue}")
