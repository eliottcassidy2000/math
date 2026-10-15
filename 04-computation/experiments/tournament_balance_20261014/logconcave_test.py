#!/usr/bin/env python3
"""Is the position distribution of a vertex, or the gap distribution pos(y) - pos(x) of a pair, log-concave over uniform
Hamiltonian paths of strong tournaments?  (Stanley's inequality gives this for linear extensions; Kahn-Saks's 3/11 bound
uses it.)  Exhaustive over strong tournaments n <= 7 (gentourng -c)."""
import subprocess, sys
from itertools import permutations

def decode(s, n):
    adj = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    return adj

def hampaths(adj, n):
    out = []
    def rec(path, used):
        if len(path) == n:
            out.append(tuple(path)); return
        last = path[-1]
        for v in range(n):
            if not used >> v & 1 and adj[last][v]:
                path.append(v); rec(path, used | 1 << v); path.pop()
    for s in range(n):
        rec([s], 1 << s)
    return out

def logconcave(a):
    nz = [i for i, x in enumerate(a) if x]
    if not nz: return True
    if any(a[i] == 0 for i in range(nz[0], nz[-1] + 1)): return False   # internal zero
    return all(a[i] * a[i] >= a[i - 1] * a[i + 1] for i in range(1, len(a) - 1))

for n in range(3, 8):
    lines = subprocess.run(['gentourng', '-q', '-c', str(n)], capture_output=True, text=True).stdout.split()
    pos_fail = gap_fail = tot = 0
    example_pos = example_gap = None
    for s in lines:
        adj = decode(s, n)
        P = hampaths(adj, n)
        tot += 1
        posd = [[0] * n for _ in range(n)]
        for p in P:
            for i, v in enumerate(p): posd[v][i] += 1
        bad = [v for v in range(n) if not logconcave(posd[v])]
        if bad:
            pos_fail += 1
            if example_pos is None: example_pos = (s, bad[0], posd[bad[0]])
        gfail = False
        for x in range(n):
            for y in range(n):
                if x == y: continue
                gd = [0] * n
                for p in P:
                    d = p.index(y) - p.index(x)
                    if d > 0: gd[d] += 1
                if not logconcave(gd[1:]):
                    gfail = True
                    if example_gap is None: example_gap = (s, x, y, gd[1:])
        gap_fail += gfail
    print(f"n={n}: strong tournaments {tot}; position distribution not log-concave in {pos_fail}; "
          f"gap distribution (pos(y)-pos(x) > 0) not log-concave in {gap_fail}")
    if example_pos: print(f"   e.g. {example_pos}")
    if example_gap: print(f"   e.g. {example_gap}")
