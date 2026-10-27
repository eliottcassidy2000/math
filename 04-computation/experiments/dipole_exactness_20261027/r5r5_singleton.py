"""HYP-9424 (1) at its first delta0 = 2 case: T = R5 => R5 (order 10, non-strong, delta0 = 2, a local minimum by
THM-4815 B4, so d(T, singleton locus) >= 3).  Find a 3-move word to a tournament with a source or a sink."""
from itertools import permutations
n = 10
R5 = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in (1, 2): R5[i][(i + d) % 5] = 1
out = [0]*n
for i in range(5):
    for j in range(5):
        if R5[i][j]: out[i] |= 1 << j; out[i+5] |= 1 << (j+5)
    for j in range(5): out[i] |= 1 << (j + 5)          # top R5 beats bottom R5
def paths(o):
    res = []
    def rec(path, used):
        if len(path) == n: res.append(path[:]); return
        nx = o[path[-1]] & ~used
        while nx:
            w = (nx & -nx).bit_length() - 1; nx &= nx - 1
            path.append(w); rec(path, used | 1 << w); path.pop()
    for s in range(n): rec([s], 1 << s)
    return res
def move(o, Q):
    full = (1 << n) - 1
    r = [(~o[u]) & full & ~(1 << u) for u in range(n)]
    for i in range(n - 1):
        u, v = Q[i], Q[i+1]; r[u] |= 1 << v; r[v] &= ~(1 << u)
    return r
def scores(o): return [bin(x).count("1") for x in o]
def d0(o): return min(min(s, n - 1 - s) for s in scores(o))
def strong(o):
    seen = 1; fr = 1
    while fr:
        nx = 0
        for v in range(n):
            if fr >> v & 1: nx |= o[v]
        nx &= ~seen; seen |= nx; fr = nx
    if seen != (1 << n) - 1: return False
    seen = 1; fr = 1
    while fr:
        nx = 0
        for v in range(n):
            if fr >> v & 1:
                for u in range(n):
                    if o[u] >> v & 1: nx |= 1 << u
        nx &= ~seen; seen |= nx; fr = nx
    return seen == (1 << n) - 1
P0 = paths(out); print("R5 => R5: delta0", d0(out), "strong", strong(out), "Hamiltonian paths", len(P0))
found = None
for Q1 in P0:
    T1 = move(out, Q1)
    if d0(T1) != 2: continue
    s1 = scores(T1)
    for Q2 in paths(T1):
        a, b = Q2[0], Q2[-1]
        if not (s1[a] == 2 or s1[b] == n - 1 - 2): continue      # delta0-descending endpoint condition (B4)
        T2 = move(T1, Q2)
        if d0(T2) == 1 and strong(T2):
            s2 = scores(T2)
            for Q3 in paths(T2):
                T3 = move(T2, Q3)
                if d0(T3) == 0: found = (Q1, Q2, Q3); break
        if found: break
    if found: break
print("3-move word to the singleton locus:", found)
