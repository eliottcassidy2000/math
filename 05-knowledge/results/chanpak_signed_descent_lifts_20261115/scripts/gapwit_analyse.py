import sys, itertools
from functools import lru_cache
def H(adj, vs):
    vs = list(vs); n = len(vs)
    if n == 0: return 1
    idx = {v: i for i, v in enumerate(vs)}
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n): dp[1 << i][i] = 1
    for S in range(1, 1 << n):
        for i in range(n):
            c = dp[S][i]
            if c:
                for j in range(n):
                    if not (S >> j) & 1 and adj[vs[i]][vs[j]]: dp[S | (1 << j)][j] += c
    return sum(dp[(1 << n) - 1])
def sccs(adj, vs):
    vs = list(vs); comp = []; rem = set(vs)
    def reach(s, fw):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for w in rem:
                if w not in seen and (adj[u][w] if fw else adj[w][u]): seen.add(w); st.append(w)
        return seen
    while rem:
        s = next(iter(rem)); c = reach(s, True) & reach(s, False); comp.append(sorted(c)); rem -= c
    # order by condensation (source first): number of vertices beaten outside
    comp.sort(key=lambda c: -sum(1 for w in vs if w not in c and adj[c[0]][w]))
    return comp
seen = {}
for line in open('gap3_wit_order8.txt'):
    t = line.split(); d, h = int(t[1]), int(t[2]); masks = [int(z) for z in t[4:12]]; X = int(t[12][2:])
    n = 9; adj = [[0] * n for _ in range(n)]
    for i in range(8):
        for j in range(8):
            if (masks[i] >> j) & 1: adj[i][j] = 1
    for j in range(8):
        if (X >> j) & 1: adj[8][j] = 1
        else: adj[j][8] = 1
    comps = sccs(adj, range(n)); big = [c for c in comps if len(c) > 1]
    xc = [c for c in comps if 8 in c][0]
    key = (d, h, tuple(len(c) for c in big), len(xc))
    if key in seen: continue
    seen[key] = 1
    desc = []
    for c in big:
        hc = H(adj, c); par = sorted(set(H(adj, [v for v in c if v != u]) for u in c))
        fb1 = any(H(adj, [v for v in c if v != u]) == 1 for u in c)
        desc.append("order %d H=%d%s parents %s%s" % (len(c), hc, " (contains x)" if 8 in c else "", par, " one-apex" if fb1 else ""))
    if len(xc) > 1:
        c = xc; sub = [v for v in c if v != 8]; cs = sccs(adj, sub)
        desc.append("core minus x: components " + ", ".join("%d:H=%d" % (len(q), H(adj, q)) for q in cs))
        desc.append("x out-degree in core %d of %d" % (sum(adj[8][v] for v in sub), len(sub)))
    print(d, h, "|", " ; ".join(desc))
