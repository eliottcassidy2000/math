"""Wave 61. Near-sink marker: Q = Q' + X with X beating only y.  Claim: a_1(Q,X) = H(Q) = beta + H(Q') - e_y and a_2(Q,X) = beta - g,
a_r = 0 for r >= 3, where beta = sum over paths P from y of H(Q' - V(P)), e_y = Hamiltonian paths of Q' ending at y, g = H(Q' - y)."""
import random, itertools
def Hsub(n, out, S):
    vs = [v for v in range(n) if S >> v & 1]
    if not vs: return 1
    from functools import lru_cache
    idx = {v: i for i, v in enumerate(vs)}; m = len(vs)
    dp = {}
    for i in range(m): dp[(1 << i, i)] = 1
    for mask in range(1, 1 << m):
        for i in range(m):
            c = dp.get((mask, i), 0)
            if not c: continue
            for j in range(m):
                if not mask >> j & 1 and out[vs[i]] >> vs[j] & 1:
                    dp[(mask | 1 << j, j)] = dp.get((mask | 1 << j, j), 0) + c
    return sum(dp.get(((1 << m) - 1, i), 0) for i in range(m))
def paths_from(n, out, y, S):
    """number of Hamiltonian paths of Q'[S] starting at y (y in S)"""
    vs = [v for v in range(n) if S >> v & 1]; m = len(vs); idx = {v: i for i, v in enumerate(vs)}
    dp = {(1 << idx[y], idx[y]): 1}
    for mask in range(1, 1 << m):
        for i in range(m):
            c = dp.get((mask, i), 0)
            if not c: continue
            for j in range(m):
                if not mask >> j & 1 and out[vs[i]] >> vs[j] & 1: dp[(mask | 1 << j, j)] = dp.get((mask | 1 << j, j), 0) + c
    return sum(dp.get(((1 << m) - 1, i), 0) for i in range(m))
def walks(N, out, X, rmax=8):
    """a_r: walks visiting X exactly r times, others once"""
    others = ((1 << N) - 1) & ~(1 << X); a = [0] * (rmax + 1)
    dp = {}
    for v in range(N):
        if v == X: dp[(v, 0, 1)] = 1
        else: dp[(v, 1 << v, 0)] = 1
    frontier = dict(dp)
    while frontier:
        nf = {}
        for (v, m, r), c in frontier.items():
            if m == others and r >= 1: a[r] += c
            for w in range(N):
                if w == v or not out[v] >> w & 1: continue
                if w == X:
                    if r < rmax: nf[(w, m, r + 1)] = nf.get((w, m, r + 1), 0) + c
                elif not m >> w & 1: nf[(w, m | 1 << w, r)] = nf.get((w, m | 1 << w, r), 0) + c
        frontier = nf
    return a
random.seed(3); bad = 0
for trial in range(300):
    n = random.randint(1, 6); out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < .5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    y = random.randrange(n); full = (1 << n) - 1
    H = Hsub(n, out, full); g = Hsub(n, out, full & ~(1 << y))
    beta = 0
    for S in range(1, 1 << n):
        if S >> y & 1: beta += paths_from(n, out, y, S) * Hsub(n, out, full & ~S)
    e_y = paths_from(n, [sum(1 << u for u in range(n) if out[u] >> v & 1) for v in range(n)], y, full)   # ending at y = starting at y in converse
    # build Q = Q' + X
    N = n + 1; o2 = out[:] + [1 << y]
    for v in range(n):
        if v != y: o2[v] |= 1 << n
    a = walks(N, o2, n)
    if a[1] != beta + H - e_y or a[2] != beta - g or any(a[3:]): bad += 1
print("near-sink formula checks: 300 random rooted tournaments (order <= 6), failures:", bad)
