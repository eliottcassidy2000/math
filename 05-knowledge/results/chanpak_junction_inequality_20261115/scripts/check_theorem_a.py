"""Checks of Lemma 1 and Theorem A of the junction-inequality note on random tournaments (literal counts),
and of Corollary B on all span-3 words with at most 8 vertices (needs lit_base.pkl)."""
import random, itertools, math, pickle, numpy as np
def ham(out, verts):
    verts = list(verts); n = len(verts)
    if n == 0: return 1
    idx = {v: i for i, v in enumerate(verts)}
    m = [sum(1 << idx[w] for w in verts if w != v and out[v][w]) for v in verts]
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n): dp[1 << i][i] = 1
    for S in range(1, 1 << n):
        for i in range(n):
            c = dp[S][i]
            if c:
                t = m[i] & ~S
                while t:
                    b = t & -t; dp[S | b][b.bit_length() - 1] += c; t ^= b
    return sum(dp[(1 << n) - 1])
def covers(out, verts, starts):
    """number of tuples (F_0, F_1..F_r): disjoint paths covering verts, F_0 possibly empty, F_i starting at starts[i-1]"""
    verts = list(verts); n = len(verts); idx = {v: i for i, v in enumerate(verts)}
    m = [sum(1 << idx[w] for w in verts if w != v and out[v][w]) for v in verts]
    full = (1 << n) - 1
    # P[S][i] = number of paths with vertex set S starting at i ; A[S] = number of paths (any start) on S, A[0] = 1
    st = [[0] * n for _ in range(1 << n)]
    for i in range(n): st[1 << i][i] = 1
    for S in range(1, 1 << n):
        for i in range(n):
            if not (S >> i) & 1 or S == 1 << i: continue
            R = S ^ (1 << i); t = m[i] & R; tot = 0
            while t:
                b = t & -t; tot += st[R][b.bit_length() - 1]; t ^= b
            st[S][i] = tot
    anyp = [1] + [sum(st[S]) for S in range(1, 1 << n)]
    f = {0: 1} if False else None
    cur = {S: anyp[S] for S in range(1 << n)}           # F_0 on S
    for y in starts:
        i = idx[y]; nxt = {}
        for S, c in cur.items():
            if not c or (S >> i) & 1: continue
            R = full ^ S; Tm = R
            sub = Tm
            while sub:
                if (sub >> i) & 1 and st[sub][i]:
                    nxt[S | sub] = nxt.get(S | sub, 0) + c * st[sub][i]
                sub = (sub - 1) & Tm
        cur = nxt
    return cur.get(full, 0)
random.seed(5); bad1 = badA = 0; n1 = nA = 0; worst = 0
for trial in range(300):
    n = random.randint(4, 8)
    out = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5: out[i][j] = True
            else: out[j][i] = True
    V = list(range(n)); H = ham(out, V)
    r = random.randint(1, 3); ys = random.sample(V, r)
    c = covers(out, V, ys); bound = H
    for y in ys: bound *= 2 ** sum(out[y])
    n1 += 1
    if c > bound: bad1 += 1
    k = random.randint(1, n - 1); A = set(random.sample(V, k)); Bs = set(V) - A
    A0 = {a for a in A if any(out[b][a] for b in Bs)}; B0 = {b for b in Bs if any(out[b][a] for a in A)}
    def kappa(S0, S, indeg):
        tot = 0
        for rr in range(len(S0) + 1):
            for R in itertools.combinations(S0, rr):
                pr = math.factorial(rr)
                for y in R:
                    d = sum(1 for w in S if w != y and (out[w][y] if indeg else out[y][w])); pr *= 2 ** d
                tot += pr
        return tot
    HA, HB = ham(out, A), ham(out, Bs); up = kappa(A0, A, False) * kappa(B0, Bs, True) * HA * HB
    nA += 1
    if not (HA * HB <= H <= up): badA += 1
    worst = max(worst, H / (HA * HB))
print("Lemma 1: %d random (tournament, targets): %d violations" % (n1, bad1))
print("Theorem A: %d random (tournament, partition): %d violations; largest H(T)/(H(A)H(B)) seen %.1f" % (nA, badA, worst))
B = pickle.load(open('lit_base.pkl', 'rb'))
Af = [M.astype(np.float64) for M in B['A']]; v0 = B['v0'].astype(np.float64); h = B['h'].astype(np.float64)
X = np.array([v0]); mx = 0
for m in range(2, 9):
    X = np.vstack([X @ Af[b & ((1 << min(3, m - 1)) - 1)].T for b in range(1 << min(3, m - 1))])
    Hs = X @ h; ok = Hs > 0
    mx = max(mx, float((X[ok].max(axis=1) / Hs[ok]).max()))
    zero = int((~ok).sum())
    print("vertices %d: %7d words; largest raw coordinate / H = %.1f ; words with H = 0: %d" % (m, len(X), float((X[ok].max(axis=1) / Hs[ok]).max()), zero))
print("Corollary B: largest ratio %.1f (bound K_3 = 26425)" % mx)
