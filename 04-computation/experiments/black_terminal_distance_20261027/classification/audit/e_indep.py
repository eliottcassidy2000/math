# Independent check of Part VII / Part VIII.4 for E_t(S,x0), |S| = 5,6,7, t = 1..3 (no shared code).
import itertools, sys
def bottom_layers(a):
    # S - x0 = R_1 => ... => R_q, blocks K1/C3, R_q = C3; x0 => R_1, R_q => x0, free on R_2..R_{q-1}
    res = []
    rest = a - 4
    def comps(r):
        if r == 0: yield []; return
        for p in (1, 3):
            if p <= r:
                for c in comps(r - p): yield [p] + c
    for c in comps(rest):
        if not c: continue
        blocks = c + [3]; q = len(blocks)
        for bits in itertools.product([0, 1], repeat=max(0, q - 2)):
            res.append((blocks, bits))
    return res
def build_S(blocks, bits):
    # vertex 0 = x0; blocks follow
    verts = [0]; B = []; nxt = 1
    for b in blocks:
        B.append(list(range(nxt, nxt + b))); nxt += b
    a = nxt; arc = [[0]*a for _ in range(a)]
    q = len(B)
    for i in range(q):
        blk = B[i]
        if len(blk) == 3:
            u, v, w = blk; arc[u][v] = arc[v][w] = arc[w][u] = 1
        for j in range(i + 1, q):
            for u in blk:
                for v in B[j]: arc[u][v] = 1
        # x0 relation
        if i == 0: xo = 1
        elif i == q - 1: xo = 0
        else: xo = bits[i - 1]
        for u in blk:
            if xo: arc[0][u] = 1
            else: arc[u][0] = 1
    return a, arc, B
def build_E(a, arcS, t):
    n = 2*a + t - 1
    layer = [0]*a + list(range(1, t)) + [t]*a
    idx_c = {i: a + i - 1 for i in range(1, t)}
    top = lambda u: a + t - 1 + u
    arc = [[0]*n for _ in range(n)]
    for u in range(a):
        for v in range(a):
            if arcS[u][v]: arc[u][v] = 1; arc[top(v)][top(u)] = 1
    jumps = [0] + [idx_c[i] for i in range(1, t)] + [top(0)]
    jumpset = set(zip(jumps, jumps[1:]))
    for u in range(n):
        for v in range(n):
            if u != v and layer[u] != layer[v]:
                if (u, v) in jumpset: arc[u][v] = 1
                elif (v, u) in jumpset: pass
                elif layer[u] > layer[v]: arc[u][v] = 1
    for u in range(n):
        for v in range(u + 1, n):
            assert arc[u][v] + arc[v][u] == 1, (u, v)
    return n, arc
def f(arc, n, X):
    return sum(1 for u in X for v in range(n) if v not in X and arc[u][v])
def all_min_cuts(arc, n):
    best = None; cuts = []
    for mask in range(1, (1 << n) - 1):
        X = {i for i in range(n) if mask >> i & 1}
        val = f(arc, n, X)
        if best is None or val < best: best = val; cuts = [X]
        elif val == best: cuts.append(X)
    return best, cuts
def is_blue(arc, n, Q):
    r = {Q[i]: Q[n - 1 - i] for i in range(n)}
    return all(arc[r[v]][r[u]] for u in range(n) for v in range(n) if u != v and arc[u][v])
def cand_paths(arc, n, X):
    # all Ham paths starting in X, ending outside X (pruned: once outside X, X must be fully visited since only lambda=1 arc leaves X)
    out = []
    def dfs(path, used):
        if len(path) == n:
            if path[-1] not in X: out.append(list(path))
            return
        u = path[-1]
        for w in range(n):
            if not (used >> w & 1) and arc[u][w]:
                if u in X and w not in X and any(not (used >> x & 1) for x in X): continue
                path.append(w); dfs(path, used | 1 << w); path.pop()
    for s in X: dfs([s], 1 << s)
    return out
def scores(arc, n): return [sum(r) for r in arc]
tot = 0
for a in (5, 6, 7):
    for blocks, bits in bottom_layers(a):
        aS, arcS, B = build_S(blocks, bits)
        for t in (1, 2, 3):
            if 2*a + t - 1 > 15: continue
            n, arc = build_E(aS, arcS, t)
            lam, cuts = all_min_cuts(arc, n)
            s = scores(arc, n); dmin = min(min(s), min(n - 1 - x for x in s))
            black = 0; ncand = 0
            for X in cuts:
                for Q in cand_paths(arc, n, X):
                    ncand += 1
                    if not is_blue(arc, n, Q): black += 1
            # Part VIII.4 certificate
            Rq = B[-1]; r, r1, r2 = Rq  # r->r1->r2->r
            mids = [aS + i - 1 for i in range(1, t)]
            topv = list(range(aS + t - 1, n))
            # Ham path of top from sigma(x0)=top(0): brute force DFS
            def hp_from(start, verts):
                vs = set(verts)
                def d(p, used):
                    if len(p) == len(vs): return list(p)
                    for w in vs:
                        if w not in used and arc[p[-1]][w]:
                            res = d(p + [w], used | {w})
                            if res: return res
                    return None
                return d([start], {start})
            Ptop = hp_from(aS + t - 1, topv)
            restv = [v for blk in B[:-1] for v in blk]
            Prest = None
            for st in restv:
                Prest = hp_from(st, restv)
                if Prest: break
            Q = [r, r1, r2, 0] + mids + Ptop + Prest
            okpath = len(set(Q)) == n and all(arc[Q[i]][Q[i + 1]] for i in range(n - 1))
            # K_Q T
            K = [[0]*n for _ in range(n)]
            pa = set(zip(Q, Q[1:]))
            for u in range(n):
                for v in range(n):
                    if u == v: continue
                    if (u, v) in pa: K[u][v] = 1
                    elif (v, u) in pa: K[u][v] = 0
                    else: K[u][v] = arc[v][u]
            lamK, _ = all_min_cuts(K, n)
            sK = scores(K, n); dK = min(min(sK), min(n - 1 - x for x in sK))
            print(f"a={a} blocks={blocks} bits={bits} t={t} n={n}: lambda={lam} delta={dmin} #mincuts={len(cuts)} candidates={ncand} black={black} | VIII.4 path ok={okpath} black={not is_blue(arc,n,Q)} lambda(K)={lamK} delta(K)={dK}")
            tot += 1
print("members checked:", tot)
