# audit-5129: product law (A1) on random non-strong tournaments of orders 8..11 (ordered joins of random blocks),
# every root; component data computed on the induced sub-tournaments; U of the whole tournament by the subset
# programme (validated literally in steps 1, 1b).
import numpy as np
import acore as A
rng = np.random.default_rng(91)
nt = 0; nr = 0; bad = 0; triv = 0
for it in range(4000):
    m = int(rng.integers(2, 5))
    sizes = [int(rng.integers(1, 6)) for _ in range(m)]
    n = sum(sizes)
    if n < 8 or n > 11:
        continue
    adj = np.zeros(n, np.int64)
    off = 0
    blocks = []
    for sz in sizes:
        code = int(rng.integers(0, 1 << (sz * (sz - 1) // 2))) if sz > 1 else 0
        sub = np.zeros(sz, np.int64)
        A.adj_from_code(code, sz, sub)
        for a in range(sz):
            adj[off + a] = int(sub[a]) << off
        blocks.append((off, sz))
        off += sz
    for i, (o1, s1) in enumerate(blocks):
        later = 0
        for (o2, s2) in blocks[i + 1:]:
            later |= ((1 << s2) - 1) << o2
        for a in range(s1):
            adj[o1 + a] |= later
    comp = np.zeros(n, np.int64)
    nc = A.components(adj, n, comp)
    hp = np.zeros((1 << n, n), np.int64); sp = np.zeros((1 << n, n), np.int64); Hs = np.zeros(1 << n, np.int64)
    st = np.zeros(4, np.int64)
    A.tables(adj, n, hp, sp, Hs)
    cH = []; cPi = []; cU = {}; csz = []
    for c in range(nc):
        vs = [v for v in range(n) if comp[v] == c]
        k = len(vs)
        sub = np.zeros(k, np.int64)
        for a_, v in enumerate(vs):
            for b_, w in enumerate(vs):
                if (adj[v] >> w) & 1:
                    sub[a_] |= (1 << b_)
        h2 = np.zeros((1 << k, k), np.int64); s2 = np.zeros((1 << k, k), np.int64); H2 = np.zeros(1 << k, np.int64)
        A.tables(sub, k, h2, s2, H2)
        cH.append(int(H2[(1 << k) - 1])); cPi.append(int(A.pi2_from_tables(k, H2))); csz.append(k)
        for a_, v in enumerate(vs):
            A.root_stats(k, h2, s2, H2, a_, st)
            cU[v] = int(st[2] + st[3])
    nt += 1
    D = int(Hs[(1 << n) - 1])
    if D != int(np.prod(np.array(cH, dtype=object))):
        bad += 1
    for y in range(n):
        A.root_stats(n, hp, sp, Hs, y, st)
        U = int(st[2] + st[3])
        k = int(comp[y])
        pred = cU[y]
        if csz[k] == 1:
            triv += 1
        for i in range(nc):
            if i < k:
                pred *= cH[i]
            elif i > k:
                pred *= cPi[i]
        nr += 1
        if pred != U:
            bad += 1
print("random non-strong tournaments of orders 8..11: %d, roots: %d (in a trivial component: %d), failures of (A1): %d" % (nt, nr, triv, bad))
