"""Probe: sole-exit-free tournaments with exactly 3 Hamiltonian cycles (order n) vs uniquely Hamiltonian classes of
order n-1 with two sole exits / two sole entries (both counted by F(2n-9))."""
import re, sys, collections, itertools

def parse(path):
    out = []
    for l in open(path):
        if l.startswith("T n="):
            n = int(re.search(r"n=(\d+)", l).group(1))
            rows = [int(x, 16) for x in l.split("rows:")[1].split()]
            out.append((n, rows))
    return out

def sub(rows, keep):
    idx = {v: i for i, v in enumerate(keep)}
    return [sum(1 << idx[w] for w in keep if rows[v] >> w & 1) for v in keep]

def hcycles(rows, cap=10**9):
    n = len(rows)
    if n < 3: return 0
    M = 1 << n
    dp = [[0] * n for _ in range(M)]
    dp[1][0] = 1
    for m in range(1, M, 2):
        for e in range(n):
            c = dp[m][e]
            if not c: continue
            nx = rows[e] & ~m & (M - 1)
            while nx:
                w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                dp[m | 1 << w][w] += c
    return sum(dp[M - 1][e] for e in range(1, n) if rows[e] & 1)

def unique_cycle(rows):
    """return the unique Hamiltonian cycle (as vertex list starting at 0) of a UH tournament"""
    n = len(rows)
    path = [0]; seen = {0}
    res = []
    def dfs(v):
        if len(path) == n:
            if rows[v] & 1: res.append(list(path))
            return
        nx = rows[v]
        for w in range(n):
            if nx >> w & 1 and w not in seen:
                seen.add(w); path.append(w); dfs(w); path.pop(); seen.discard(w)
                if len(res) > 1: return
    dfs(0)
    assert len(res) == 1
    return res[0]

def canon_uh(rows):
    """canonical form of a UH tournament: min over rotations of the adjacency in cycle order"""
    n = len(rows); c = unique_cycle(rows); best = None
    for r in range(n):
        order = c[r:] + c[:r]
        pos = {v: i for i, v in enumerate(order)}
        code = tuple(sum(1 << pos[w] for w in range(n) if rows[v] >> w & 1) for v in order)
        if best is None or code < best: best = code
    return best

def degs(rows):
    n = len(rows); od = [bin(r).count("1") for r in rows]
    return od, [n - 1 - d for d in od]

data = parse(sys.argv[1])
by_n = collections.defaultdict(list)
for n, rows in data: by_n[n].append(rows)
for n in sorted(by_n):
    Ts = by_n[n]
    images = collections.Counter(); stat = collections.Counter(); per_T = []
    for rows in Ts:
        hits = []
        for v in range(n):
            keep = [u for u in range(n) if u != v]
            s = sub(rows, keep)
            hc = hcycles(s)
            if hc != 1: continue
            od, idg = degs(s)
            e1 = od.count(1); i1 = idg.count(1)
            hits.append((v, e1, i1, canon_uh(s), bin(rows[v]).count("1")))
        per_T.append(hits)
        stat[tuple(sorted((h[1], h[2]) for h in hits))] += 1
    # candidate canonical rules: deletion giving UH with two sole exits; with two sole entries
    for label, sel in [("two sole exits", lambda h: h[1] == 2), ("two sole entries", lambda h: h[2] == 2)]:
        img = []
        ok = True
        for hits in per_T:
            hs = [h for h in hits if sel(h)]
            cs = set(h[3] for h in hs)
            if len(cs) != 1: ok = False
            img.extend(cs)
        print(f"n={n}: #T={len(Ts)}  rule '{label}': each T has exactly one image class: {ok};  distinct images: {len(set(img))}")
    print(f"   UH-deletion profiles (#sole exits,#sole entries) per T: {dict(stat)}")
