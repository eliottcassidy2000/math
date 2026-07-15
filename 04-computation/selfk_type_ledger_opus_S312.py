# opus-2026-07-15-S312 -- HYP-6895: the Burnside type ledger.
# Both sides of weighted-black-qf = SC as sums over permutation cycle types:
#   anti side: (1/n!) Sum_pi #{Y: pi Y = rev Y}                      = SC
#   qf side:   (1/n!) Sum_pi #{(X,p): p dir Ham path, pi X = X^p, black} = wbqf
# The witnesses have ODD order; the anti perms have specific involutive-ish types.
# Ledger: per cycle type, both masses, at n = 5, 6, 7.
import sys, itertools, time
from collections import defaultdict
sys.path.insert(0, '04-computation')
from smith_diagram_of_the_metagraph_opus_S307 import build

def cycle_type(pi, n):
    seen = [False]*n; ct = []
    for v in range(n):
        if not seen[v]:
            l = 0; w = v
            while not seen[w]:
                seen[w] = True; w = pi[w]; l += 1
            ct.append(l)
    return tuple(sorted(ct))

for n in (5, 6, 7):
    t0 = time.time()
    m = n*(n-1)//2 - (n-1)
    B = build(n)
    cls_of, H_of, rcls = B['cls_of'], B['H_of'], B['rcls']
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
    tidx = {t: i for i, t in enumerate(tiles)}
    sig = [tidx[(n+1-y, n+1-x)] for (x, y) in tiles]
    FULL = (1 << m) - 1
    def sig_t(t):
        s = 0
        for i in range(m):
            if (t >> i) & 1: s |= 1 << sig[i]
        return s

    # arc adjacency of tiling t (vertices 1..n; arc matrix)
    def adj_of(t):
        adj = [[0]*(n+1) for _ in range(n+1)]
        for k in range(2, n+1): adj[k][k-1] = 1
        for i, (x, y) in enumerate(tiles):
            if (t >> i) & 1: adj[x][y] = 1
            else: adj[y][x] = 1
        return adj

    perms = list(itertools.permutations(range(1, n+1)))

    # qf side: witnesses of black qf tilings, weighted count per type
    qf_mass = defaultdict(int)     # type -> #(t, pi) black qf pairs (t a tiling)
    nqf_black = 0
    for t in range(1 << m):
        if cls_of[t] != cls_of[t ^ FULL]: continue
        if sig_t(t) == t: continue          # black only
        nqf_black += 1
        A, Bm = adj_of(t), adj_of(t ^ FULL)
        for pi in perms:
            # pi: V -> V; check pi(A) == Bm:  Bm[pi[u]][pi[v]] == A[u][v]
            ok = True
            for u in range(1, n+1):
                pu = pi[u-1]
                Au = A[u]
                for v in range(1, n+1):
                    if u != v and Bm[pu][pi[v-1]] != Au[v]: ok = False; break
                if not ok: break
            if ok:
                qf_mass[cycle_type([p-1 for p in pi], n)] += 1

    # anti side: for each SC class rep, anti-automorphisms; mass per type
    # anti-pairs total = Sum_{SC classes} |Aut| * (n!/|Aut|)/... per class:
    # labeled copies n!/Aut, each with Aut anti-autos -> n! per class; type profile
    # from ONE labeled rep: its anti-auto set, then conjugation spreads types...
    # SAFER: per class, anti-autos of the rep; each labeled copy is conjugate ->
    # same TYPE multiset. mass(type) per class = (n!/Aut) * #anti-autos-of-rep of
    # that type.
    rep = {}
    fiber = defaultdict(int)
    for t in range(1 << m):
        c = cls_of[t]
        fiber[c] += 1
        if c not in rep: rep[c] = t
    anti_mass = defaultdict(int)
    SCcls = [c for c in fiber if rcls[c] == c]
    fact = 1
    for i in range(2, n+1): fact *= i
    for c in SCcls:
        A = adj_of(rep[c])
        R = [list(col) for col in zip(*A)]  # rev: R[a][b] = A[b][a]
        aut = H_of[c] // fiber[c]
        antis = []
        for pi in perms:
            ok = True
            for u in range(1, n+1):
                pu = pi[u-1]
                for v in range(1, n+1):
                    if u != v and R[pu][pi[v-1]] != A[u][v]: ok = False; break
                if not ok: break
            if ok: antis.append(pi)
        assert len(antis) == aut, (c, len(antis), aut)
        for pi in antis:
            anti_mass[cycle_type([p-1 for p in pi], n)] += fact // aut

    print(f"\n===== n={n} (black qf tilings: {nqf_black}, SC: {len(SCcls)}) "
          f"[{time.time()-t0:.0f}s]")
    print(f"qf-side (t,pi) mass by witness type (total = weighted = "
          f"{sum(qf_mass.values())}):")
    for ct, v in sorted(qf_mass.items()): print(f"    {ct}: {v}")
    print(f"anti-side mass by type (total/n! = SC = "
          f"{sum(anti_mass.values())//fact}):")
    for ct, v in sorted(anti_mass.items()): print(f"    {ct}: {v} = n!*{v/fact:g}")
