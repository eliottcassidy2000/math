#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S84 -- CONJECTURE ATLAS for the merged-metagraph bucket/parity structure + the
Z2xZ2 tiling symmetry (flip = complement-tiling, sigma = complement-fold/transpose = the half-tiling mirror).

Generates a data table across n and tests many small conjectures about:
 - category counts (#pure_blue, #mixed, #pure_black), #SC, #NS-pairs
 - tiling-count parity (SC odd / NS even) and divisibility
 - blue/black line counts (closed forms), self-loops, NS-NS sea onset
 - the <flip, sigma> = Z2xZ2 symmetry orbits on tilings (which tilings share structure)
 - the half-tiling quarter-square shape floor((n-1)^2/4)
 - the structure=constraint UNIQUENESS test (is the colored line-assignment forced?)
"""
import itertools
from math import comb
from collections import defaultdict, Counter

def build(n):
    VERTS = list(range(n, 0, -1))
    TILES = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
    m = len(TILES); tileIdx = {t: i for i, t in enumerate(TILES)}
    TRANS = [tileIdx[(n-y+1, n-x+1)] for (x, y) in TILES]
    vpos = {v: i for i, v in enumerate(VERTS)}
    perms = list(itertools.permutations(range(n)))
    def adj(bits):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1] = 1
        for i, (xL, yL) in enumerate(TILES):
            xi, yi = vpos[xL], vpos[yL]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A
    def canon(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best
    T = []
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        gs = all(TRANS[i] == i or bits[i] == bits[TRANS[i]] for i in range(m))
        flip = mask ^ ((1 << m) - 1)
        tb = [0]*m
        for i in range(m): tb[TRANS[i]] = bits[i]
        sigma = sum(b << k for k, b in enumerate(tb))
        T.append(dict(mask=mask, canon=canon(adj(bits)), gs=gs, flip=flip, sigma=sigma))
    sigs = sorted(set(t['canon'] for t in T)); cidx = {s: i for i, s in enumerate(sigs)}
    for t in T: t['ci'] = cidx[t['canon']]
    m2 = {t['mask']: t for t in T}
    tgt = {}
    for t in T: tgt.setdefault(t['ci'], m2[t['sigma']]['ci'])
    merged = lambda ci: min(ci, tgt[ci])
    return n, m, T, m2, merged, len(sigs), tgt

def stats(n):
    n, m, T, m2, merged, nclasses, tgt = build(n)
    node_t = defaultdict(list)
    for t in T: node_t[merged(t['ci'])].append(t)
    def cat(v):
        g = set(t['gs'] for t in node_t[v])
        return 'BLUE' if g == {True} else ('BLACK' if g == {False} else 'MIXED')
    cats = {v: cat(v) for v in node_t}
    nSC = sum(1 for ci in range(nclasses) if tgt[ci] == ci)
    # lines
    seen = set(); lines = []
    for t in T:
        k = min(t['mask'], t['flip'])
        if k in seen: continue
        seen.add(k)
        a, b = merged(t['ci']), merged(m2[t['flip']]['ci'])
        lines.append((a, b, 'blue' if t['gs'] else 'black'))
    # per-node degrees
    deg = defaultdict(lambda: dict(bo=0, ko=0, bs=0, ks=0, cnt=0))
    for v in node_t: deg[v]['cnt'] = len(node_t[v])
    for a, b, c in lines:
        col = 'b' if c == 'blue' else 'k'
        if a == b: deg[a][col+'s'] += 1
        else: deg[a][col+'o'] += 1; deg[b][col+'o'] += 1
    catcount = Counter(cats.values())
    ngs = sum(1 for t in T if t['gs'])
    nblue = sum(1 for *_, c in lines if c == 'blue')
    # NS-NS sea (black between two BLACK) and pure-black self-loops
    sea = sum(1 for a, b, c in lines if c == 'black' and a != b and cats[a] == 'BLACK' and cats[b] == 'BLACK')
    pbself = sum(deg[v]['ks'] for v in node_t if cats[v] == 'BLACK')
    # Z2xZ2 <flip,sigma> orbits on tilings; sigma-fixed = grid-sym; flip has no fixed pt
    sigfix = sum(1 for t in T if t['sigma'] == t['mask'])
    # tiling-count parity check
    sc_odd = all(deg[v]['cnt'] % 2 == 1 for v in node_t if cats[v] in ('BLUE', 'MIXED'))
    ns_even = all(deg[v]['cnt'] % 2 == 0 for v in node_t if cats[v] == 'BLACK')
    return dict(n=n, m=m, classes=nclasses, nodes=len(node_t), SC=nSC,
                blue=catcount.get('BLUE', 0), mixed=catcount.get('MIXED', 0), black=catcount.get('BLACK', 0),
                gridsym=ngs, bluelines=nblue, blacklines=len(lines)-nblue, sea=sea, pbself=pbself,
                sigfix=sigfix, quarter=((n-1)**2)//4, sc_odd=sc_odd, ns_even=ns_even,
                tcounts={c: sorted(deg[v]['cnt'] for v in node_t if cats[v] == c) for c in ['BLUE','MIXED','BLACK']})

print("DATA TABLE (merged metagraph, n=3..6):\n")
data = {}
for n in range(3, 7):
    s = stats(n); data[n] = s
    print(f"n={n}: classes(A000568)={s['classes']} nodes(V_merged)={s['nodes']} SC={s['SC']} | "
          f"cats [blue={s['blue']} mixed={s['mixed']} black={s['black']}]")
    print(f"     grid-sym tilings={s['gridsym']} sigma-fixed={s['sigfix']} (equal? {s['gridsym']==s['sigfix']}) | "
          f"blue lines={s['bluelines']} black lines={s['blacklines']} sea={s['sea']} pb-selfloops={s['pbself']}")
    print(f"     quarter-square floor((n-1)^2/4)={s['quarter']} vs m=C(n-1,2)={s['m']} | "
          f"PARITY sc_odd={s['sc_odd']} ns_even={s['ns_even']}")
    print(f"     tiling counts: BLUE={s['tcounts']['BLUE']} MIXED={s['tcounts']['MIXED']} BLACK={s['tcounts']['BLACK']}")
    print()

print("SEQUENCES to mine (n=3..6):")
for key in ['nodes','SC','blue','mixed','black','bluelines','blacklines','sea','pbself','gridsym','sigfix']:
    print(f"  {key:12s}: {[data[n][key] for n in range(3,7)]}")
