#!/usr/bin/env python3
"""analyze_mg.py n  -- exact native-move distances on isomorphism classes versus the dipole potentials.
Reads mg{n}.txt from movegraph.c.  Checks: blue moves are loops on classes; distance to the non-strong
locus = lambda; distance to the singleton locus vs delta0; distance to the transitive class vs L and the
other potentials; black-only distances; lists the exactness failures with their symmetry data."""
import sys, collections, itertools
import numpy as np
n = int(sys.argv[1]); V = {}; E = collections.defaultdict(list); Eb = collections.defaultdict(list)
for line in open(f'mg{n}.txt'):
    t = line.split()
    if t[0] == 'V':
        idx = int(t[1]); key = int(t[2], 16); sc = list(map(int, t[3].split(',')))
        V[idx] = dict(key=key, s=sc, strong=int(t[4]), sc=int(t[5]), aut=int(t[6]), H=int(t[7]))
    else:
        a, b, blue, c = map(int, t[1:]); E[a].append((b, blue, c))
N = len(V)
def adj_from_key(key):
    A = np.zeros((n, n), dtype=int); bp = 0
    for i in range(n):
        for j in range(i+1, n):
            if key >> bp & 1: A[i, j] = 1
            else: A[j, i] = 1
            bp += 1
    return A
# ---- blue loops
blue_nonloop = sum(c for a in E for (b, blue, c) in E[a] if blue and b != a)
blue_total = sum(c for a in E for (b, blue, c) in E[a] if blue)
print(f'n={n} classes={N} blue moves={blue_total} blue non-loops={blue_nonloop}')
# ---- potentials
pot = {}
for idx, v in V.items():
    s = sorted(v['s']); A = adj_from_key(v['key'])
    L = sum(abs(s[i] - i) for i in range(n)) // 2
    d0 = min(min(x, n-1-x) for x in v['s'])
    # lambda = min over nonempty proper X of e(X -> V\X)
    lam = 10**9
    for mask in range(1, (1 << n) - 1):
        X = [i for i in range(n) if mask >> i & 1]; Y = [i for i in range(n) if not mask >> i & 1]
        f = int(A[np.ix_(X, Y)].sum()); lam = min(lam, f)
    K = A - A.T; ev = np.linalg.eigvalsh(1j * K.astype(complex)); nu = float(np.abs(ev).sum()); s1 = float(np.abs(ev).max())
    pot[idx] = dict(L=L, d0=d0, lam=lam, nu=nu, s1=s1, maxd=max(abs(2*x-(n-1)) for x in v['s'])/2)
# ---- BFS distances from target sets (graph is undirected on classes: the move is an involution)
def bfs(targets, black_only=False):
    dist = {t: 0 for t in targets}; q = collections.deque(targets)
    while q:
        a = q.popleft()
        for (b, blue, c) in E[a]:
            if black_only and blue: continue
            if b not in dist: dist[b] = dist[a] + 1; q.append(b)
    return dist
TT = [i for i, v in V.items() if sorted(v['s']) == list(range(n))]
NS = [i for i, v in V.items() if not v['strong']]
SING = [i for i, v in V.items() if 0 in v['s'] or (n-1) in v['s']]
dTT = bfs(TT); dNS = bfs(NS); dSING = bfs(SING); dTTb = bfs(TT, True); dNSb = bfs(NS, True); dSINGb = bfs(SING, True)
def report(name, dist, potname, fn):
    fails = []; unreachable = 0
    for i in V:
        if i not in dist: unreachable += 1; continue
        b = fn(i)
        if dist[i] < b: print('  BOUND VIOLATED', name, i, dist[i], b)
        if dist[i] > b: fails.append(i)
    print(f'  {name}: potential {potname}: exact on {N-len(fails)-unreachable}/{N}, failures {len(fails)}, unreachable {unreachable}')
    return fails
print('both colours:')
fTT = report('dist to TT', dTT, 'L', lambda i: pot[i]['L'])
fNS = report('dist to non-strong', dNS, 'lambda', lambda i: pot[i]['lam'])
fS = report('dist to singleton', dSING, 'delta0', lambda i: pot[i]['d0'])
print('black only:')
fTTb = report('dist to TT', dTTb, 'L', lambda i: pot[i]['L'])
fNSb = report('dist to non-strong', dNSb, 'lambda', lambda i: pot[i]['lam'])
fSb = report('dist to singleton', dSINGb, 'delta0', lambda i: pot[i]['d0'])
def show(fails, dist, potname):
    rows = collections.Counter()
    for i in fails:
        v = V[i]; p = pot[i]
        rows[(dist[i], p[potname], v['sc'], v['aut'], v['strong'])] += 1
    for (d, b, sc, aut, st), c in sorted(rows.items()):
        print(f'    dist={d} bound={b} sc={sc} aut={aut} strong={st}: {c} classes')
print('failures to TT (both colours) by (dist, L, sc, aut, strong):'); show(fTT, dTT, 'L')
print('failures to singleton (both colours):'); show(fS, dSING, 'd0')
print('failures to non-strong, black only:'); show(fNSb, dNSb, 'lam')
print('failures to singleton, black only:'); show(fSb, dSINGb, 'd0')
# how far does TT distance exceed L overall; maximum distances
print('max dist to TT', max(dTT.values()), ' classes at max:', [ (i, V[i]['s'], pot[i]['L']) for i in V if dTT.get(i) == max(dTT.values())][:6])
# self-converse classes among failures and whether every L-tight move is blue
if n <= 8:
    for i in fTT[:12]:
        v = V[i]; p = pot[i]
        tight = [(b, blue, c) for (b, blue, c) in E[i] if pot[b]['L'] == p['L'] - 1]
        descend = [(b, blue, c) for (b, blue, c) in E[i] if dTT[b] == dTT[i] - 1]
        print(f'  class {i} scores={v["s"]} sc={v["sc"]} aut={v["aut"]} H={v["H"]} dist={dTT[i]} L={p["L"]} L-tight moves: {sum(c for _,_,c in tight)} (blue {sum(c for _,bl,c in tight if bl)}) geodesic moves: {sum(c for _,_,c in descend)}')
