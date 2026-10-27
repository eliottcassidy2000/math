#!/usr/bin/env python3
"""singleton_failures.py n -- the classes whose native distance to the singleton locus exceeds delta0,
tested against the condensation criterion: a legal word starts in the initial strong component and
ends in the terminal one, so a vertex v with s(v) = delta0 can only gain by being a start (v initial)
and a vertex with n-1-s(v) = delta0 only by being an end (v terminal).  Also: L-local-minima census."""
import sys, collections
import numpy as np
n = int(sys.argv[1]); V = {}; E = collections.defaultdict(list)
for line in open(f'mg{n}.txt'):
    t = line.split()
    if t[0] == 'V': V[int(t[1])] = dict(key=int(t[2], 16), s=list(map(int, t[3].split(','))), strong=int(t[4]), sc=int(t[5]), aut=int(t[6]))
    else: a, b, blue, c = map(int, t[1:]); E[a].append((b, blue))
def adj_from_key(key):
    A = np.zeros((n, n), dtype=int); bp = 0
    for i in range(n):
        for j in range(i+1, n):
            if key >> bp & 1: A[i, j] = 1
            else: A[j, i] = 1
            bp += 1
    return A
def components(A):
    """strong components in condensation order (first dominates all later ones)"""
    rem = set(range(n)); comps = []
    while rem:
        # the initial component of the remaining tournament: vertices not beaten by anyone outside their own closure
        # compute reachability closure from each vertex within rem
        best = None
        for v in rem:
            reach = {v}; stack = [v]
            while stack:
                x = stack.pop()
                for y in rem:
                    if y not in reach and A[x, y]: reach.add(y); stack.append(y)
            if reach == rem:  # v reaches everything: v is in the initial component
                # initial component = vertices mutually reachable with v
                comp = {u for u in rem if u == v or (A[u, v] and any(True for _ in [0]) )}
                best = v; break
        # mutual reachability
        v = best
        back = {v}; stack = [v]
        while stack:
            x = stack.pop()
            for y in rem:
                if y not in back and A[y, x]: back.add(y); stack.append(y)
        comp = back  # vertices that reach v; since v reaches all, these are mutually reachable with v
        comps.append(sorted(comp)); rem -= comp
    return comps
SING = [i for i, v in V.items() if 0 in v['s'] or (n-1) in v['s']]
dist = {t: 0 for t in SING}; q = collections.deque(SING)
while q:
    a = q.popleft()
    for (b, blue) in E[a]:
        if b not in dist: dist[b] = dist[a] + 1; q.append(b)
explained = unexplained = 0; crit_but_exact = 0
for i, v in V.items():
    d0 = min(min(x, n-1-x) for x in v['s'])
    A = adj_from_key(v['key']); comps = components(A); init, term = set(comps[0]), set(comps[-1])
    sc_ = [int(x) for x in A.sum(axis=1)]  # scores in the canonical labelling (the V-line scores use the input labelling)
    eligible = [u for u in range(n) if (sc_[u] == d0 and u in init) or (n-1-sc_[u] == d0 and u in term)]
    if dist[i] > d0:
        ok = (len(eligible) == 0)
        explained += ok; unexplained += (not ok)
        print(f'  class {i} scores={v["s"]} comps={[len(c) for c in comps]} sc={v["sc"]} aut={v["aut"]} dist={dist[i]} delta0={d0} eligible={eligible} {"EXPLAINED" if ok else "UNEXPLAINED"}')
    elif len(eligible) == 0 and d0 > 0:
        crit_but_exact += 1; print('  criterion fires but distance = delta0 ?!', i, v['s'])
print(f'n={n}: singleton failures {explained+unexplained}, explained by the condensation criterion {explained}, unexplained {unexplained}, criterion-fires-but-exact {crit_but_exact}')
# L local minima
def L_pot(s):
    s = sorted(s); return sum(abs(s[i] - i) for i in range(n)) // 2
Lv = {i: L_pot(v['s']) for i, v in V.items()}
mins = [i for i in V if Lv[i] > 0 and all(Lv[b] >= Lv[i] for (b, _) in E[i])]
rows = collections.Counter((V[i]['sc'], V[i]['aut'], V[i]['strong']) for i in mins)
print(f'n={n}: L-local-minima classes {len(mins)} of {len(V)}; by (sc, aut, strong):', dict(rows))
