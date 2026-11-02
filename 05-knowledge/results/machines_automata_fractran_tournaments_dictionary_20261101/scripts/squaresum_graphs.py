#!/usr/bin/env python3
"""squaresum_graphs.py -- the square-sum graphs G_n on {1..n} (i ~ j iff i + j is a perfect square) for n = 1..25:
edges, components, degree extremes, Hamiltonian path and cycle existence (exact backtracking), plus the partition of
the 25 primes below 100 into Ogg's supersingular layers (Fi22, B\\Fi22, M\\B) and the two non-supersingular pentads,
and the role of those primes as vertices of G_25.  Standard library only."""
import itertools, sys
sys.setrecursionlimit(10000)
def graph(n):
    sq = {k * k for k in range(2, 20)}
    adj = {i: set() for i in range(1, n + 1)}
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if i + j in sq: adj[i].add(j); adj[j].add(i)
    return adj
def components(adj):
    seen = set(); comps = []
    for v in adj:
        if v in seen: continue
        stack = [v]; comp = set()
        while stack:
            x = stack.pop()
            if x in comp: continue
            comp.add(x); stack.extend(adj[x] - comp)
        seen |= comp; comps.append(sorted(comp))
    return comps
def ham_path(adj, cycle=False):
    n = len(adj); verts = sorted(adj, key=lambda v: len(adj[v]))
    for start in (verts if not cycle else verts[:1]):
        path = [start]; used = {start}
        def rec():
            if len(path) == n:
                return (not cycle) or (path[0] in adj[path[-1]])
            last = path[-1]
            for w in sorted(adj[last] - used, key=lambda v: len(adj[v] - used)):
                path.append(w); used.add(w)
                if rec(): return True
                path.pop(); used.discard(w)
            return False
        if rec(): return list(path)
    return None
print('n | edges | comps (sizes) | min deg | isolated | Ham path | Ham cycle')
rows = {}
for n in range(1, 26):
    adj = graph(n); e = sum(len(s) for s in adj.values()) // 2
    comps = components(adj); sizes = sorted((len(c) for c in comps), reverse=True)
    iso = sum(1 for v in adj if not adj[v]); mind = min(len(s) for s in adj.values())
    hp = ham_path(adj) if (len(comps) == 1 and mind >= 1) else None
    hc = ham_path(adj, cycle=True) if (hp and mind >= 2) else None
    rows[n] = (e, sizes, mind, iso, hp is not None, hc is not None, hp)
    print('%2d | %3d | %s | %d | %d | %s | %s' % (n, e, sizes, mind, iso, 'yes' if hp else 'no', 'yes' if hc else 'no'))
for n in (15, 16, 17, 23, 25):
    print('Hamiltonian path for n=%d: %s' % (n, rows[n][6]))
# the 25 primes below 100 and Ogg's layers
primes = [p for p in range(2, 100) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
Fi22 = [2, 3, 5, 7, 11, 13]; B_minus = [17, 19, 23, 31, 47]; M_minus = [29, 41, 59, 71]
supers = Fi22 + B_minus + M_minus; nons = [p for p in primes if p not in supers]
print('primes below 100:', len(primes), primes)
print('supersingular (divide |M|):', sorted(supers), 'count', len(supers), '; non-supersingular:', nons, 'count', len(nons))
print('owner pentads:', [37, 43, 53, 61, 67], [73, 79, 83, 89, 97], 'equal the two size-halves of the non-supersingular primes:', nons[:5] == [37, 43, 53, 61, 67] and nons[5:] == [73, 79, 83, 89, 97])
adj25 = graph(25)
print('primes as vertices of G_25 with their degrees:', {p: len(adj25[p]) for p in primes if p <= 25})
print('square sums among primes <= 25 (edges of G_25 between primes):', sorted((p, q) for p in primes if p <= 25 for q in primes if p < q <= 25 and q in adj25[p]))
