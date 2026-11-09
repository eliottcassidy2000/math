"""AUDIT D1: independent census of the CONNECTED graphs with d in {3, 5, 7, 11, 13, 23, 25, 37, 47}, from aw_dump
(exhaustive over all graphs on <= 9 vertices; by Ore's bound d >= 2^(n/2 + 1) - 1 every connected realiser of
d <= 62 has <= 9 vertices).  Canonical forms by brute force over all permutations; parity by brute force over all
automorphisms (a third method, independent of aw_verify and of the census code).
usage: python3 aw_smalld.py aw_dump_smalld.out"""
import sys, itertools
from collections import defaultdict
def connected(n, E):
    adj = defaultdict(set)
    for a, b in E: adj[a].add(b); adj[b].add(a)
    seen = {0}; st = [0]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen: seen.add(w); st.append(w)
    return len(seen) == n
def canon(n, E):
    S = set(E) | set((b, a) for a, b in E); best = None
    for p in itertools.permutations(range(n)):
        code = tuple(1 if (p[u], p[v]) in S else 0 for u in range(n) for v in range(u + 1, n))
        if best is None or code < best: best = code
    return best
def parity_odd(n, E):
    S = set(E) | set((b, a) for a, b in E)
    for p in itertools.permutations(range(n)):
        if all(((p[a], p[b]) in S) for a, b in E):
            if sum(1 for a, b in E if p[min(a, b)] > p[max(a, b)]) % 2: return True
    return False
classes = defaultdict(dict)    # d -> canon -> (n, E)
for line in open(sys.argv[1]):
    p = line.split(); n, d = int(p[0]), int(p[1])
    E = [tuple(map(int, e.split('-'))) for e in p[2:]]
    if n == 1 or not connected(n, E): continue
    c = canon(n, E)
    classes[d].setdefault((n, c), (n, E))
for d in sorted(classes):
    rows = []
    for (n, c), (_, E) in sorted(classes[d].items()):
        deg = [0] * n
        for a, b in E: deg[a] += 1; deg[b] += 1
        rows.append((n, len(E), 'ODD' if parity_odd(n, E) else 'EVEN', 'pendant' if 1 in deg else 'no-pendant', sorted(deg)))
    ne = sum(1 for r in rows if r[2] == 'EVEN')
    print(f"d={d}: {len(rows)} connected classes, {ne} even; even with a pendant: "
          f"{sum(1 for r in rows if r[2] == 'EVEN' and r[3] == 'pendant')}")
    for r in rows: print("   n=%d m=%d %s %s degrees=%s" % r)
