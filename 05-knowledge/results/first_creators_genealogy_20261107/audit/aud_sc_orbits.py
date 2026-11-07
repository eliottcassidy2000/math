"""Auditor's check of the orbit step in the proof of HYP-9512 (D): for n = 5, 7 the orbits of the centralizer
G = C(sigma) (order 2^m m!) on X (sigma-anti-invariant tournaments) are in bijection with the self-converse classes,
and they split evenly by H mod 4.  Reuses build/hpaths/keys from aud_sc_mod4.py (same folder)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = [sys.argv[0], '0']
import importlib.util
spec = importlib.util.spec_from_file_location('m4', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'aud_sc_mod4.py'))
src = open(spec.origin).read().split('samples = int')[0]          # definitions only
ns = {}; exec(src, ns)
build, hpaths, keys = ns['build'], ns['hpaths'], ns['keys']

def code(adj, perm):
    n = len(adj); return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i + 1, n))

def canon(adj):
    n = len(adj); sc = [sum(r) for r in adj]
    cells = sorted(set(sc)); groups = [[v for v in range(n) if sc[v] == s] for s in cells]
    best = None
    for parts in itertools.product(*[itertools.permutations(g) for g in groups]):
        perm = [v for p in parts for v in p]
        c = code(adj, perm)
        if best is None or c < best: best = c
    return best

for m in (2, 3):
    n = 2 * m + 1; ks = keys(m)
    X = []
    for v in itertools.product((0, 1), repeat=len(ks)):
        adj = build(m, dict(zip(ks, v)))
        X.append(tuple(tuple(r) for r in adj))
    Xset = set(X)
    # centralizer of sigma: permute pairs, optionally swap inside each pair, fix 0
    G = []
    for pp in itertools.permutations(range(1, m + 1)):
        for sw in itertools.product((0, 1), repeat=m):
            g = [0] * n
            for i in range(1, m + 1):
                j = pp[i - 1]; a, b = 2 * i - 1, 2 * i; c, d = 2 * j - 1, 2 * j
                g[a], g[b] = (c, d) if not sw[i - 1] else (d, c)
            G.append(g)
    seen = set(); orbits = []
    for T in X:
        if T in seen: continue
        orb = set()
        for g in G:
            gi = [0] * n
            for u in range(n): gi[g[u]] = u
            T2 = tuple(tuple(T[gi[u]][gi[w]] for w in range(n)) for u in range(n))
            assert T2 in Xset
            orb.add(T2)
        seen |= orb; orbits.append(T)
    classes = set(canon([list(r) for r in T]) for T in X)
    res = {1: 0, 3: 0}
    for T in orbits: res[hpaths([list(r) for r in T]) % 4] += 1
    print('n=%d: |X| = %d, |G| = %d, G-orbits on X = %d, isomorphism classes among X = %d, orbits by H mod 4: %s'
          % (n, len(X), len(G), len(orbits), len(classes), res))
