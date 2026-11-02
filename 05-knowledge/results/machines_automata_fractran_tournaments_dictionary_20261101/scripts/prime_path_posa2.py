#!/usr/bin/env python3
"""prime_path_posa2.py (klein-2026-11-01) -- Hamiltonian path through the primes below N in the square-sum graph, second
solver: (i) forced chains: a degree-1 prime is an endpoint and a degree-2 prime forces both its edges, so maximal forced
chains are contracted and kept contiguous; (ii) Posa rotation-extension with Warnsdorff extension (fewest unused
neighbours first), rotations only across non-forced edges, rotation targets scored by the scarcity of the unvisited
primes they expose (sum of 1/deg over unvisited neighbours of the new endpoint); (iii) restarts.  Output verified and
written to data/prime_path_N<N>.txt.  usage: python3 prime_path_posa2.py N [seed=1] [time_limit_s=1500]"""
import math, random, sys, time, os
N = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
tlimit = float(sys.argv[3]) if len(sys.argv) > 3 else 1500.0
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
P = primes_below(N); Ps = set(P); top = P[-1]
adj = {p: [] for p in P}
for p in P:
    s = math.isqrt(p + 2)
    while s * s <= p + top:
        q = s * s - p
        if q > p and q in Ps: adj[p].append(q); adj[q].append(p)
        s += 1
n = len(P); deg = {p: len(adj[p]) for p in P}
hist = {}
for p in P: hist[deg[p]] = hist.get(deg[p], 0) + 1
print('N = %d: %d primes, %d edges; degree histogram (low end): %s' % (N, n, sum(deg.values()) // 2, sorted(hist.items())[:8]), flush=True)
# forced edges: from degree-1 and degree-2 primes
forced = {p: set() for p in P}
for p in P:
    if deg[p] <= 2:
        for q in adj[p]: forced[p].add(q); forced[q].add(p)
bad = [p for p in P if len(forced[p]) > 2]
if bad:
    print('NO HAMILTONIAN PATH: primes with three forced edges:', bad[:5]); sys.exit(2)
# chains: components of the forced graph (paths; a forced cycle would also be fatal)
chain_of = {}; chains = []
for p in P:
    if p in chain_of or not forced[p]: continue
    # walk to one end
    end = p; prev = None
    while True:
        nxt = [q for q in forced[end] if q != prev]
        if not nxt or nxt[0] == p: break
        prev, end = end, nxt[0]
        if end == p: break
    if end == p and len(forced[p]) == 2 and prev is not None:
        print('NO HAMILTONIAN PATH: forced cycle through', p); sys.exit(2)
    # collect from the end
    seq = [end]; prev = None; cur = end
    while True:
        nxt = [q for q in forced[cur] if q != prev]
        if not nxt: break
        prev, cur = cur, nxt[0]; seq.append(cur)
    for v in seq: chain_of[v] = len(chains)
    chains.append(seq)
print('forced chains: %d (lengths %s); longest %d; leaf chains (must be endpoints): %d' % (
    len(chains), sorted({len(c) for c in chains}), max((len(c) for c in chains), default=0),
    sum(1 for c in chains if deg[c[0]] == 1 or deg[c[-1]] == 1)), flush=True)
leaf_chains = [c for c in chains if deg[c[0]] == 1 or deg[c[-1]] == 1]
if len(leaf_chains) > 2: print('NO HAMILTONIAN PATH: more than two leaves'); sys.exit(2)
rng = random.Random(seed)
scarcity = {p: 1.0 / deg[p] for p in P}

def attempt(seed_run, budget_s):
    rng = random.Random(seed_run)
    path = []; pos = {}
    unused = dict(deg)
    def append_vertex(w):
        pos[w] = len(path); path.append(w)
        for x in adj[w]: unused[x] -= 1
    def append_unit(w):
        """append w; if w starts a forced chain, append the whole chain from w's end."""
        if w in chain_of:
            c = chains[chain_of[w]]
            if c[0] == w: seq = c
            elif c[-1] == w: seq = c[::-1]
            else: return False
            for v in seq: append_vertex(v)
        else: append_vertex(w)
        return True
    # start: a leaf chain, oriented leaf first
    if leaf_chains:
        c = leaf_chains[0]; seq = c if deg[c[0]] == 1 else c[::-1]
        for v in seq: append_vertex(v)
    else:
        start = min(P, key=lambda p: (deg[p], rng.random())); append_unit(start)
    t0 = time.time(); rotations = 0; best = len(path)
    def enterable(w):
        return (w not in chain_of) or chains[chain_of[w]][0] == w or chains[chain_of[w]][-1] == w
    while len(path) < n:
        v = path[-1]
        cand = [w for w in adj[v] if w not in pos and enterable(w)]
        if cand:
            w = min(cand, key=lambda w: (unused[w], rng.random()))
            append_unit(w); best = max(best, len(path))
            continue
        rotations += 1
        if time.time() - t0 > budget_s:
            return None, rotations, best
        # rotation at u in N(v) with pos[u] <= len-3 and edge (u, succ u) not forced
        opts = []
        for u in adj[v]:
            i = pos[u]
            if i >= len(path) - 2: continue
            nxt = path[i + 1]
            if nxt in forced[u]: continue
            opts.append((u, nxt))
        if not opts: return None, rotations, best
        sample = opts if len(opts) <= 32 else rng.sample(opts, 32)
        def score(o):
            u, nxt = o
            sc = sum(scarcity[x] for x in adj[nxt] if x not in pos and enterable(x))
            return (sc, rng.random())
        u, nxt = max(sample, key=score)
        i = pos[u] + 1
        path[i:] = path[i:][::-1]
        for j in range(i, len(path)): pos[path[j]] = j
        if rotations % 20000 == 0: print('  ... rotations %d, length %d, %.0f s' % (rotations, len(path), time.time() - t0), flush=True)
    return path, rotations, len(path)

T0 = time.time(); run = 0
while time.time() - T0 < tlimit:
    run += 1
    path, rot, best = attempt(seed * 1000 + run, min(600.0, tlimit - (time.time() - T0)))
    if path:
        assert len(path) == n and len(set(path)) == n and all(math.isqrt(a + b) ** 2 == a + b for a, b in zip(path, path[1:]))
        print('HAMILTONIAN PATH FOUND (run %d): %d primes, %d rotations, %.0f s total; starts %s ends %s' % (run, n, rot, time.time() - T0, path[:5], path[-5:]), flush=True)
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'prime_path_N%d.txt' % N)
        open(out, 'w').write(' '.join(map(str, path)) + '\n'); print('written', os.path.relpath(out)); sys.exit(0)
    print('run %d failed: best length %d of %d after %d rotations' % (run, best, n, rot), flush=True)
print('gave up after %.0f s' % (time.time() - T0)); sys.exit(1)
