#!/usr/bin/env python3
"""prime_path_posa4.py (klein-2026-11-01) -- Hamiltonian path through the primes below N in the square-sum graph, fourth
solver = prime_path_posa3 plus a DEPTH-2 ROTATION LOOKAHEAD when no single rotation exposes an unvisited prime (the
endpoint after two rotations is computed in O(1) from the original indices, so all ~deg^2 pairs are scanned), i.e. (forced chains, Posa rotation-extension, scarcity-scored rotations) plus an endgame of
INSERTION moves: an unvisited unit U = (first..last) (a prime, or a forced chain) with path vertices x ~ first,
y ~ last, pos x < pos y, is inserted as  ... x, U, y, y-, ..., x+, y+, ...  whenever y is the endpoint or x+ ~ y+
(a 3-opt move needing one extra edge; the reversed segment keeps forced chains contiguous).  Output verified and
written to data/prime_path_N<N>.txt.  usage: python3 prime_path_posa4.py N [seed=1] [time_limit_s=1500]"""
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
adjs = {p: set(adj[p]) for p in P}
n = len(P); deg = {p: len(adj[p]) for p in P}
forced = {p: set() for p in P}
for p in P:
    if deg[p] <= 2:
        for q in adj[p]: forced[p].add(q); forced[q].add(p)
if any(len(forced[p]) > 2 for p in P): print('NO HAMILTONIAN PATH: a prime with three forced edges'); sys.exit(2)
chain_of = {}; chains = []
for p in P:
    if p in chain_of or not forced[p]: continue
    end = p; prev = None
    while True:
        nxt = [q for q in forced[end] if q != prev]
        if not nxt: break
        prev, end = end, nxt[0]
        if end == p: print('NO HAMILTONIAN PATH: forced cycle'); sys.exit(2)
    seq = [end]; prev = None; cur = end
    while True:
        nxt = [q for q in forced[cur] if q != prev]
        if not nxt: break
        prev, cur = cur, nxt[0]; seq.append(cur)
    for v in seq: chain_of[v] = len(chains)
    chains.append(seq)
leaf_chains = [c for c in chains if deg[c[0]] == 1 or deg[c[-1]] == 1]
if len(leaf_chains) > 2: print('NO HAMILTONIAN PATH: more than two leaves'); sys.exit(2)
print('N = %d: %d primes, %d edges, min degree %d, %d forced chains, %d leaf chain(s)' % (N, n, sum(deg.values()) // 2, min(deg.values()), len(chains), len(leaf_chains)), flush=True)
scarcity = {p: 1.0 / deg[p] for p in P}
def unit_of(w):
    if w in chain_of:
        c = chains[chain_of[w]]
        return c if c[0] == w else (c[::-1] if c[-1] == w else None)
    return [w]

def attempt(seed_run, budget_s):
    rng = random.Random(seed_run)
    path = []; pos = {}; unused = dict(deg); unvisited = set(P)
    def append_vertex(w):
        pos[w] = len(path); path.append(w); unvisited.discard(w)
        for x in adj[w]: unused[x] -= 1
    def reindex(i):
        for j in range(i, len(path)): pos[path[j]] = j
    if leaf_chains:
        c = leaf_chains[0]; seq = c if deg[c[0]] == 1 else c[::-1]
    else:
        seq = unit_of(min(P, key=lambda p: (deg[p], rng.random())))
    for v in seq: append_vertex(v)
    t0 = time.time(); rotations = 0; insertions = 0; last_try = -10**9
    def enterable(w): return unit_of(w) is not None
    def try_insertions():
        nonlocal insertions
        units = []
        seen = set()
        for w in unvisited:
            if w in chain_of:
                cid = chain_of[w]
                if cid in seen: continue
                seen.add(cid); units.append(chains[cid])
            else: units.append([w])
        units.sort(key=lambda u: min(deg[u[0]], deg[u[-1]]))
        L = len(path); done = False
        for u in units:
            for orient in (u, u[::-1]):
                first, last = orient[0], orient[-1]
                xs = [x for x in adj[first] if x in pos]
                ys = [y for y in adj[last] if y in pos]
                ysp = {pos[y]: y for y in ys}
                for x in xs:
                    i = pos[x]
                    if i + 1 < L and path[i + 1] in forced[x]: continue
                    xp = path[i + 1] if i + 1 < L else None
                    if xp is None: continue   # x is the endpoint: plain extension would have caught first ~ x
                    for j, y in ysp.items():
                        if j <= i: continue
                        if j == L - 1 or (path[j + 1] not in forced[y] and path[j + 1] in adjs[xp]):
                            # insert
                            seg = path[i + 1:j + 1][::-1]
                            tail = path[j + 1:]
                            del path[i + 1:]
                            for v in orient: pos[v] = len(path); path.append(v); unvisited.discard(v)
                            for v in orient:
                                for z in adj[v]: unused[z] -= 1
                            path.extend(seg); path.extend(tail); reindex(i + 1)
                            insertions += 1; done = True; break
                    if done: break
                if done: break
            if done: return True
        return False
    while len(path) < n:
        v = path[-1]
        cand = [w for w in adj[v] if w not in pos and enterable(w)]
        if cand:
            w = min(cand, key=lambda w: (unused[w], rng.random()))
            for x in unit_of(w): append_vertex(x)
            continue
        if len(unvisited) <= 400 and rotations - last_try >= 10:
            last_try = rotations
            if try_insertions(): continue
        rotations += 1
        if time.time() - t0 > budget_s: return None, rotations, insertions, len(path)
        opts = []
        for u in adj[v]:
            i = pos[u]
            if i >= len(path) - 2: continue
            nxt = path[i + 1]
            if nxt in forced[u]: continue
            opts.append((u, nxt))
        if not opts: return None, rotations, insertions, len(path)
        def expo(e): return sum(scarcity[x] for x in adj[e] if x not in pos and enterable(x))
        scored = [(expo(o[1]), rng.random(), o) for o in opts]
        best1 = max(scored)
        if best1[0] > 0:
            u, nxt = best1[2]; i = pos[u] + 1
            path[i:] = path[i:][::-1]; reindex(i)
        else:
            # depth-2 lookahead: endpoint after rotating at u (index i) then at u2 (original index j)
            L = len(path); best2 = (0.0, 0.0, None)
            for u, nxt in opts:
                i = pos[u]
                for u2 in adj[nxt]:
                    j = pos.get(u2)
                    if j is None or j == i + 1: continue
                    if j < i:
                        e2 = path[j + 1]
                        if e2 in forced[u2]: continue
                    elif j == i:
                        e2 = v
                        if e2 in forced[u2]: continue
                    elif j >= i + 3:
                        e2 = path[j - 1]
                        if e2 in forced[u2]: continue
                    else: continue
                    sc = expo(e2)
                    if sc > best2[0]: best2 = (sc, rng.random(), (u, u2))
            if best2[2] is not None:
                u, u2 = best2[2]; i = pos[u] + 1
                path[i:] = path[i:][::-1]; reindex(i)
                k = pos[u2] + 1
                path[k:] = path[k:][::-1]; reindex(k)
                rotations += 1
            else:
                u, nxt = rng.choice(opts); i = pos[u] + 1
                path[i:] = path[i:][::-1]; reindex(i)
        if rotations % 20000 == 0: print('  ... rotations %d, insertions %d, length %d, %.0f s' % (rotations, insertions, len(path), time.time() - t0), flush=True)
    return path, rotations, insertions, len(path)

T0 = time.time(); run = 0
while time.time() - T0 < tlimit:
    run += 1
    path, rot, ins, best = attempt(seed * 1000 + run, min(900.0, tlimit - (time.time() - T0)))
    if path:
        assert len(path) == n and len(set(path)) == n and all(math.isqrt(a + b) ** 2 == a + b for a, b in zip(path, path[1:]))
        print('HAMILTONIAN PATH FOUND (run %d): %d primes, %d rotations, %d insertions, %.0f s; starts %s ends %s' % (run, n, rot, ins, time.time() - T0, path[:5], path[-5:]), flush=True)
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'prime_path_N%d.txt' % N)
        open(out, 'w').write(' '.join(map(str, path)) + '\n'); print('written', os.path.relpath(out)); sys.exit(0)
    print('run %d failed: length %d of %d after %d rotations, %d insertions' % (run, best, n, rot, ins), flush=True)
print('gave up after %.0f s' % (time.time() - T0)); sys.exit(1)
