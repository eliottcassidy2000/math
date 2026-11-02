#!/usr/bin/env python3
"""prime_path_posa.py (klein-2026-11-01) -- search a Hamiltonian path through the primes below N in the square-sum graph
(consecutive sums perfect squares) by Posa rotation-extension with Warnsdorff ordering, the leaf prime (if any) fixed as
the start.  The result is verified (distinct primes, every consecutive sum a square) and written to
data/prime_path_N<N>.txt.  usage: python3 prime_path_posa.py N [seed=1] [max_rotations=2000000] [time_limit_s=1500]"""
import math, random, sys, time, os
N = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
max_rot = int(sys.argv[3]) if len(sys.argv) > 3 else 2_000_000
tlimit = float(sys.argv[4]) if len(sys.argv) > 4 else 1500.0
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
n = len(P); degs = {p: len(adj[p]) for p in P}
leaves = [p for p in P if degs[p] == 1]
A = sum(1 for p in P if p == 2 or p % 4 == 1); B = n - A
print('N = %d: %d primes, %d edges, min degree %d, leaves %s, |A| - |B| = %d' % (N, n, sum(degs.values()) // 2, min(degs.values()), leaves, A - B), flush=True)
assert len(leaves) <= 2 and abs(A - B) <= 1 and min(degs.values()) >= 1
rng = random.Random(seed)
start = leaves[0] if leaves else min(P, key=lambda p: degs[p])
path = [start]; pos = {start: 0}
unused_deg = dict(degs)           # number of unused neighbours of each vertex
def use(w):
    pos[w] = len(path); path.append(w)
    for x in adj[w]: unused_deg[x] -= 1
t0 = time.time(); rotations = 0; best_len = 1
while len(path) < n:
    v = path[-1]
    cand = [w for w in adj[v] if w not in pos]
    if cand:
        w = min(cand, key=lambda w: (unused_deg[w], rng.random()))
        use(w)
        if len(path) > best_len: best_len = len(path)
    else:
        rotations += 1
        if rotations > max_rot or time.time() - t0 > tlimit:
            print('gave up: rotations %d, best length %d of %d, %.1f s' % (rotations, best_len, n, time.time() - t0)); sys.exit(1)
        nbrs = [u for u in adj[v] if pos[u] < len(path) - 2]
        if not nbrs:
            print('stuck with no rotation at length %d' % len(path)); sys.exit(1)
        sample = nbrs if len(nbrs) <= 24 else rng.sample(nbrs, 24)
        u = max(sample, key=lambda u: (unused_deg[path[pos[u] + 1]], rng.random()))
        i = pos[u] + 1
        path[i:] = path[i:][::-1]
        for j in range(i, len(path)): pos[path[j]] = j
    if rotations and rotations % 20000 == 0 and cand == []:
        print('  ... rotations %d, length %d, %.0f s' % (rotations, len(path), time.time() - t0), flush=True)
assert len(path) == n and len(set(path)) == n and all(math.isqrt(a + b) ** 2 == a + b for a, b in zip(path, path[1:]))
elapsed = time.time() - t0
print('HAMILTONIAN PATH FOUND: %d primes, %d rotations, %.1f s; starts %s ends %s' % (n, rotations, elapsed, path[:6], path[-6:]))
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data'), exist_ok=True)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'prime_path_N%d.txt' % N)
open(out, 'w').write(' '.join(map(str, path)) + '\n')
print('written', os.path.relpath(out))
