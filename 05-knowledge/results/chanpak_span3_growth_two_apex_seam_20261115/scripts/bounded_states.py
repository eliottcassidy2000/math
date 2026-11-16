#!/usr/bin/env python3
"""bounded_states.py: breadth-first search over the exact transfer states x_W of literal span-3 words with H(W) <= M.
If the search closes (no new state with H <= M), every word with H(W) <= M has a state in the finite set, and the
envelope env[H] = max H(Wl) over words with H(W) = H is exact for all H <= M (not only bounded lengths).
Prints the number of states, the envelope and the maximum ratio, with the two-fragment split."""
import sys, numpy as np
M = int(sys.argv[1])
A = np.load('A.npy'); h = np.load('h.npy'); v0 = np.load('v0.npy'); N = len(h)
nfrag = list(map(int, open('frag.txt').read().split()))
hA = [h @ A[b] for b in range(8)]
two = [np.array([hA[b][s] if nfrag[s] == 2 else 0 for s in range(N)]) for b in range(8)]
start = tuple(int(x) for x in v0)
seen = {start}; frontier = [start]; depth = 0
env = {}; envZ = {}; best = (0, None, None); total_words_seen = 0
while frontier:
    depth += 1; nf = []
    for st in frontier:
        x = np.array(st, dtype=np.int64); Hx = int(h @ x)
        for b in range(8):
            Hy = int(hA[b] @ x); Z = int(two[b] @ x)
            if Hy > env.get(Hx, 0): env[Hx] = Hy
            if Z > envZ.get(Hx, 0): envZ[Hx] = Z
            r = Hy / Hx
            if r > best[0]: best = (r, Hx, Hy)
            if Hy <= M:
                y = tuple(int(v) for v in (A[b] @ x))
                if y not in seen: seen.add(y); nf.append(y)
    frontier = nf
    if depth % 5 == 0 or not nf: print(f"depth {depth}: states with H <= {M}: {len(seen)}, frontier {len(nf)}", flush=True)
    if len(seen) > 3_000_000: print("cap reached"); break
print(f"closed: {not frontier}; states with H <= {M}: {len(seen)}")
print("max ratio over all words with H(W) <= M:", best)
print("exact envelope env[H] (max H(Wl) over ALL words with H(W) = H), and the two-fragment part Z, for H <= M:")
for H in sorted(env):
    if H <= M: print(f"  H={H} env={env[H]} ratio={env[H]/H:.3f} maxZ={envZ[H]} Z/H={envZ[H]/H:.3f}")
