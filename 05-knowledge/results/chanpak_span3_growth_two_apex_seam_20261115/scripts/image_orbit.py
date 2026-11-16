#!/usr/bin/env python3
"""image_orbit.py: is the set of image subspaces {col(A_W)} of literal span-3 prefix products finite?
Works in the 119-state representation modulo the prime P with exact row-echelon canonical forms of column spaces.
BFS: start from the full space (identity); apply each letter to a basis; canonicalise; stop when closed or at CAP."""
import sys, numpy as np
A = np.load('A.npy'); N = A.shape[1]; P = 2147483629; CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
def rref(M):
    """row-reduced echelon basis (tuple of tuples) of the row space of M mod P; rows = basis vectors"""
    M = [[int(x) % P for x in row] for row in M]; rk = 0; cols = N; out = []
    for c in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], P - 2, P); M[rk] = [x * inv % P for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % P for x, y in zip(M[i], M[rk])]
        rk += 1
        if rk == len(M): break
    return tuple(tuple(r) for r in M[:rk])
# column space of A_W: start with v0-reachable? use the full identity (all states) as the start image
start = rref(np.identity(N, dtype=np.int64).tolist())
seen = {start: 0}; frontier = [start]; depth = 0; from collections import Counter
while frontier and len(seen) < CAP:
    depth += 1; nf = []
    for S in frontier:
        B = np.array(S, dtype=object)  # basis rows (each a vector in the 119-space): image = span of columns of A_W = A_W * span
        for b in range(8):
            # image of A_b applied to span(S): vectors A_b s for s in basis
            img = [[int(sum(int(A[b][i][j]) * int(s[j]) for j in range(N)) % P) for i in range(N)] for s in S]
            key = rref(img)
            if key and key not in seen: seen[key] = depth; nf.append(key)
    frontier = nf
    ranks = Counter(len(k) for k in seen)
    print(f"depth {depth}: distinct image subspaces {len(seen)} (new {len(nf)}); ranks {sorted(ranks.items())[:12]}", flush=True)
    if not nf: print("CLOSED: the image orbit is finite"); break
