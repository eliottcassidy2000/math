#!/usr/bin/env python3
"""(SS) at k = 1 cannot be realised by local moves (HYP-9280; session opus-2026-10-14-S36, post-audit; seeded random sample).

For a vertex w of in-degree d >= 3, X = {Hamiltonian paths with N <= 1}, Y = {N >= d - 1} (|Y| >= |X| in every tournament
tested).  Edges join P in X to the paths of Y obtained from P by relocating one vertex (test 1), or one vertex or a block
of two consecutive vertices (test 2).  A maximum bipartite matching that misses some P in X shows that no injection X -> Y
built from such moves exists for that (T, w).  Random strong tournaments, n = 6, 7, 8, random.seed(1).
"""
import random, sys, itertools

def ham_paths(n, adj):
    out = []
    def rec(path, used):
        if len(path) == n:
            out.append(tuple(path)); return
        v = path[-1]
        for x in range(n):
            if not used >> x & 1 and adj[v] >> x & 1:
                path.append(x); rec(path, used | 1 << x); path.pop()
    for s in range(n):
        rec([s], 1 << s)
    return out

def is_path(adj, seq):
    return all(adj[seq[i]] >> seq[i + 1] & 1 for i in range(len(seq) - 1))

def strong(n, adj):
    def reach(fwd):
        seen = {0}; st = [0]
        while st:
            v = st.pop()
            for x in range(n):
                if x not in seen and ((adj[v] >> x & 1) if fwd else (adj[x] >> v & 1)):
                    seen.add(x); st.append(x)
        return len(seen) == n
    return reach(True) and reach(False)

def moves(seq, block):
    n = len(seq)
    res = set()
    for L in range(1, block + 1):
        for i in range(n - L + 1):
            blk = seq[i:i + L]; rest = seq[:i] + seq[i + L:]
            for j in range(len(rest) + 1):
                if j == i:
                    continue
                res.add(rest[:j] + blk + rest[j:])
    return res

def max_matching(X, edges):
    match = {}
    def aug(x, seen):
        for y in edges[x]:
            if y in seen:
                continue
            seen.add(y)
            if y not in match or aug(match[y], seen):
                match[y] = x
                return True
        return False
    return sum(1 for x in X if aug(x, set()))

random.seed(1)
for n, trials in ((6, 400), (7, 300), (8, 120)):
    stats = {1: [0, 0], 2: [0, 0]}
    done = 0
    while done < trials:
        adj = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
        if not strong(n, adj):
            continue
        P = ham_paths(n, adj)
        Pset = set(P)
        for w in range(n):
            U = {u for u in range(n) if u != w and adj[u] >> w & 1}
            d = len(U)
            if d < 3:
                continue
            Nof = {p: sum(1 for u in p[:p.index(w)] if u in U) for p in P}
            X = [p for p in P if Nof[p] <= 1]
            Yset = {p for p in P if Nof[p] >= d - 1}
            for block in (1, 2):
                edges = {x: [y for y in moves(x, block) if y in Yset] for x in X}
                m = max_matching(X, edges)
                stats[block][0] += 1
                if m < len(X):
                    stats[block][1] += 1
        done += 1
    print(f'n={n}: (T,w) pairs with d>=3: {stats[1][0]}; single-vertex relocation matching fails in {stats[1][1]}; '
          f'vertex-or-2-block relocation fails in {stats[2][1]}', flush=True)
