#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, session 2 part 2 — THE CHANG SHADOW: tower-Ramsey numbers at m=1.
mac-mini-2026-06-09-S2  (THM-460, HYP-2366)

Ambient: vectors [C]^s, plain lex (index 0 most significant) — a truncation of
omega^omega in its omega^s presentation.  Split(u,v) = first differing index.

BINARY h-GRID (THM-460 B3(ii)): 2^h elements; sorted halves A < B; all cross pairs
split at one common position p; all internal splits at positions > p (less
significant); halves recursively binary (h-1)-grids.  (Tails free.)

BINARY M-TOWER (m=1): B_1 < B_2 < ... < B_M (order-separated), B_j a binary j-grid.
THM-460 B2/B3: every subset of omega^omega of full type contains one for every M;
THM-460 C1/C2: if Q(M; s,C) := "exists triangle-free graph on [C]^s with no
independent binary M-tower" were SAT for ALL (s,C), omega^omega -/-> (omega^omega,3)
would follow — contradicting CHANG's theorem.  So Chang FORCES finite cutoffs in
(s,C) for each fixed M.  This lab finds the first ones ("Chang numbers").

Verifier: complete structural DFS for an independent tower (top part first),
CROSS-VALIDATED against brute-force subset enumeration on small ambients.
CEGAR with Glucose3; triangle clauses upfront (lazy above a vertex threshold).
"""

import itertools, time, sys
from pysat.solvers import Glucose3

def ambient(s, C):
    return list(itertools.product(range(C), repeat=s))

def split(u, v):
    for i in range(len(u)):
        if u[i] != v[i]:
            return i
    return None

def is_binary_grid(S):
    """S: sorted tuple of vectors. Binary h-grid check per THM-460 B3(ii)."""
    n = len(S)
    if n == 1:
        return True
    if n % 2:
        return False
    A, B = S[:n // 2], S[n // 2:]
    p = split(A[0], B[0])
    for a in A:
        for b in B:
            if split(a, b) != p:
                return False
    for half in (A, B):
        for x, y in itertools.combinations(half, 2):
            sp = split(x, y)
            if sp is None or sp <= p:
                return False
    return is_binary_grid(A) and is_binary_grid(B)

# ----------------------------------------------------------------------------------
# complete independent-tower search (structural DFS, top part first)
# ----------------------------------------------------------------------------------

class TowerFinder:
    def __init__(self, s, C, M, L, idx):
        self.s, self.C, self.M = s, C, M
        self.L = L                      # sorted ambient
        self.idx = idx
        self.byprefix = {}              # cache: elements grouped for pruning (unused now)

    def find_independent_tower(self, adj):
        """Return list of parts [B_1..B_M] (each a list of vectors) independent as a
        whole, or None.  Complete search."""
        L, idx, M = self.L, self.idx, self.M
        N = len(L)

        # enumerate independent binary j-grids with min element index > lo,
        # all elements non-adjacent to 'taken'
        def grids(j, lo, taken):
            """Yield sorted tuples of 2^j vectors forming an independent binary j-grid,
            min index > lo, independent from taken (list of indices)."""
            if j == 0:
                for i in range(lo + 1, N):
                    if all(not adj[i][t] for t in taken):
                        yield (i,)
                return
            # a binary j-grid = two (j-1)-grids A < B with common cross split p,
            # internal splits > p.  Structural recursion: choose A, then B compatible.
            for A in grids(j - 1, lo, taken):
                amax = A[-1]
                takenA = taken + list(A)
                for B in grids(j - 1, amax, takenA):
                    S = tuple(sorted(A + B))
                    # quick reject: A,B halves must be exactly A then B in order
                    if S[:len(A)] != tuple(sorted(A)) or S[len(A):] != tuple(sorted(B)):
                        continue
                    vecs = tuple(L[i] for i in S)
                    if is_binary_grid(vecs):
                        yield S
        # build tower from TOP part (most constrained) downward? simpler: bottom-up
        # with completeness: recursive over parts j=1..M, each above the previous.
        def build(j, lo, taken):
            if j > M:
                return []
            for G in grids(j, lo, taken):
                rest = build(j + 1, G[-1], taken + list(G))
                if rest is not None:
                    return [G] + rest
            return None

        res = build(1, -1, [])
        if res is None:
            return None
        return [[self.L[i] for i in part] for part in res]

# brute-force reference for cross-validation
def brute_independent_tower(s, C, M, L, idx, adj):
    N = len(L)
    sizes = [2 ** j for j in range(1, M + 1)]
    total = sum(sizes)

    def indep(ids):
        return all(not adj[a][b] for a, b in itertools.combinations(ids, 2))

    for combo in itertools.combinations(range(N), total):
        if not indep(combo):
            continue
        # parts in order
        pos = 0
        ok = True
        for sz in sizes:
            part = combo[pos:pos + sz]
            if not is_binary_grid(tuple(L[i] for i in part)):
                ok = False
                break
            pos += sz
        if ok:
            return combo
    return None

# ----------------------------------------------------------------------------------

def solve_chang(M, s, C, tlimit=1800, verbose=True):
    L = ambient(s, C)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)
    evar = {}
    cnt = 0
    for i in range(N):
        for j in range(i + 1, N):
            cnt += 1
            evar[(i, j)] = cnt
    sol = Glucose3()
    lazy_tri = N > 100
    if not lazy_tri:
        for i, j, k in itertools.combinations(range(N), 3):
            sol.add_clause([-evar[(i, j)], -evar[(i, k)], -evar[(j, k)]])
    tf = TowerFinder(s, C, M, L, idx)
    t0 = time.time()
    rounds = added_t = added_tri = 0
    while True:
        rounds += 1
        if time.time() - t0 > tlimit:
            if verbose: print(f"   TIMEOUT M={M} s={s} C={C} ({added_t} tower cls, {time.time()-t0:.0f}s)", flush=True)
            return None
        if not sol.solve():
            if verbose: print(f"   UNSAT  M={M} s={s} C={C}  (towers={added_t}, tri={added_tri}, {time.time()-t0:.1f}s)", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        adj = [[False] * N for _ in range(N)]
        edges = []
        for (i, j), v in evar.items():
            if v in model:
                adj[i][j] = adj[j][i] = True
                edges.append((i, j))
        if lazy_tri:
            tri = None
            for (i, j) in edges:
                for k in range(N):
                    if adj[i][k] and adj[j][k]:
                        tri = (i, j, k); break
                if tri: break
            if tri:
                i, j, k = tri
                a, b, c = sorted((i, j, k))
                sol.add_clause([-evar[(a, b)], -evar[(a, c)], -evar[(b, c)]])
                added_tri += 1
                continue
        bad = tf.find_independent_tower(adj)
        if bad is None:
            if verbose: print(f"   SAT    M={M} s={s} C={C}  ({len(edges)} edges, towers={added_t}, {time.time()-t0:.1f}s)", flush=True)
            return True
        ids = sorted(idx[v] for part in bad for v in part)
        sol.add_clause([evar[(a, b)] for a, b in itertools.combinations(ids, 2)])
        added_t += 1

def crossvalidate():
    import random
    rng = random.Random(11)
    print("=== cross-validation: structural tower finder vs brute force ===", flush=True)
    for (s, C, M) in ((2, 3, 2), (3, 2, 2), (2, 2, 2), (3, 3, 2)):
        L = ambient(s, C); idx = {v: i for i, v in enumerate(L)}; N = len(L)
        if sum(2 ** j for j in range(1, M + 1)) > N:
            continue
        tf = TowerFinder(s, C, M, L, idx)
        bad = 0
        trials = 25 if N <= 16 else 12
        for tr in range(trials):
            adj = [[False] * N for _ in range(N)]
            for i in range(N):
                for j in range(i + 1, N):
                    if rng.random() < 0.25 * rng.random() * 2:
                        adj[i][j] = adj[j][i] = True
            r1 = tf.find_independent_tower(adj)
            r2 = brute_independent_tower(s, C, M, L, idx, adj)
            if (r1 is None) != (r2 is None):
                bad += 1
                print(f"   DISAGREE s={s} C={C} M={M} trial={tr}: struct={r1 is not None} brute={r2 is not None}", flush=True)
            if r1 is not None:
                ids = [idx[v] for part in r1 for v in part]
                assert all(not adj[a][b] for a, b in itertools.combinations(ids, 2)), "non-independent!"
        print(f"   s={s} C={C} M={M}: {trials} trials, {'OK' if bad == 0 else f'{bad} DISAGREEMENTS'}", flush=True)
    print(flush=True)

def main():
    crossvalidate()
    print("=== M=2 Chang-shadow sweep (Chang's theorem forces UNSAT at some (s,C)) ===", flush=True)
    for (s, C) in ((2, 3), (3, 2), (2, 4), (3, 3), (2, 5), (4, 2), (2, 6), (3, 4), (4, 3), (2, 8)):
        r = solve_chang(2, s, C, tlimit=1200)
        if r is False:
            print(f"   *** CHANG CUTOFF at M=2: first UNSAT (s,C)=({s},{C}) in sweep order ***", flush=True)
    print("=== M=3 spot checks ===", flush=True)
    for (s, C) in ((3, 3), (4, 3), (3, 4)):
        r = solve_chang(3, s, C, tlimit=1200)

if __name__ == "__main__":
    main()
