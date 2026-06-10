#!/usr/bin/env python3
# Chang M=2 tower game at (3,3) and (2,6) with BATCHED CEGAR (up to 64 violated towers
# per round — every tower clause is valid, so batching is sound).
# mac-mini-2026-06-10-S1 (HYP-2375, THM-460 D; POKE Steering Task 1.1)
import itertools, time
from pysat.solvers import Glucose3
from erdos592_chang_towers_macmini_s2 import ambient, is_binary_grid
from erdos592_chang_towers_v2_macmini_s2 import all_2grids

def solve_batched(s, C, tlimit=14000, K=64):
    L = ambient(s, C); N = len(L)
    grids2 = all_2grids(L)
    print(f"[{s},{C}] N={N}, 2-grids={len(grids2)}", flush=True)
    evar = {}; cnt = 0
    for i in range(N):
        for j in range(i + 1, N):
            cnt += 1; evar[(i, j)] = cnt
    sol = Glucose3()
    for i, j, k in itertools.combinations(range(N), 3):
        sol.add_clause([-evar[(i, j)], -evar[(i, k)], -evar[(j, k)]])
    def ep(i, j): return evar[(min(i, j), max(i, j))]
    t0 = time.time(); added = 0; rounds = 0
    while True:
        rounds += 1
        if time.time() - t0 > tlimit:
            print(f"   TIMEOUT (s,C)=({s},{C}) rounds={rounds} clauses={added} {time.time()-t0:.0f}s", flush=True)
            return None
        if not sol.solve():
            print(f"   UNSAT  (s,C)=({s},{C})  rounds={rounds} clauses={added} {time.time()-t0:.1f}s", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        nb = [0] * N
        for (i, j), v in evar.items():
            if v in model:
                nb[i] |= 1 << j; nb[j] |= 1 << i
        batch = 0
        for g in grids2:
            if batch >= K: break
            if any((nb[a] >> b) & 1 for a, b in itertools.combinations(g, 2)): continue
            mask = (1 << g[0]) - 1
            for x in g: mask &= ~nb[x]; mask &= ~(1 << x)
            ids = []
            m = mask
            while m:
                b = m & -m; ids.append(b.bit_length() - 1); m ^= b
            pr = None
            for a, b in itertools.combinations(ids, 2):
                if not (nb[a] >> b) & 1: pr = (a, b); break
            if pr:
                tower = list(pr) + list(g)
                sol.add_clause([ep(a, b) for a, b in itertools.combinations(tower, 2)])
                added += 1; batch += 1
        if batch == 0:
            edges = sum(1 for v in evar.values() if v in model)
            print(f"   SAT    (s,C)=({s},{C})  {edges} edges, rounds={rounds} clauses={added} {time.time()-t0:.1f}s", flush=True)
            return True

if __name__ == "__main__":
    for (s, C) in ((3, 3), (2, 6), (4, 2)):
        solve_batched(s, C)
