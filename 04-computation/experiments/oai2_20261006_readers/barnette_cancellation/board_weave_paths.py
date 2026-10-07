#!/usr/bin/env python3
"""8x8 board: the weave law (THM-4550 W) says every closed tour has exactly k moves inside the
outer ring 3 and k moves inside rings 1-2, with 1 <= k <= 8 (only 8 moves lie inside ring 3).
So any prescribed knight path with >= 9 moves inside rings 1-2 lies in no closed tour.
Search random such paths that PASS the local forcing closure (so the obstruction is global),
and find the fewest moves of a locally-consistent non-extendable path found this way.
Also: random locally-consistent paths with m = 6..12 moves, testing extendability, to locate
the first non-local failure."""
import random, sys, time
import board_paths as BP

N = 8
eid, E, adj = BP.board(N)
nV = N * N
ring = lambda v: max(abs(2 * (v // N) - (N - 1)), abs(2 * (v % N) - (N - 1))) // 2
inner = {v for v in range(nV) if ring(v) in (1, 2)}


def rand_path(m, rng, within=None):
    for _ in range(1000):
        v = rng.choice(sorted(within) if within else range(nV))
        path = [v]
        ok = True
        while len(path) < m + 1:
            nb = [w for w in adj[path[-1]] if w not in path and (within is None or w in within)]
            if not nb:
                ok = False; break
            path.append(rng.choice(nb))
        if ok:
            return path
    return None


def edges_of(path):
    return [eid[(min(a, b), max(a, b))] for a, b in zip(path, path[1:])]


rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
T0 = time.time()
# (1) paths inside rings 1-2 with 9 moves (all inner): locally consistent? extendable?
cnt = {"local_fail": 0, "local_ok": 0, "ext": 0}
for trial in range(60):
    P = rand_path(9, rng, within=inner)
    if P is None:
        continue
    EP = edges_of(P)
    if not BP.closure_ok(nV, E, adj, eid, EP):
        cnt["local_fail"] += 1; continue
    cnt["local_ok"] += 1
    r = BP.find_hc(E, nV, EP)
    if r:
        cnt["ext"] += 1
print(f"(1) 9-move paths inside rings 1-2: {cnt}  [weave law predicts ext = 0] ({time.time()-T0:.0f}s)", flush=True)
# (2) random paths (anywhere) with m moves passing the local closure: fraction non-extendable
for m in range(5, 13):
    tot = nonext = 0
    ex = None
    for trial in range(150):
        P = rand_path(m, rng)
        if P is None:
            continue
        EP = edges_of(P)
        if not BP.closure_ok(nV, E, adj, eid, EP):
            continue
        tot += 1
        r = BP.find_hc(E, nV, EP)
        if r is False:
            nonext += 1
            inn = sum(1 for a, b in zip(P, P[1:]) if a in inner and b in inner)
            ex = ex or ([divmod(x, N) for x in P], inn)
    print(f"(2) m={m}: {tot} locally-consistent random paths, {nonext} not in any closed tour; example {ex} "
          f"({time.time()-T0:.0f}s)", flush=True)
