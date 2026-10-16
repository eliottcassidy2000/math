#!/usr/bin/env python3
"""S38 audit (independent referee): plain backtracking (no SAT) cross-check of the finding of
audit_nineteen_vertex_exhaustive.py that NO graph on 19 vertices is triangle-free, has all codegrees <= 2 and
minimum degree >= 5 -- i.e. no 19-vertex five-core of the THM-4661 type exists even combinatorially.

Structure (forced by double counting, see the report): a degree-6 root z, its neighbours u_1..u_6 of degree 5,
and twelve second neighbours w_P indexed by the twelve edges P of the octahedron K_{2,2,2} on {1..6} with parts
{1,2}, {3,4}, {5,6}; w_P is adjacent to u_i for i in P.  Remaining unknown: the graph on the twelve w's.

The script enumerates every subset of the 66 possible w-w pairs by backtracking over pairs in a fixed order with
the constraints: W-degree in {3,4}; triangle-free (w_P ~ w_Q forces P, Q disjoint, and no triangle inside W);
codegree <= 2 for every pair of the 19 vertices.  It also reports the solution counts of relaxations that drop one
family of constraints, to locate the obstruction, and prints a human-readable obstruction count.
"""
import itertools, json, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_nineteen_vertex_backtrack.json')

match = [(1, 2), (3, 4), (5, 6)]
pairs = [p for p in itertools.combinations(range(1, 7), 2) if p not in match]
assert len(pairs) == 12
W = list(range(12))


def run(check_ww_codegree=True, check_uw_codegree=True, deg_range=(3, 4), allow_shared=False):
    cand = [(a, b) for a, b in itertools.combinations(W, 2)
            if allow_shared or not (set(pairs[a]) & set(pairs[b]))]
    adj = [set() for _ in W]
    deg = [0] * 12
    sols = []
    lo, hi = deg_range
    # remaining candidate count per vertex (for pruning the lower degree bound)
    remaining = [sum(1 for e in cand if v in e) for v in W]

    def ok_add(a, b):
        if adj[a] & adj[b]:
            return False                      # triangle inside W
        if set(pairs[a]) & set(pairs[b]):
            return False                      # triangle with a shared u (only reachable if allow_shared)
        if deg[a] >= hi or deg[b] >= hi:
            return False
        # codegree after adding a-b: pairs (x, b) for x in adj[a] gain common neighbour a, and (x, a) gain b
        if check_ww_codegree:
            for x in adj[a]:
                if x == b:
                    continue
                c = len(adj[x] & adj[b]) + 1 + len(set(pairs[x]) & set(pairs[b]))
                if c > 2:
                    return False
            for x in adj[b]:
                if x == a:
                    continue
                c = len(adj[x] & adj[a]) + 1 + len(set(pairs[x]) & set(pairs[a]))
                if c > 2:
                    return False
        if check_uw_codegree:
            # u_i (i not in P_b) and w_b: common neighbours = W-neighbours of w_b containing i
            for i in pairs[a]:
                if i not in pairs[b] and sum(1 for x in adj[b] if i in pairs[x]) + 1 > 2:
                    return False
            for i in pairs[b]:
                if i not in pairs[a] and sum(1 for x in adj[a] if i in pairs[x]) + 1 > 2:
                    return False
        return True

    def rec(k):
        # prune: every vertex must still be able to reach degree lo
        for v in W:
            if deg[v] + remaining[v] < lo:
                return
        if k == len(cand):
            if all(lo <= d <= hi for d in deg):
                sols.append(frozenset(e for e in cand if e[1] in adj[e[0]]))
            return
        a, b = cand[k]
        remaining[a] -= 1; remaining[b] -= 1
        if ok_add(a, b):
            adj[a].add(b); adj[b].add(a); deg[a] += 1; deg[b] += 1
            rec(k + 1)
            adj[a].discard(b); adj[b].discard(a); deg[a] -= 1; deg[b] -= 1
        rec(k + 1)
        remaining[a] += 1; remaining[b] += 1

    t0 = time.time()
    rec(0)
    return len(sols), round(time.time() - t0, 2), sols


out = {}
for label, kw in [('full constraints', {}),
                  ('drop w-w codegree', dict(check_ww_codegree=False)),
                  ('drop u-w codegree', dict(check_uw_codegree=False)),
                  ('drop both codegree families', dict(check_ww_codegree=False, check_uw_codegree=False))]:
    cnt, sec, sols = run(**kw)
    out[label] = dict(solutions=cnt, seconds=sec)
    print(f'{label:32s}: {cnt} labelled W-graphs ({sec}s)', flush=True)

# Human-readable obstruction: for each w, the five admissible partners {a'b', a'c, a'c', b'c, b'c'}; u-w codegree
# allows at most two of the three containing a' and at most two of the three containing b'.
adm = {}
for a, P in enumerate(pairs):
    part = {1: 2, 2: 1, 3: 4, 4: 3, 5: 6, 6: 5}
    x, y = P
    xp, yp = part[x], part[y]
    rest = [i for i in range(1, 7) if i not in (x, y, xp, yp)]
    adm[P] = sorted(tuple(sorted(Q)) for Q in pairs if not (set(Q) & set(P)))
    assert len(adm[P]) == 5
out['admissible_partners'] = {str(k): v for k, v in adm.items()}
json.dump(out, open(OUTJ, 'w'), indent=1)
print('wrote', OUTJ)
