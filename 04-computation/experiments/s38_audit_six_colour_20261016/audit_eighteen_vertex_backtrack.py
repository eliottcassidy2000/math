#!/usr/bin/env python3
"""S38 audit (independent referee): plain backtracking (no SAT) cross-check of audit_eighteen_vertex_five_cores.py:
no 5-regular graph on 18 vertices is triangle-free with all codegrees <= 2.

Same case split (root z, N(z) = u_1..u_5, twelve further vertices with 2/1/0 neighbours in N(z): cases A, B, C1, C2);
the backtracking adds edges among the twelve further vertices in a fixed pair order, maintaining exact degree 5
(upper bound while adding, lower bound by remaining capacity), triangle-freeness and codegree <= 2 for every pair
of the 18 vertices (the fixed part included).  As a positive control the same routine, run on the 16-vertex tight
structure (N(z) = 5, ten pair vertices), must find the Clebsch graph.
"""
import itertools, json, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_eighteen_vertex_backtrack.json')


def case_att(name):
    allp = list(itertools.combinations(range(1, 6), 2))
    table = {'A': ([], [], 2), 'B': ([(1, 2)], [1, 2], 1), 'C1': ([(1, 2), (3, 4)], [1, 2, 3, 4], 0),
             'C2': ([(1, 2), (1, 3)], [1, 1, 2, 3], 0), 'clebsch16': ([], [], 0)}
    miss, singles, far = table[name]
    return [set(p) for p in allp if p not in miss] + [{s} for s in singles] + [set() for _ in range(far)]


def search(att, deg_target=5):
    m = len(att)
    n = 6 + m
    adj = [set() for _ in range(n)]

    def add(a, b):
        adj[a].add(b); adj[b].add(a)

    for u in range(1, 6):
        add(0, u)
    for i, A in enumerate(att):
        for u in A:
            add(u, 6 + i)
    others = list(range(6, n))
    # sanity of the fixed part
    for a, b in itertools.combinations(range(n), 2):
        assert not (b in adj[a] and adj[a] & adj[b])
        assert len(adj[a] & adj[b]) <= 2
    free = list(itertools.combinations(others, 2))
    remaining = {v: sum(1 for e in free if v in e) for v in others}
    sols = []

    def can_add(a, b):
        if len(adj[a]) >= deg_target or len(adj[b]) >= deg_target:
            return False
        if adj[a] & adj[b]:
            return False                                # triangle
        for x in adj[a]:                                # x and b gain common neighbour a
            if len(adj[x] & adj[b]) + 1 > 2:
                return False
        for x in adj[b]:
            if len(adj[x] & adj[a]) + 1 > 2:
                return False
        return True

    def rec(k):
        for v in others:
            if len(adj[v]) + remaining[v] < deg_target:
                return
        if k == len(free):
            if all(len(adj[v]) == deg_target for v in others):
                sols.append(sorted((a, b) for a, b in free if b in adj[a]))
            return
        a, b = free[k]
        remaining[a] -= 1; remaining[b] -= 1
        if can_add(a, b):
            add(a, b)
            rec(k + 1)
            adj[a].discard(b); adj[b].discard(a)
        rec(k + 1)
        remaining[a] += 1; remaining[b] += 1

    rec(0)
    return sols


out = {}
for name in ('clebsch16', 'A', 'B', 'C1', 'C2'):
    t0 = time.time()
    sols = search(case_att(name))
    out[name] = dict(labelled_solutions=len(sols), seconds=round(time.time() - t0, 2))
    print(f'{name:10s}: {len(sols)} labelled completions ({out[name]["seconds"]}s)', flush=True)
assert out['clebsch16']['labelled_solutions'] >= 1
json.dump(out, open(OUTJ, 'w'), indent=1)
print('wrote', OUTJ)
