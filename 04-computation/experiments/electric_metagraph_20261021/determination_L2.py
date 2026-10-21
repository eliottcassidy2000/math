#!/usr/bin/env python3
"""Second-order determination: port types carry the PARTNER's invariants computed from the tiling itself -- level L1: the
partner's score multiset (THM-646); level L2: partner scores + partner curl profile (sorted vector of cyclic-triangle counts
per vertex); level L3: L2 + the partner's sorted out-neighbour score lists (second scores). The compatibility of a closed port
with a class requires the class invariants to equal the port's partner invariants, and the open port of the forced type
(n-a, n-2-b, sym) to exist. Counts CP-SAT solutions (cap 2e5) and whether the actual merged metagraph is among them."""
import sys, json, itertools
from collections import Counter, defaultdict
from ortools.sat.python import cp_model
from determination import tiles_of, tournament_from_mask, ham_paths, is_anti
def c3_profile(n, adjm):
    prof = [0] * n
    for a, b, c in itertools.combinations(range(n), 3):
        if (adjm[a][b] and adjm[b][c] and adjm[c][a]) or (adjm[a][c] and adjm[c][b] and adjm[b][a]):
            prof[a] += 1; prof[b] += 1; prof[c] += 1
    return tuple(sorted(prof))
def second_scores(n, adjm):
    sc = [sum(r) for r in adjm]
    return tuple(sorted(tuple(sorted(sc[v] for v in range(n) if adjm[u][v])) for u in range(n)))
def invariants(n, adjm, level):
    sc = tuple(sorted(sum(r) for r in adjm))
    if level == 1: return (sc,)
    if level == 2: return (sc, c3_profile(n, adjm))
    return (sc, c3_profile(n, adjm), second_scores(n, adjm))
def relabel_along(n, adjm, P):
    # tournament relabelled so that P becomes the base path n -> ... -> 1: new vertex label of P[i] is n-1-i (0-based: index n-1-i)
    pos = {P[i]: n - 1 - i for i in range(n)}
    B = [[0] * n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if adjm[u][v]: B[pos[u]][pos[v]] = 1
    return B
def complement_nonpath(n, B):
    C = [[0] * n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u == v: continue
            path_pair = abs(u - v) == 1
            if path_pair: C[u][v] = B[u][v]
            else: C[u][v] = B[v][u]
    return C
for level in [int(a) for a in sys.argv[2:]] or [1, 2, 3]:
    for n in [int(sys.argv[1])] if len(sys.argv) > 1 else [5, 6, 7]:
        data = json.load(open(f'electric_n{n}.json')); nodes = data['nodes']; N = len(nodes)
        blue = data['blue_closed_to_open']; black = data['black_closed_to_open']
        classes = []
        for u in nodes:
            classes.append((u['id'], u['rep_mask']))
            if not u['SC'] and u['conv_mask'] is not None: classes.append((u['id'], u['conv_mask']))
        cinv = {}; cports = {}; oports = {}
        for ci, (uid, mk) in enumerate(classes):
            adjm = tournament_from_mask(n, mk); A = nodes[uid]['Aut']; cinv[ci] = invariants(n, adjm, level)
            sc = [sum(r) for r in adjm]; cc = Counter(); oo = Counter()
            for P in ham_paths(n, adjm):
                sym = is_anti(n, adjm, P); closed = bool(adjm[P[-1]][P[0]])
                if closed:
                    B = relabel_along(n, adjm, P); Cm = complement_nonpath(n, B)
                    cc[(sc[P[0]], sc[P[-1]], sym, invariants(n, Cm, level))] += 1
                else:
                    oo[(sc[P[0]], sc[P[-1]], sym)] += 1
            cports[ci] = Counter({k: v // A for k, v in cc.items()}); oports[ci] = Counter({k: v // A for k, v in oo.items()})
        by_inv = defaultdict(list)
        for ci in range(len(classes)): by_inv[cinv[ci]].append(ci)
        model = cp_model.CpModel(); X = {}; forced = 0; free = 0
        for ci in range(len(classes)):
            for (a, b, sym, tinv), k in cports[ci].items():
                otype = (n - a, n - 2 - b, sym)
                cands = [cj for cj in by_inv.get(tinv, []) if oports[cj].get(otype, 0) > 0]
                assert cands, (n, level, ci, (a, b, sym))
                vs = []
                for cj in cands:
                    X[(ci, (a, b, sym, tinv), cj)] = model.NewIntVar(0, k, f"x_{len(X)}"); vs.append(X[(ci, (a, b, sym, tinv), cj)])
                model.Add(sum(vs) == k); forced += (len(cands) == 1); free += (len(cands) > 1)
        for cj in range(len(classes)):
            for otype, k in oports[cj].items():
                inc = [x for (ci, (a, b, sym, tinv), cj2), x in X.items() if cj2 == cj and (n - a, n - 2 - b, sym) == otype]
                model.Add(sum(inc) == k)
        conv_of = {}
        for ci, (uid, mk) in enumerate(classes):
            conv_of[ci] = ci if nodes[uid]['SC'] else [cj for cj, (u2, m2) in enumerate(classes) if u2 == uid and cj != ci][0]
        agg = defaultdict(list)
        for (ci, (a, b, sym, tinv), cj), x in X.items():
            if not sym: agg[(ci, cj)].append(x)
        for (ci, cj), xs in agg.items():
            ck, cl = conv_of[ci], conv_of[cj]
            if (ck, cl) in agg and (ck, cl) != (ci, cj): model.Add(sum(xs) == sum(agg[(ck, cl)]))
            if (ck, cl) == (ci, cj):
                h = model.NewIntVar(0, 10 ** 6, f"h_{ci}_{cj}"); model.Add(sum(xs) == 2 * h)
        solver = cp_model.CpSolver(); solver.parameters.enumerate_all_solutions = True; solver.parameters.num_search_workers = 1; solver.parameters.max_time_in_seconds = 600
        node_of_class = [uid for uid, mk in classes]
        class Cb(cp_model.CpSolverSolutionCallback):
            def __init__(s): super().__init__(); s.count = 0; s.matched = 0
            def on_solution_callback(s):
                s.count += 1; ab = Counter(); ak = Counter()
                for (ci, (a, b, sym, tinv), cj), x in X.items(): (ab if sym else ak)[(node_of_class[ci], node_of_class[cj])] += s.Value(x)
                if all(ab[(i, j)] == blue[i][j] for i in range(N) for j in range(N)) and all(ak[(i, j)] == black[i][j] for i in range(N) for j in range(N)): s.matched += 1
                if s.count >= 200000: s.StopSearch()
        cb = Cb(); status = solver.Solve(model, cb)
        print(f"level L{level}, n = {n}: classes {len(classes)}, port classes forced {forced} / free {free}, variables {len(X)}; solutions {cb.count} (status {solver.StatusName(status)}); actual merged metagraph among them: {cb.matched}", flush=True)
