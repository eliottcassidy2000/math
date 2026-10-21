#!/usr/bin/env python3
"""Determination problem (HYP-9345 (g)) in electrical form. Node-level data only: for each merged node u, its one or two
classes with representative tournaments; from each class the CLOSED tilings (Hamiltonian paths P with last -> first) and OPEN
tilings, each with the port type (score of start, score of end, grid-symmetric?). The complement line sends a closed port
(a, b, sym) of a node with score multiset S to an open port (n-a, n-2-b, sym) of a node with score multiset S' = {n-1-s} with
start n-a and end n-2-b (THM-646). Unknown: the transportation table m[(u,port) -> v] over nodes v with the right score multiset
and open ports of the forced type. Count the integer solutions with CP-SAT (cap 10^6); compare with the actual multigraph.
If exactly one solution exists, the node-level port data determine the merged metagraph (lines by colour, closed->open)."""
import sys, json, itertools
from collections import Counter, defaultdict
from ortools.sat.python import cp_model
def tiles_of(n): return [(a, b) for b in range(1, n - 1) for a in range(b + 2, n + 1)]
def tournament_from_mask(n, mask):
    tiles = tiles_of(n); adjm = [[0] * n for _ in range(n)]
    for i in range(1, n): adjm[i][i - 1] = 1
    for k, (a, b) in enumerate(tiles):
        if (mask >> k) & 1: adjm[b - 1][a - 1] = 1
        else: adjm[a - 1][b - 1] = 1
    return adjm
def ham_paths(n, adjm):
    adj = [[v for v in range(n) if adjm[u][v]] for u in range(n)]; out = []; full = (1 << n) - 1
    def rec(v, mask, path):
        if mask == full: out.append(tuple(path)); return
        for w in adj[v]:
            if not (mask >> w) & 1: path.append(w); rec(w, mask | (1 << w), path); path.pop()
    for v in range(n): rec(v, 1 << v, [v])
    return out
def is_anti(n, adjm, P):
    # positional reversal of P is an anti-automorphism iff for all i<j: arc(P_i,P_j) == arc(P_{n-1-j}, P_{n-1-i})
    for i in range(n):
        for j in range(i + 1, n):
            if adjm[P[i]][P[j]] != adjm[P[n - 1 - j]][P[n - 1 - i]]: return False
    return True
def main():
    for n in [int(a) for a in sys.argv[1:]] or [4, 5, 6, 7]:
        data = json.load(open(f'electric_n{n}.json')); nodes = data['nodes']; N = len(nodes)
        blue = data['blue_closed_to_open']; black = data['black_closed_to_open']
        # CLASS level: every node contributes its representative class and, if merged, the converse class
        classes = []   # (node id, mask)
        for u in nodes:
            classes.append((u['id'], u['rep_mask']))
            if not u['SC'] and u['conv_mask'] is not None: classes.append((u['id'], u['conv_mask']))
        cports = {}; oports = {}; cscores = {}
        for ci, (uid, mk) in enumerate(classes):
            adjm = tournament_from_mask(n, mk); sc = [sum(r) for r in adjm]; A = nodes[uid]['Aut']
            cc = Counter(); oo = Counter()
            for P in ham_paths(n, adjm):
                typ = (sc[P[0]], sc[P[-1]], is_anti(n, adjm, P))
                (cc if adjm[P[-1]][P[0]] else oo)[typ] += 1
            cports[ci] = Counter({k: v // A for k, v in cc.items()}); oports[ci] = Counter({k: v // A for k, v in oo.items()})
            for v in list(cc.values()) + list(oo.values()): assert v % A == 0
            cscores[ci] = tuple(sorted(sc))
        by_scores = defaultdict(list)
        for ci in range(len(classes)): by_scores[cscores[ci]].append(ci)
        model = cp_model.CpModel(); X = {}; free = 0; forced = 0
        for ci in range(len(classes)):
            for (a, b, sym), k in cports[ci].items():
                S2 = Counter(n - 1 - s for s in cscores[ci]); S2[n - 1 - a] -= 1; S2[n - a] += 1; S2[n - 1 - b] -= 1; S2[n - 2 - b] += 1
                target = tuple(sorted(S2.elements())); otype = (n - a, n - 2 - b, sym)
                cands = [cj for cj in by_scores.get(target, []) if oports[cj].get(otype, 0) > 0]
                assert cands, (n, ci, (a, b, sym), target)
                vs = []
                for cj in cands:
                    X[(ci, (a, b, sym), cj)] = model.NewIntVar(0, k, f"x_{ci}_{a}_{b}_{int(sym)}_{cj}"); vs.append(X[(ci, (a, b, sym), cj)])
                model.Add(sum(vs) == k)
                if len(cands) == 1: forced += 1
                else: free += 1
        for cj in range(len(classes)):
            for otype, k in oports[cj].items():
                inc = [x for (ci, (a, b, sym), cj2), x in X.items() if cj2 == cj and (n - a, n - 2 - b, sym) == otype]
                model.Add(sum(inc) == k)
        # black class-to-class multiplicities: the mirror sends a closed black tiling of class C to a closed black tiling of C^rev with partner D^rev,
        # so m_black[C -> D] = m_black[C^rev -> D^rev] (no parity constraint at class level in general); impose this symmetry
        conv_of = {}
        for ci, (uid, mk) in enumerate(classes):
            if nodes[uid]['SC']: conv_of[ci] = ci
            else:
                others = [cj for cj, (u2, m2) in enumerate(classes) if u2 == uid and cj != ci]; conv_of[ci] = others[0]
        agg = defaultdict(list)
        for (ci, (a, b, sym), cj), x in X.items():
            if not sym: agg[(ci, cj)].append(x)
        for (ci, cj), xs in agg.items():
            ck, cl = conv_of[ci], conv_of[cj]
            if (ck, cl) in agg and (ck, cl) != (ci, cj): model.Add(sum(xs) == sum(agg[(ck, cl)]))
            if (ck, cl) == (ci, cj):
                # a closed black tiling and its mirror image are both in C with the same partner class: even multiplicity
                h = model.NewIntVar(0, 10 ** 6, f"h_{ci}_{cj}"); model.Add(sum(xs) == 2 * h)
        solver = cp_model.CpSolver(); solver.parameters.enumerate_all_solutions = True; solver.parameters.num_search_workers = 1
        node_of_class = [uid for uid, mk in classes]
        class Cb(cp_model.CpSolverSolutionCallback):
            def __init__(s): super().__init__(); s.count = 0; s.matched = 0
            def on_solution_callback(s):
                s.count += 1
                agg_b = Counter(); agg_k = Counter()
                for (ci, (a, b, sym), cj), x in X.items():
                    (agg_b if sym else agg_k)[(node_of_class[ci], node_of_class[cj])] += s.Value(x)
                if all(agg_b[(i, j)] == blue[i][j] for i in range(N) for j in range(N)) and all(agg_k[(i, j)] == black[i][j] for i in range(N) for j in range(N)): s.matched += 1
                if s.count >= 200000: s.StopSearch()
        cb = Cb(); solver.parameters.max_time_in_seconds = 900; status = solver.Solve(model, cb)
        print(f"n = {n}: classes {len(classes)}, nodes {N}, port variables {len(X)} (forced port classes {forced}, free {free}); solutions found {cb.count} (cap 2e5, status {solver.StatusName(status)}); merged actual multigraph among them: {cb.matched}", flush=True)

if __name__ == '__main__':
    main()
