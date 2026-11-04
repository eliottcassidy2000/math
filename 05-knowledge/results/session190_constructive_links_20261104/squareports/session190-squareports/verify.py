#!/usr/bin/env python3
"""Exact endpoint-pair transfer, direct DFS replay, and deletion response tests.
No dependencies. python3 verify.py > RESULTS.json
"""
from collections import Counter
from itertools import combinations
from math import isqrt
import json
from direct import direct


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def square_graph(n):
    return {v: {w for w in range(1, n + 1)
                if v != w and isqrt(v + w) ** 2 == v + w}
            for v in range(1, n + 1)}


def edges_graph(vertices, edges):
    G = {v: set() for v in vertices}
    for a, b in edges:
        G[a].add(b)
        G[b].add(a)
    return G


def delete(G, v):
    return {w: nb - {v} for w, nb in G.items() if w != v}


def transitions(P, z, S):
    """Canonical pairs encode components; (v,v) is a singleton, not an edge."""
    yield tuple(sorted(P + ((z, z),))), ()
    available = []
    for i, (a, b) in enumerate(P):
        if a in S:
            available.append((a, b, i))
        if b != a and b in S:
            available.append((b, a, i))
    for a, b, i in available:
        R = tuple(sorted(P[:i] + P[i+1:] + (tuple(sorted((b, z))),)))
        yield R, ((a, z),)
    for (a, b, i), (c, d, j) in combinations(available, 2):
        if i == j:
            continue
        R = tuple(sorted(tuple(p for k, p in enumerate(P) if k != i and k != j)
                         + (tuple(sorted((b, d))),)))
        yield R, ((a, z), (c, z))


def full_transfer(C, z, S):
    T = Counter()
    for P, count in C.items():
        for R, _ in transitions(P, z, S):
            T[R] += count
    return T


def full_profile(G):
    C = Counter({(): 1})
    for z in sorted(G):
        C = full_transfer(C, z, G[z])
    return C


def solve(G, retain_penultimate=False):
    """Exact target-aware transfer. Only completion-impossible states are pruned."""
    V = sorted(G)
    N = len(V)
    if not N:
        return Counter(), 0, Counter()
    pos = {v: i for i, v in enumerate(V)}
    last = {v: max([pos[v]] + [pos[w] for w in G[v]]) for v in V}
    C = Counter({(): 1})
    peak = 1
    penultimate = Counter()
    for i, z in enumerate(V):
        if retain_penultimate and i == N - 1:
            penultimate = C.copy()
        T = Counter()
        remaining = N - i - 1
        for P, count in C.items():
            for R, _ in transitions(P, z, G[z]):
                if len(R) > remaining + 1:
                    continue
                dead = 0
                bad = False
                for a, b in R:
                    da, db = last[a] <= i, last[b] <= i
                    if a == b:
                        if da and N > 1:
                            bad = True
                            break
                    else:
                        dead += da + db
                        if da and db and remaining:
                            bad = True
                            break
                if bad or dead > 2:
                    continue
                T[R] += count
        C = T
        peak = max(peak, len(C))
    return Counter({P[0]: c for P, c in C.items()}), peak, penultimate


def subset_profile(G):
    """Independent edge-subset enumeration, including all component numbers."""
    V = sorted(G)
    E = [(a, b) for a in V for b in sorted(G[a]) if a < b]
    out = Counter()
    for mask in range(1 << len(E)):
        parent = {v: v for v in V}
        deg = Counter()
        def root(v):
            while parent[v] != v:
                v = parent[v]
            return v
        valid = True
        for j, (a, b) in enumerate(E):
            if not (mask >> j) & 1:
                continue
            x, y = root(a), root(b)
            deg[a] += 1
            deg[b] += 1
            if x == y or deg[a] > 2 or deg[b] > 2:
                valid = False
                break
            parent[x] = y
        if not valid:
            continue
        comps = {}
        for v in V:
            comps.setdefault(root(v), []).append(v)
        P = []
        for comp in comps.values():
            if len(comp) == 1:
                P.append((comp[0], comp[0]))
            else:
                ends = tuple(v for v in comp if deg[v] == 1)
                need(len(ends) == 2, 'subset component is not a path')
                P.append(tuple(sorted(ends)))
        out[tuple(sorted(P))] += 1
    return out


def graph_family(n):
    E = list(combinations(range(n), 2))
    for mask in range(1 << len(E)):
        yield edges_graph(range(n), (e for j, e in enumerate(E) if (mask >> j) & 1))


def encoded(C):
    return [{'endpoints': list(p), 'count': c} for p, c in sorted(C.items())]


def port_encoded(C):
    return [{'component_pairs': [list(p) for p in P], 'count': c}
            for P, c in sorted(C.items())]


def check_witness(G, pair, path):
    need(set(path) == set(G) and len(path) == len(G), 'witness vertex set')
    need(tuple(sorted((path[0], path[-1]))) == pair, 'witness endpoints')
    need(all(b in G[a] for a, b in zip(path, path[1:])), 'witness edge')


def main():
    R = {'status': 'PROVED transfer and FINITE-EXACT experiments; not an all-order square-sum theorem'}
    counts = {'full_profile_graphs': 0, 'endpoint_profile_graphs': 0,
              'deletion_update_cases': 0, 'four_port_response_cases': 0}
    for n in range(1, 6):
        for G in graph_family(n):
            F = full_profile(G)
            need(F == subset_profile(G), ('full profile', n, G))
            counts['full_profile_graphs'] += 1
            C, _, _ = solve(G)
            D, W = direct(G)
            need(C == D, ('endpoint profile', n, G))
            counts['endpoint_profile_graphs'] += 1
            if n > 1:
                z = n - 1
                old = delete(G, z)
                for v in old:
                    left = full_profile(delete(G, v))
                    right = full_transfer(full_profile(delete(old, v)), z, G[z] - {v})
                    need(left == right, ('deletion update', n, v))
                    counts['deletion_update_cases'] += 1
            if n == 4:
                pairs = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
                X = [F[p] for p in pairs]
                Y = []
                for selected, output in pairs:
                    H = edges_graph(range(5), [(a, b) for a in G for b in G[a] if a < b]
                                    + [(a, 4) for a in selected])
                    y = direct(H)[0][output]
                    Y.append(y)
                need(Y == [X[1]+X[2], X[0]+X[2], X[0]+X[1]], 'four-port response')
                need(sum(Y) % 2 == 0, 'response lattice parity')
                counts['four_port_response_cases'] += 3
    R['exhaustive_tests'] = counts
    H1 = edges_graph(range(5), [(0,1),(1,2),(2,3),(3,4)])
    H2 = edges_graph(range(5), [(0,1),(1,3),(3,2),(2,4)])
    old1, old2 = direct(H1)[0], direct(H2)[0]
    need(old1 == old2 == Counter({(0,4):1}), 'hostile old profiles')
    new = []
    for G in (H1,H2):
        K = {v:set(nb) for v,nb in G.items()}
        K[5] = {1,2}
        K[1].add(5)
        K[2].add(5)
        new.append(direct(K)[0])
    need(new == [Counter({(0,4):1}),Counter()], 'hostile new profiles')
    R['endpoint_matrix_hostile'] = {'old_path_orders': [[0,1,2,3,4],[0,1,3,2,4]],
        'old_profile': encoded(old1), 'new_vertex': 5, 'new_neighbors': [1,2],
        'new_profiles': [encoded(C) for C in new], 'positive_path': [0,1,5,2,3,4]}
    R['four_port_probe'] = {'matrix': [[0,1,1],[1,0,1],[1,1,0]],
        'determinant': 2, 'rational_rank': 3, 'mod2_rank': 2,
        'K4_two_cover_vector': [1,1,1], 'K4_response': [2,2,2]}
    rows = []
    for n in range(1,31):
        G = square_graph(n)
        C, peak, P = solve(G, retain_penultimate=True)
        D, W = direct(G)
        need(C == D, ('square endpoint matrix', n))
        for pair,path in W.items():
            check_witness(G,pair,path)
        extension = sum(c for (a,b),c in C.items() if n in (a,b)) if n > 1 else 0
        row = {'n':n, 'path_count':sum(C.values()), 'endpoint_pairs':encoded(C),
               'endpoint_pair_count':len(C), 'maximum_transfer_states':peak,
               'arrival_neighbors':sorted(G[n]), 'endpoint_extensions':extension,
               'two_component_joins':sum(C.values())-extension if n>1 else 0,
               'witnesses':[{'endpoints':list(pair),'path':path} for pair,path in sorted(W.items())]}
        if 23 <= n <= 30:
            row['completion_relevant_previous_profile'] = port_encoded(P)
            ext = join = 0
            for state,c in P.items():
                for target, extra in transitions(state,n,G[n]):
                    if len(target) == 1:
                        if len(extra) == 1: ext += c
                        if len(extra) == 2: join += c
            need(ext == extension and join == row['two_component_joins'], 'arrival decomposition')
        rows.append(row)
    R['square_endpoint_census'] = rows
    deletions=[]
    for n in (24,25,26,29,30):
        G=square_graph(n)
        cards=[]
        for v in G:
            H=delete(G,v)
            C,peak,_=solve(H)
            D,W=direct(H)
            need(C==D, ('square deletion matrix',n,v))
            for pair,path in W.items():check_witness(H,pair,path)
            cards.append({'deleted':v,'path_count':sum(C.values()),'endpoint_pair_count':len(C),
                          'endpoint_pairs':encoded(C),'maximum_transfer_states':peak,
                          'one_path':W[min(W)] if W else None})
        deletions.append({'n':n,'cards':cards})
    R['square_deletion_census']=deletions
    card_lookup={(r['n'],c['deleted']):c['path_count'] for r in deletions for c in r['cards']}
    need(card_lookup[25,3]==1 and card_lookup[26,3]==0,'deletion nonmonotonicity')
    need(rows[27]['path_count']==52 and rows[28]['path_count']==19,'scalar nonmonotonicity')
    need(rows[28]['endpoint_extensions']==0 and rows[28]['two_component_joins']==19,'Q29 compulsory internal arrival')
    R['finite_consequences']={'Q25_minus3_paths':1,'Q26_minus3_paths':0,
        'Q28_paths':52,'Q29_paths':19,'Q29_endpoint_extensions':0,
        'Q29_internal_arrival_paths':19}
    print(json.dumps(R,sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
