#!/usr/bin/env python3
"""Exact bounded certificates for the square-sum cycle threshold 31 -> 32.

No third-party dependencies. G31 uses only C(15,8)=6435 subsets;
G32 and the specified-endpoint G31 path use degree-obligation branching.
Checks use explicit exceptions and remain active under python -O.
"""
import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations, product
from math import isqrt
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def edges(n):
    return {(u, v) for u in range(1, n + 1) for v in range(u + 1, n + 1)
            if isqrt(u + v) ** 2 == u + v}


def incident(es, v):
    return {e for e in es if v in e}


def propagate(n, original, target, forced=None, allowed=None, trace=None):
    forced = set() if forced is None else set(forced)
    allowed = set(original) if allowed is None else set(allowed)
    while True:
        old = len(forced), len(allowed)
        for v in range(1, n + 1):
            a, f = incident(allowed, v), incident(forced, v)
            if len(a) < target[v] or len(f) > target[v]:
                return None
            if len(a) == target[v] and a - forced:
                if trace is not None:
                    trace.append({'action': 'force', 'vertex': v,
                                  'edges': sorted(a - forced),
                                  'reason': 'allowed degree equals target'})
                forced |= a
                f = incident(forced, v)
            if len(f) == target[v] and a - f:
                if trace is not None:
                    trace.append({'action': 'delete', 'vertex': v,
                                  'edges': sorted(a - f),
                                  'reason': 'forced degree equals target'})
                allowed -= a - f
        if old == (len(forced), len(allowed)):
            return forced, allowed


def components(n, selected):
    adj = {v: [] for v in range(1, n + 1)}
    for u, v in sorted(selected):
        adj[u].append(v)
        adj[v].append(u)
    unseen, answer = set(adj), []
    while unseen:
        seed = min(unseen)
        todo, seen = [seed], set()
        while todo:
            x = todo.pop()
            if x in seen:
                continue
            seen.add(x)
            todo.extend(adj[x])
        unseen -= seen
        ends = sorted(v for v in seen if len(adj[v]) < 2)
        start = ends[0] if ends else min(seen)
        path, previous, current = [start], None, start
        while True:
            candidates = sorted(v for v in adj[current] if v != previous)
            if not candidates:
                break
            nxt = candidates[0]
            if nxt == start:
                break
            require(nxt not in path, 'unexpected repeated vertex')
            path.append(nxt)
            previous, current = current, nxt
        require(set(path) == seen, 'component is not a path or cycle')
        answer.append(path)
    return sorted(answer, key=lambda p: (len(p), p))


def degree_ok(n, selected, target):
    degree = Counter(v for e in selected for v in e)
    return all(degree[v] == target[v] for v in range(1, n + 1))


def subset_completions(n, forced, allowed, target):
    optional = sorted(allowed - forced)
    need = sum(target.values()) // 2 - len(forced)
    solutions, trials = [], 0
    for subset in combinations(optional, need):
        trials += 1
        selected = forced | set(subset)
        if degree_ok(n, selected, target):
            solutions.append(selected)
    return solutions, trials


def obligation_completions(n, original, target):
    """Choose every way to satisfy a least-flexible vertex, then propagate.

    A solution chooses exactly target[v] - forced_degree[v] of its remaining
    incident edges; these choices partition all completions without overlap.
    No connectivity or proper-cycle pruning is used, so all 2-factors survive.
    """
    solutions, nodes = [], 0

    def visit(forced, allowed):
        nonlocal nodes
        nodes += 1
        result = propagate(n, original, target, forced, allowed)
        if result is None:
            return
        forced, allowed = result
        if forced == allowed:
            require(degree_ok(n, forced, target), 'bad terminal degrees')
            solutions.append(forced)
            return
        candidates = []
        for v in range(1, n + 1):
            undecided = sorted(incident(allowed - forced, v))
            need = target[v] - len(incident(forced, v))
            if undecided:
                ways = list(combinations(undecided, need))
                candidates.append((len(ways), v, undecided, ways))
        _, _, undecided, ways = min(candidates, key=lambda x: (x[0], x[1]))
        for chosen in ways:
            chosen = set(chosen)
            visit(forced | chosen, allowed - (set(undecided) - chosen))

    visit(set(), set(original))
    require(len({tuple(sorted(s)) for s in solutions}) == len(solutions),
            'duplicate obligation completion')
    return solutions, nodes


def cycle_record(n, selected, forced):
    cs = components(n, selected)
    return {'component_lengths': sorted(map(len, cs)), 'cycles': cs,
            'optional_selected': sorted(selected - forced)}


def edge_set_of_path(path, close=False):
    pairs = list(zip(path, path[1:]))
    if close:
        pairs.append((path[-1], path[0]))
    return {tuple(sorted(e)) for e in pairs}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='certificate.json')
    args = parser.parse_args()
    e31, e32 = edges(31), edges(32)
    t31 = {v: 2 for v in range(1, 32)}
    t32 = {v: 2 for v in range(1, 33)}
    trace31, trace32 = [], []
    f31, a31 = propagate(31, e31, t31, trace=trace31)
    f32, a32 = propagate(32, e32, t32, trace=trace32)
    s31, trials31 = subset_completions(31, f31, a31, t31)
    s32, nodes32 = obligation_completions(32, e32, t32)
    r31 = sorted((cycle_record(31, s, f31) for s in s31),
                 key=lambda x: x['component_lengths'])
    r32 = sorted((cycle_record(32, s, f32) for s in s32),
                 key=lambda x: x['component_lengths'])
    require((len(e31), len(f31), len(a31), trials31) == (44, 23, 38, 6435),
            'G31 propagation changed')
    require([r['component_lengths'] for r in r31] == [[11, 20], [14, 17]],
            'G31 factor types changed')
    require((len(e32), len(f32), len(a32)) == (46, 23, 42),
            'G32 propagation changed')
    require([r['component_lengths'] for r in r32] ==
            [[3, 6, 23], [3, 29], [5, 27], [9, 23], [32]],
            'G32 factor types changed')

    # Forced components retain every vertex and both distinct endpoint ports.
    chains31 = components(31, f31)
    require(sorted(map(len, chains31)) == [1, 1, 1, 3, 3, 5, 7, 10],
            'wrong forced-chain lengths')
    block_of = {v: i for i, p in enumerate(chains31) for v in p}
    kernel_edges = [{'original_edge': [u, v], 'blocks': [block_of[u], block_of[v]],
                     'ports': [u, v]} for u, v in sorted(a31 - f31)]
    require(all(block_of[u] != block_of[v] for u, v in a31 - f31),
            'unexpected kernel loop')

    # Optional graph is a subdivided K4, with demand 2 at 1,4,15 and 1 elsewhere.
    corridors = [[1, 3], [1, 15], [1, 24, 12, 4], [3, 6, 10, 15],
                 [3, 22, 14, 11, 5, 4], [4, 21, 15]]
    require(set().union(*(edge_set_of_path(p) for p in corridors)) == a31 - f31,
            'K4 corridors do not partition optional edges')
    require(sum(len(p) - 1 for p in corridors) == 15, 'overlapping corridors')
    binary_solutions = []
    binary_edge_sets = []
    for a, b, c, d, e, f in product((0, 1), repeat=6):
        if not (a+b+c == 2 and a+d+e == 1 and c+e+f == 2
                and b+d+1-f == 2):
            continue
        binary_solutions.append([a, b, c, d, e, f])
        indicators = {(1,3):a, (1,15):b, (1,24):c, (4,12):c,
                      (12,24):1-c, (3,6):d, (10,15):d, (6,10):1-d,
                      (3,22):e, (11,14):e, (4,5):e,
                      (14,22):1-e, (5,11):1-e, (4,21):f, (15,21):1-f}
        binary_edge_sets.append(f31 | {edge for edge, chosen in indicators.items()
                                      if chosen})
    require(binary_solutions == [[0,1,1,0,1,0], [0,1,1,1,0,1]],
            'binary proof has different solutions')
    require({tuple(sorted(s)) for s in binary_edge_sets} ==
            {tuple(sorted(s)) for s in s31}, 'binary proof/subset disagreement')

    # A bare quotient does have a spurious spanning cycle: its 5 and 11 ports
    # would each be used twice, while their partner ports 22 and 24 go unused.
    representative_cycle = [1, 3, 4, 5, 11, 10, 6, 15]
    bare_cycle_blocks = [block_of[v] for v in representative_cycle]
    require(len(set(bare_cycle_blocks)) == len(chains31) == 8,
            'hostile does not span the kernel')
    bare_edges = {tuple(sorted((block_of[u], block_of[v]))) for u, v in a31-f31}
    require(edge_set_of_path(bare_cycle_blocks, close=True) <= bare_edges,
            'hostile is not a bare-kernel cycle')

    # The new vertex 32 is precisely the compulsory ear 4--32--17.
    require(e32 - e31 == {(4,32), (17,32)}, 'wrong arrival32 neighbours')
    path_target = dict(t31)
    path_target[4] = path_target[17] = 1
    path_completions, path_nodes = obligation_completions(31, e31, path_target)
    connected_paths = [s for s in path_completions if len(components(31, s)) == 1]
    h32 = [s for s in s32 if len(components(32, s)) == 1]
    require(len(connected_paths) == len(h32) == 1, 'uniqueness failed')
    require(connected_paths[0] | {(4,32), (17,32)} == h32[0],
            'ear bijection failed')
    specified_path = components(31, connected_paths[0])[0]
    require(specified_path[0] == 4 and specified_path[-1] == 17,
            'wrong prescribed endpoints')
    cycle32 = components(32, h32[0])[0]

    # Explicit local surgery from the (11,20) factor to the 32-cycle.
    cover_11_20 = next(s for s in s31 if sorted(map(len, components(31,s))) == [11,20])
    removed, added = cover_11_20 - h32[0], h32[0] - cover_11_20
    require(removed == {(1,24), (4,12), (8,17)}, 'wrong surgery deletion')
    require(added == {(1,8), (12,24), (4,32), (17,32)}, 'wrong surgery addition')

    # Lower orders need no Hamiltonian census: each has a vertex of degree <=1.
    low_witnesses = []
    for n in range(1,31):
        en = edges(n)
        v = min(v for v in range(1,n+1) if len(incident(en,v)) <= 1)
        low_witnesses.append({'n': n, 'vertex': v,
                              'neighbours': sorted(y if x==v else x
                                                   for x,y in incident(en,v))})
    require(all(incident(edges(n),18) == {(7,18)} for n in range(18,31)),
            'leaf18 window changed')
    require(min(len(incident(e31,v)) for v in range(1,32)) == 2,
            'G31 minimum degree not two')

    certificate = {
        'status': 'FINITE-EXACT, with an elementary structural G31 proof',
        'graph': 'Simple undirected G_n on 1..n; distinct vertices adjacent iff their sum is a square.',
        'scope': 'First Hamiltonian cycle is at n=32; this does not prove existence for every n>=32.',
        'g31': {'edges': sorted(e31), 'forced_edges': sorted(f31),
                'allowed_edges': sorted(a31), 'optional_edges': sorted(a31-f31),
                'forcing_trace': trace31, 'subset_trials': trials31,
                'two_factors': r31, 'forced_chains': chains31,
                'port_kernel_edges': kernel_edges, 'optional_K4_corridors': corridors,
                'binary_variable_order': ['a','b','c','d','e','f'],
                'binary_solutions': binary_solutions,
                'bare_kernel_hostile_representatives': representative_cycle},
        'g32': {'edges': sorted(e32), 'forced_edges': sorted(f32),
                'allowed_edges': sorted(a32), 'forcing_trace': trace32,
                'two_factors': r32, 'degree_solver_nodes': nodes32,
                'hamiltonian_cycle': cycle32},
        'ear_bijection': {'new_edges': [[4,32],[17,32]],
                         'g31_prescribed_endpoint_path': specified_path,
                         'degree_feasible_endpoint_subgraphs': len(path_completions),
                         'connected_endpoint_paths': len(connected_paths),
                         'degree_solver_nodes': path_nodes},
        'surgery_from_11_20_cover': {'delete': sorted(removed), 'add': sorted(added)},
        'lower_order_degree_witnesses': low_witnesses,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(args.output).write_text(json.dumps(certificate, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'g31_two_factor_types': [r['component_lengths'] for r in r31],
                      'g32_two_factor_types': [r['component_lengths'] for r in r32],
                      'g31_subset_trials': trials31, 'g32_solver_nodes': nodes32,
                      'prescribed_path_solver_nodes': path_nodes,
                      'output': args.output}, sort_keys=True))


if __name__ == '__main__':
    main()
