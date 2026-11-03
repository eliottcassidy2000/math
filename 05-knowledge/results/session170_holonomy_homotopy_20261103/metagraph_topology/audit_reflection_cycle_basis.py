#!/usr/bin/env python3
"""Independent signed action on actual complement-graph cycle bases.

Uses check.py only for complete tournament isomorphism orbits and comparison
data. The H1 trace is obtained from fundamental cycles, not the fixed-cell
character formula. All directions are physical masks: bit 1 means i -> j.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import json
import sys


def audit(checker, n):
    pairs, canonical, representatives = checker.tournament_orbits(n)
    class_index = {code: i for i, code in enumerate(representatives)}
    chords = [pairs.index((i, j)) for i, j in pairs if j >= i + 2]
    chord_pairs = [pairs[i] for i in chords]
    all_arcs = (1 << len(pairs)) - 1
    full = (1 << len(chords)) - 1
    # S151 flag bit 1 reverses a chord from the forward transitive tournament.
    flags = [class_index[canonical[all_arcs ^ sum(
        ((t >> i) & 1) << arc for i, arc in enumerate(chords))]]
        for t in range(full + 1)]
    lines = [(t, t ^ full) for t in range(full + 1) if t < (t ^ full)]
    coordinate_reflection = [chord_pairs.index((n - 1 - j, n - 1 - i))
                             for i, j in chord_pairs]

    def reflect(t):
        return sum(((t >> i) & 1) << target
                   for i, target in enumerate(coordinate_reflection))

    # Each physical line is oriented from its smaller numerical flag.
    line_of_flag = {t: (e, 1) for e, (t, u) in enumerate(lines)}
    line_of_flag.update({u: (e, -1) for e, (t, u) in enumerate(lines)})
    reflected_edges = [line_of_flag[reflect(t)] for t, u in lines]
    parent = list(range(len(representatives)))
    adjacency = [[] for _ in representatives]
    tree_edges = set()

    def root(v):
        while v != parent[v]:
            v = parent[v]
        return v

    for e, (t, u) in enumerate(lines):
        a, b = flags[t], flags[u]
        if root(a) != root(b):
            parent[root(a)] = root(b)
            tree_edges.add(e)
            adjacency[a].append((b, e, 1))
            adjacency[b].append((a, e, -1))

    trace = 0
    for e, (t, u) in enumerate(lines):
        if e in tree_edges:
            continue
        a, b = flags[t], flags[u]
        # The fundamental cycle traverses e from a to b, then tree b to a.
        pending = [(b, {e: 1})]
        seen = {b}
        cycle = None
        while pending:
            v, candidate = pending.pop()
            if v == a:
                cycle = candidate
                break
            for w, tree_edge, sign in adjacency[v]:
                if w not in seen:
                    seen.add(w)
                    pending.append((w, candidate | {tree_edge: sign}))
        if cycle is None:
            raise RuntimeError('Failed to close a non-tree edge in its component')
        moved = {}
        for edge, coefficient in cycle.items():
            target, sign = reflected_edges[edge]
            moved[target] = moved.get(target, 0) + coefficient * sign
        # In a fundamental-cycle basis, coefficient at non-tree e is its
        # own basis coefficient. Thus this is the diagonal matrix entry.
        trace += moved.get(e, 0)

    expected_dimensions = checker.corner_test(n)['reflection']['graph_eigenspace_dimensions']
    if trace != expected_dimensions[0] - expected_dimensions[1]:
        raise RuntimeError(f'Reflection trace mismatch at n={n}')
    return {'n': n, 'cycle_rank': len(lines) - len(tree_edges),
            'direct_signed_h1_trace': trace}


if __name__ == '__main__':
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('check.py')
    spec = spec_from_file_location('s170_checker', path)
    checker = module_from_spec(spec)
    spec.loader.exec_module(checker)
    print(json.dumps([audit(checker, n) for n in range(3, 7)], indent=2))
