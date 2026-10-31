#!/usr/bin/env python3
"""Discover complete final lifts preserving a native path and the inherited seed.

The solver is a discovery dependency only. verify_native_forest.py consumes the
saved words without OR-Tools, independently checking all universal port states.
"""
import argparse
import json
from math import isqrt
from pathlib import Path
from ortools.sat.python import cp_model

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'seeded_final_split_table_20261019/final_split_glue_certificates_20261019.json'


def endpoints(state):
    p, d = state['parity'], state['residue']
    chosen = int(d >= 57)
    pairs = {}
    for c in range(-24, 25):
        i = int(49 + (-1) ** (p + 1) * c <= d)
        end = 10 if i == chosen else 8
        a, b = (98-c, 343+c) if state['order'+str(i)] == 1 else (343+c, 98-c)
        pairs[(c, i, 0)] = (3381-c, a)
        pairs[(c, i, 1)] = (b, 49*end + (-1)**(p+i)*c)
    return pairs


def solve(state, pattern):
    pairs = endpoints(state)
    nodes = list(range(1, 25)) + [x for pair in pairs.values() for x in pair]
    if len(nodes) != len(set(nodes)):
        raise RuntimeError('port collision')
    indexes = {v: i for i, v in enumerate(nodes)}
    owner = {x: (key, pair) for key, pair in pairs.items() for x in pair}
    internals = {frozenset(pair) for pair in pairs.values()}
    model = cp_model.CpModel()
    arcs = []
    variables = {}
    for ia, a in enumerate(nodes):
        for ib, b in enumerate(nodes):
            if ia == ib:
                continue
            if frozenset((a, b)) in internals or isqrt(a+b)**2 == a+b:
                variable = model.NewBoolVar(f'e{a}_{b}')
                arcs.append((ia, ib, variable))
                variables[a, b] = variable
    model.AddCircuit(arcs)
    for a, b in pairs.values():
        model.Add(variables[a, b] + variables[b, a] == 1)
    for a, b in list(zip(pattern, pattern[1:])) + [(482, 3362)]:
        model.Add(variables[a, b] + variables[b, a] == 1)
    # Original exact cycle is only a hint; the requested new edge may invalidate it.
    old = []
    for atom in state['word']:
        if atom[0] == 'single':
            old.append(atom[1])
        else:
            _, c, i, segment, orientation = atom
            pair = pairs[c, i, segment]
            old.extend(pair if orientation == 1 else pair[::-1])
    old_arcs = set(zip(old, old[1:]+old[:1]))
    for edge, variable in variables.items():
        model.AddHint(variable, int(edge in old_arcs))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 20261031
    solver.parameters.max_time_in_seconds = 15
    status = solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None, solver.StatusName(status)
    nxt = {a: b for (a, b), variable in variables.items() if solver.Value(variable)}
    word, seen = [], set()
    vertex = 1
    while vertex not in seen:
        seen.add(vertex)
        if vertex <= 24:
            word.append(['single', vertex])
            vertex = nxt[vertex]
        else:
            key, pair = owner[vertex]
            other = pair[1] if vertex == pair[0] else pair[0]
            if nxt[vertex] != other:
                raise RuntimeError('internal segment was not traversed')
            seen.add(other)
            word.append(['segment', *key, 1 if vertex == pair[0] else -1])
            vertex = nxt[other]
    if vertex != 1 or len(seen) != 220:
        raise RuntimeError('not a spanning port cycle')
    if word[1] == ['single', 8]:
        def reverse(a):
            return a if a[0] == 'single' else a[:-1]+[-a[-1]]
        word = word[:1]+[reverse(a) for a in word[:0:-1]]
    if word[-2:] != [['single', 17], ['single', 8]]:
        raise RuntimeError('native path normalization failed')
    return word, 'FEASIBLE'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pattern', choices=('short', 'n24'), default='n24')
    args = parser.parse_args()
    pattern = [1, 8, 17] if args.pattern == 'short' else [1,8,17,19,6,10,15,21,4,12,24]
    stem = 'native_forest' if args.pattern == 'short' else 'n24_forest'
    source = json.loads(SOURCE.read_text())
    output = []
    for index, state in enumerate(source):
        word, status = solve(state, pattern)
        if word is None:
            (HERE / (stem+'_partial.json')).write_text(json.dumps(output))
            raise RuntimeError(('incomplete search', index, status,
                                {k: state[k] for k in ('parity','residue','order0','order1')}))
        output.append({**{k: state[k] for k in ('parity','residue','order0','order1')},
                       'word': word})
        if len(output) % 20 == 0:
            print(json.dumps({'certified_words_found': len(output)}), flush=True)
    (HERE / (stem+'_table.json')).write_text(json.dumps(output, separators=(',', ':'))+'\n')
    print(json.dumps({'status': 'discovery complete', 'states': len(output)}), flush=True)


if __name__ == '__main__':
    main()
