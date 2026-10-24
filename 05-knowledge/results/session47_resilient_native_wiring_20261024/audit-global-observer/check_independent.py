#!/usr/bin/env python3
"""Independent global-observer compiler audit. No producer imports."""
import argparse
from collections import Counter
import hashlib
import itertools as it
import json
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
COUNTS = Counter()


def need(value, section, reason):
    COUNTS[section] += 1
    if not value:
        raise RuntimeError(section + ': ' + reason)


def pair(a, b):
    return frozenset((a, b))


def all_sources(n):
    pairs = list(it.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        yield frozenset((a, b) if bits >> k & 1 else (b, a) for k, (a, b) in enumerate(pairs))


def adjacency(path):
    return frozenset(pair(a, b) for a, b in zip(path, path[1:]))


def legal_paths(arcs, n):
    return [p for p in it.permutations(range(n)) if all(edge in arcs for edge in zip(p, p[1:]))]


def literal(legal, required, forbidden, pins):
    result = set()
    for p in legal:
        adj = adjacency(p)
        if required <= adj and not forbidden & adj and all(p.index(v) == clock for v, clock in pins.items()):
            result.add(p)
    return result


def compile_blocks(arcs, n, required, forbidden):
    """Independent weak-component + topological-order implementation."""
    if required & forbidden:
        return None
    marked = {(a, b) for a, b in arcs if pair(a, b) in required}
    incoming, outgoing = Counter(b for a, b in marked), Counter(a for a, b in marked)
    if any(incoming[v] > 1 or outgoing[v] > 1 for v in range(n)):
        return None
    components = [{v} for v in range(n)]
    for a, b in marked:
        left = next(g for g in components if a in g)
        right = next(g for g in components if b in g)
        if left is not right:
            components.remove(left)
            components.remove(right)
            components.append(left | right)
    blocks = []
    for component in components:
        unused, order = set(component), []
        while unused:
            ready = [v for v in unused if not any((u, v) in marked for u in unused)]
            if len(ready) != 1:
                return None
            v = ready[0]
            order.append(v)
            unused.remove(v)
        block = tuple(order)
        if not all(edge in marked for edge in zip(block, block[1:])):
            return None
        blocks.append(block)
    blocks.sort()
    transition = {(i, j) for i, b in enumerate(blocks) for j, c in enumerate(blocks)
                  if i != j and (b[-1], c[0]) in arcs and pair(b[-1], c[0]) not in forbidden}
    return blocks, transition


def expand(compiled, n, pins):
    if compiled is None:
        return set()
    blocks, transition = compiled
    starts = {}
    for i, block in enumerate(blocks):
        requested = {clock - block.index(v) for v, clock in pins.items() if v in block}
        if len(requested) > 1:
            return set()
        if requested:
            start = next(iter(requested))
            if start < 0 or start + len(block) > n:
                return set()
            starts[i] = start
    result = set()
    # Dynamic extension of a block prefix, not permutations as in the producer.
    stack = [((), frozenset(range(len(blocks))), 0)]
    while stack:
        order, available, elapsed = stack.pop()
        if not available:
            result.add(tuple(v for i in order for v in blocks[i]))
            continue
        for i in available:
            if order and (order[-1], i) not in transition:
                continue
            if i in starts and starts[i] != elapsed:
                continue
            stack.append((order + (i,), available - {i}, elapsed + len(blocks[i])))
    return result


def test_compiler():
    rows = []
    pin_cases = converse_cases = contradiction_cases = 0
    for n in range(1, 5):
        pairs = [pair(a, b) for a, b in it.combinations(range(n), 2)]
        cases = 0
        for source_index, arcs in enumerate(all_sources(n)):
            legal = legal_paths(arcs, n)
            for flags in it.product(range(4), repeat=len(pairs)):
                required = frozenset(e for e, f in zip(pairs, flags) if f & 1)
                forbidden = frozenset(e for e, f in zip(pairs, flags) if f & 2)
                direct = literal(legal, required, forbidden, {})
                compiled = compile_blocks(arcs, n, required, forbidden)
                got = expand(compiled, n, {})
                need(got == direct, 'compiler', 'all flags, including contradictions')
                if required & forbidden:
                    contradiction_cases += 1
                # Multiple pins and out-of-range clocks go beyond the producer.
                if source_index % 16 == 0 and cases % 19 == 0:
                    assignments = [{v: c} for v in range(n) for c in range(-1, n + 1)]
                    assignments += [{u: c, v: d} for u, v in it.combinations(range(n), 2)
                                    for c in range(n) for d in range(n)]
                    converse = frozenset((b, a) for a, b in arcs)
                    other = compile_blocks(converse, n, required, forbidden)
                    for pins in assignments:
                        need(expand(compiled, n, pins) == literal(legal, required, forbidden, pins),
                             'pins', 'absolute clocks including conflicting and impossible pins')
                        reflected = {v: n - 1 - clock for v, clock in pins.items()}
                        need(expand(other, n, reflected) == {p[::-1] for p in expand(compiled, n, pins)},
                             'converse', 'reflected absolute clocks and reversed block paths')
                        pin_cases += 1
                        converse_cases += 1
                if source_index % 8 == 0 and not required & forbidden:
                    converse = frozenset((b, a) for a, b in arcs)
                    need(expand(compile_blocks(converse, n, required, forbidden), n, {}) == {p[::-1] for p in direct},
                         'converse', 'unclocked converse')
                    converse_cases += 1
                cases += 1
        rows.append({'owners': n, 'all_tournament_four_flag_instances': cases})
    rng = random.Random(47003)
    extra_cases = 0
    for n in (5, 6):
        pairs = list(it.combinations(range(n), 2))
        for _ in range(256):
            arcs = frozenset((a, b) if rng.randrange(2) else (b, a) for a, b in pairs)
            flags = [rng.randrange(3) for _ in pairs]
            required = frozenset(pair(a, b) for (a, b), f in zip(pairs, flags) if f == 1)
            forbidden = frozenset(pair(a, b) for (a, b), f in zip(pairs, flags) if f == 2)
            pins = {v: rng.randrange(-1, n + 1) for v in rng.sample(range(n), 2)}
            compiled = compile_blocks(arcs, n, required, forbidden)
            need(expand(compiled, n, pins) == literal(legal_paths(arcs, n), required, forbidden, pins),
                 'compiler', 'additional five/six-owner two-pin cases')
            extra_cases += 1
    return {'rows': rows, 'contradictory_flag_instances': contradiction_cases,
            'global_pin_instances': pin_cases, 'converse_instances': converse_cases,
            'extra_five_six_owner_cases': extra_cases, 'extra_seed': 47003}


def insertion(arcs, order):
    path = []
    for v in order:
        index = next((i for i, w in enumerate(path) if (v, w) in arcs), len(path))
        path.insert(index, v)
    return tuple(path)


def semicomplete_controls():
    graph_count = insertion_orders = 0
    for n in range(1, 5):
        pairs = list(it.combinations(range(n), 2))
        for states in it.product((0, 1, 2), repeat=len(pairs)):
            arcs = set()
            for (a, b), state in zip(pairs, states):
                if state in (0, 2):
                    arcs.add((a, b))
                if state in (1, 2):
                    arcs.add((b, a))
            for order in it.permutations(range(n)):
                path = insertion(arcs, order)
                need(all(edge in arcs for edge in zip(path, path[1:])),
                     'semicomplete', 'insertion constructs a Hamiltonian path')
                insertion_orders += 1
            graph_count += 1
    # A linear forest need not contract to a tournament or even a semicomplete digraph.
    required = frozenset((pair(0, 1), pair(2, 3)))
    neither = frozenset(((0, 1), (2, 3), (2, 1), (0, 3), (0, 2), (1, 3)))
    both = frozenset(((0, 1), (2, 3), (1, 2), (3, 0), (0, 2), (1, 3)))
    cn, cb = (compile_blocks(a, 4, required, frozenset()) for a in (neither, both))
    need(cn is not None and len(cn[1]) == 0 and expand(cn, 4, {}) == set(),
         'hostiles', 'linear forest with neither block transition')
    need(cb is not None and len(cb[1]) == 2 and len(expand(cb, 4, {})) == 2,
         'hostiles', 'linear forest with both block transitions')
    transitive = frozenset(it.combinations(range(3), 2))
    c = compile_blocks(transitive, 3, frozenset(), frozenset())
    need(expand(c, 3, {}) == {(0, 1, 2)} and expand(c, 3, {2: 0}) == set(),
         'hostiles', 'semicomplete sufficiency does not erase absolute clocks')
    c = compile_blocks(transitive, 3, frozenset(), frozenset((pair(0, 2),)))
    need(expand(c, 3, {}) == {(0, 1, 2)} and len(c[1]) == 2,
         'hostiles', 'semicompleteness is sufficient, not necessary')
    cycle = frozenset(((0, 1), (1, 2), (2, 0)))
    need(compile_blocks(cycle, 3, frozenset(pair(a, b) for a, b in cycle), frozenset()) is None,
         'hostiles', 'required directed cycle is rejected')
    return {'all_semicomplete_graphs_through_four': graph_count,
            'all_insertion_orders': insertion_orders}


def native_controls():
    raw = (HERE / 'inputs/shared_owner.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest() == '5266665af6f013ff54a64a5b4b7ef83ac35a8b80c01507ef05cc9ed26d1b6880',
         'native', 'actual inherited source digest')
    data = json.loads(raw)
    host = data['two_cell_host_adjacency']
    arcs = frozenset((i, j) for i, row in enumerate(host) for j, value in enumerate(row) if value)
    need(len(host) == 14 and all(host[i][i] == 0 for i in range(14)) and
         all(host[i][j] + host[j][i] == 1 for i in range(14) for j in range(i + 1, 14)),
         'native', 'full named source host is a tournament')
    records = data['two_cell_records']
    for label in ('176', '368'):
        embedding = data['embeddings'][label]
        for row in records[label]:
            p = tuple(row['global_path'])
            need(p == tuple(embedding[v] for v in row['local_path']), 'native', 'literal owner embedding')
            need(len(set(p)) == 8 and set(p) == set(embedding), 'native', 'local Hamiltonian owner census')
            need(all(edge in arcs for edge in zip(p, p[1:])), 'native', 'local path arcs in complete host')
    allowed = [(r['core_bit'], r['helper_bit']) for r in data['all_bit_pairs'] if r['allowed']]
    need(allowed == [(0, 1), (1, 0)], 'native', 'exact allowed local states')
    obstructions = []
    for c, h in allowed:
        left = tuple(records['176'][c]['global_path'])
        right = tuple(records['368'][h]['global_path'])
        required_arcs = set(zip(left, left[1:])) | set(zip(right, right[1:]))
        branch = 2 if c == 0 else 6
        need({b for a, b in required_arcs if a == branch} == {3, 12}, 'native', 'literal successor obstruction')
        required = frozenset(pair(a, b) for a, b in required_arcs)
        need(compile_blocks(arcs, 14, required, frozenset()) is None,
             'native', 'no global observer independently of missing negative flags')
        # Shared-pair partner flags still agree. This verifies the positive
        # local-host claim is not what the branch obstruction disproves.
        shared = set(left) & set(right)
        need(shared == {2, 6} and all((e in adjacency(left)) == (e in adjacency(right))
                                     for e in (pair(a, b) for a, b in it.combinations(shared, 2))),
             'native', 'compatible shared-pair flags coexist with global observer failure')
        obstructions.append({'core_bit': c, 'helper_bit': h, 'branch': branch, 'successors': [3, 12]})
    return {'source_sha256': hashlib.sha256(raw).hexdigest(), 'obstructions': obstructions,
            'source_host_order': 14, 'local_paths_authenticated': 4}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    pins = json.loads((HERE / 'INPUT_PINS.json').read_text())
    for row in pins['inputs']:
        p = HERE / row['snapshot']
        need(p.stat().st_size == row['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest() == row['sha256'],
             'authentication', row['source'])
    result = {'status': 'PASS', 'python_optimization': sys.flags.optimize,
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'input_pins_sha256': hashlib.sha256((HERE / 'INPUT_PINS.json').read_bytes()).hexdigest()}
    result['compiler'] = test_compiler()
    result['semicomplete'] = semicomplete_controls()
    result['native'] = native_controls()
    result['checks_by_section'] = dict(sorted(COUNTS.items()))
    result['explicit_checks'] = sum(COUNTS.values())
    result['elapsed_seconds'] = round(time.monotonic() - started, 3)
    output = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
