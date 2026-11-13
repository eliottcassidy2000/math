#!/usr/bin/env python3
"""Independent small audit plus differential testing of the producer cover solver.

The concrete witnesses and partition universal property are rebuilt here.
Only the differential small-family test loads the producer exact_cover function;
its reference answer enumerates every column subset. No producer main is run.
No writes occur; deterministic JSON goes to stdout.
"""
import ast
from fractions import Fraction
from hashlib import sha256
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
CHECKS = 0


def need(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(msg)


def union(values):
    out = 0
    for value in values:
        out |= value
    return out


def brute_cover(masks):
    target = union(masks)
    for k in range(len(masks)+1):
        for cc in combinations(range(len(masks)), k):
            if union(masks[i] for i in cc) == target:
                return k
    raise RuntimeError("finite universe not covered by its own columns")


def differential_cover():
    path = PACKET/'closure/analyze.py'
    spec = importlib.util.spec_from_file_location('audited_cover_producer', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    families = 0
    for size in range(6):
        for masks in combinations(range(1, 16), size):
            expected = brute_cover(masks)
            result = module.exact_cover(masks)
            need(result['minimum'] == expected, "differential cover optimum")
            need(union(masks[i] for i in result['columns']) == union(masks),
                 "differential witness")
            families += 1
    # Zeros, repeated columns, column order, holes in bit indices, strong dominance.
    extra = [(), (0,), (0, 0), (1, 1), (0, 3, 1, 2, 3),
             (2, 8, 10), (1, 2, 4, 8, 15), (5, 10, 3, 12),
             (12, 3, 10, 5), (1 << 17, (1 << 17) | 4, 4)]
    for masks in extra:
        need(module.exact_cover(masks)['minimum'] == brute_cover(masks),
             "differential edge case")
    return {'distinct_nonzero_four_bit_families': families,
            'extra_cases': len(extra), 'reference': 'exhaustive column combinations',
            'producer_function_loaded_only_for_differential_test': True}


def concrete_witnesses(data, analysis):
    rows = []
    for row, claimed in zip(data['rows'], analysis['rows']):
        need(row['n'] == claimed['n'], 'aligned universe')
        nodes = row['merged_nodes']
        codes = [x['code'] for x in nodes]
        index = {c: i for i, c in enumerate(codes)}
        n = len(nodes)
        A = [[0]*n for _ in range(n)]
        B = [[0]*n for _ in range(n)]
        for i, j, value in row['merged_A']:
            A[index[i]][index[j]] = value
        for i, j, value in row['merged_blue']:
            B[index[i]][index[j]] = value
        d = list(map(sum, A))
        need(d == [x['degree'] for x in nodes], 'source degree from actual row')
        need(all(d), 'strictly positive normalization mass')
        # Use ordered covers rather than the producer monomial Q vector.
        pairs = [(i, j) for i, j in combinations(range(n), 2)
                 if nodes[i]['ordered_covers'] == nodes[j]['ordered_covers']]
        need({(codes[i], codes[j]) for i, j in pairs} ==
             {tuple(p) for p in claimed['Q_collision_pairs']}, 'full slot equivalence')
        outputs = {}
        for name in ('A', 'P', 'A_colored', 'P_colored'):
            normalized = name.startswith('P')
            colored = name.endswith('colored')
            def observed(i, c):
                div = d[i] if normalized else 1
                if colored:
                    return (Fraction(B[i][c], div), Fraction(A[i][c]-B[i][c], div))
                return (Fraction(A[i][c], div),)
            masks = [sum(1 << k for k, (i, j) in enumerate(pairs)
                         if observed(i, c) != observed(j, c)) for c in range(n)]
            need(union(masks) == (1 << len(pairs))-1, 'complete literal columns separate')
            certificate = claimed['anchor_optimization'][name]
            cols = [index[c] for c in certificate['columns']]
            need(len(cols) == len(set(cols)) == certificate['minimum'], 'distinct anchor witness')
            need(union(masks[c] for c in cols) == union(masks), 'anchor witness covers every pair')
            # Validate packing with ALL original columns, independently of dominance reduction.
            neighborhoods = []
            for ci, cj in certificate['packing']:
                i, j = index[ci], index[cj]
                candidates = {c for c in range(n) if observed(i, c) != observed(j, c)}
                need(bool(candidates), 'packing demand is satisfiable')
                need(all(not (candidates & old) for old in neighborhoods), 'disjoint anchor neighborhoods')
                neighborhoods.append(candidates)
            if row['n'] <= 6:
                need(brute_cover(masks) == certificate['minimum'], 'small-order exhaustive minimum')
            else:
                need([x['budget'] for x in certificate['failed_budgets']] == [5, 6, 7, 8],
                     'all budgets between packing bound and witness were recorded')
                need(len(neighborhoods) == 5 and certificate['minimum'] == 9,
                     'n7 lower and upper certificates')
            outputs[name] = {'minimum_claimed': certificate['minimum'],
                             'witness_verified': True,
                             'packing_lower_bound': len(neighborhoods),
                             'brute_force_lower_bound_here': row['n'] <= 6}
        rows.append({'n': row['n'], 'same_slot_pairs': len(pairs), 'palettes': outputs})
    return rows


def partitions(n):
    def visit(labels):
        if len(labels) == n:
            return [tuple(tuple(i for i, x in enumerate(labels) if x == k)
                          for k in range(max(labels, default=-1)+1))]
        out = []
        for k in range(max(labels, default=-1)+2):
            out.extend(visit(labels+[k]))
        return out
    return visit([])


def refines(fine, coarse):
    return all(any(set(c) <= set(d) for d in coarse) for c in fine)


def stable(cells, matrices):
    return all(len({tuple(sum(M[i][j] for j in target) for M in matrices
                          for target in cells) for i in source}) <= 1 for source in cells)


def close(cells, matrices):
    trace = [cells]
    while True:
        groups = []
        for source in cells:
            buckets = {}
            for i in source:
                key = tuple(sum(M[i][j] for j in target) for target in cells for M in matrices)
                buckets.setdefault(key, []).append(i)
            groups.extend(tuple(x) for x in buckets.values())
        result = tuple(sorted(groups))
        if result == cells:
            return trace
        cells = result
        trace.append(cells)


def universal_property_controls():
    cases = 0
    for n in range(1, 5):
        pp = partitions(n)
        edges = list(combinations(range(n), 2))
        for bits in product((0, 1), repeat=len(edges)):
            A = [[0]*n for _ in range(n)]
            B = [[0]*n for _ in range(n)]
            for (i, j), value in zip(edges, bits):
                A[i][j] = A[j][i] = value
                B[i][j] = B[j][i] = value * ((i+j) % 2)
            K = [[A[i][j]-B[i][j] for j in range(n)] for i in range(n)]
            for matrices in ([A], [B, K]):
                stable_partitions = [theta for theta in pp if stable(theta, matrices)]
                for initial in pp:
                    trace = close(initial, matrices)
                    result = trace[-1]
                    need(stable(result, matrices) and refines(result, initial), 'closure is stable refinement')
                    need(len(trace)-1 <= n-len(initial), 'strict-round finite bound')
                    need(all(refines(theta, result) for theta in stable_partitions
                             if refines(theta, initial)), 'coarsest stable universal property')
                    cases += 1
    return {'cases': cases, 'universe': 'all simple undirected graphs on 1..4 labelled vertices; all initial partitions; raw and endpoint-parity two-color palettes'}


def pins():
    names = ['closure/analyze.py', 'closure/check.py', 'closure/orbits.cpp',
             'closure/PROOF.md', 'closure/data.json', 'closure/analysis.json',
             'audit/independent.py']
    out = {name: sha256((PACKET/name).read_bytes()).hexdigest() for name in names}
    text = (PACKET/'audit/independent.py').read_text()
    tree = ast.parse(text)
    out['reviewed_root_functions'] = {
        node.name: sha256(ast.get_source_segment(text, node).encode()).hexdigest()
        for node in tree.body if isinstance(node, ast.FunctionDef)
        and node.name in ('literal', 'refine', 'initial')}
    out['audit_closure/check.py'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return out


def main():
    before = pins()
    data = json.loads((PACKET/'closure/data.json').read_text())
    analysis = json.loads((PACKET/'closure/analysis.json').read_text())
    result = {'status': 'FINITE-EXACT; audit controls, not independent n7 cover exhaustion',
              'differential_cover': differential_cover(),
              'concrete_anchor_certificates': concrete_witnesses(data, analysis),
              'partition_universal_property': universal_property_controls()}
    after = pins()
    need(all(before[name] == after[name] for name in before if name != 'audit/independent.py'),
         'reviewed producer files unchanged during checks')
    result['source_pins'] = before
    result['checks'] = CHECKS
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
