#!/usr/bin/env python3
"""Complete affine review of inherited endpoint-four base cases at n=12,13,14.

No tournament or Hamiltonian-path enumeration is used. The two imported
independent-audit files are pinned below; every other control is local.
"""
from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent / 'endpoint_four_repairs_20261021' / 'audit'
PINS = {
    'independent_audit.py': '376782e9c9eab3f3dae563a22735271d064fe73f4ef0361de044beb56cdc85b0',
    'audit_core.py': '9e0f86d30773b46dad359c04ead55dc07e606cad54156fb4eb0025b42a7339a4',
}

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

for name, digest in PINS.items():
    need(sha256((AUDIT / name).read_bytes()).hexdigest() == digest,
         ('inherited audit source changed', name))
sys.path.insert(0, str(AUDIT))
spec = spec_from_file_location('inherited_endpoint_audit',
                              AUDIT / 'independent_audit.py')
bank = module_from_spec(spec)
old_argv = sys.argv
try:
    sys.argv = [str(AUDIT / 'independent_audit.py')]
    spec.loader.exec_module(bank)
finally:
    sys.argv = old_argv


def edge_set(q):
    return {frozenset((u, v)) for u, v in zip(q, q[1:])}


def literal_arc(rows, u, v):
    return (rows[u] >> v) & 1


def literal_directed(rows, q):
    return all(literal_arc(rows, u, v) for u, v in zip(q, q[1:]))


def literal_blue(rows, q):
    r = dict(zip(q, q[::-1]))
    return all(literal_arc(rows, u, v) ==
               literal_arc(rows, r[v], r[u])
               for u in range(len(rows)) for v in range(len(rows)) if u != v)


def literal_k(rows, q):
    edges = edge_set(q)
    return tuple(sum((literal_arc(rows, u, v) ^
                      int(frozenset((u, v)) not in edges)) << v
                     for v in range(len(rows)) if u != v)
                 for u in range(len(rows)))


def literal_check(rows, ports, opened, chord, words, expected):
    n = len(rows)
    P = list(range(n))
    need(all(rows[u] >> u & 1 == 0 for u in range(n)), 'literal loops')
    need(all(literal_arc(rows, u, v) ^ literal_arc(rows, v, u) == 1
             for u in range(n) for v in range(u+1, n)), 'literal tournament')
    need(literal_directed(rows, P) and literal_blue(rows, P)
         and literal_arc(rows, n-1, 0), 'literal closed blue source')
    need(rows[0].bit_count() == n-5 and rows[-1].bit_count() == 4,
         'literal endpoint scores')
    need({v for v in range(1, n-1) if literal_arc(rows, v, 0)} == set(ports),
         'literal port set')
    need(literal_arc(rows, 1, n-2) == opened
         and literal_arc(rows, 2, n-2) == chord, 'literal branch')
    target = literal_k(rows, P)
    xor = set()
    blue = []
    for q in words:
        need(sorted(q) == P and literal_directed(rows, q),
             'literal directed Hamiltonian word')
        blue.append(literal_blue(rows, q))
        xor ^= edge_set(q)
        rows = literal_k(rows, q)
    need(xor == edge_set(P) and rows == target, 'literal target and XOR')
    need(blue == expected, ('literal blue flags', blue, expected))
    return blue


def affine_model(source, equations):
    basis = bank.eliminate(equations)
    need(basis is not None, 'blue hostile must be feasible')
    assignment = 0
    for pivot, equation in sorted(basis.items()):
        bit = ((equation & 1) ^
               ((equation & ~(1 << pivot) & assignment).bit_count() & 1))
        if bit:
            assignment |= 1 << pivot
    need(all((e & 1) ^ ((e & assignment).bit_count() & 1) == 0
             for e in equations), 'blue model solves every equation')
    return tuple(sum(((source[u][v] & 1) ^
                      ((source[u][v] & assignment).bit_count() & 1)) << v
                     for v in range(len(source)) if u != v)
                 for u in range(len(source)))


def main():
    totals = {}
    hostiles = []
    for n in (12, 13, 14):
        cases = Counter()
        tested = skipped = feasible = 0
        for ports in combinations(range(2, n-1), 3):
            I = set(ports)
            for opened, chord in ((True, False), (True, True),
                                  (False, False), (False, True)):
                # Precisely the residual left by base cases 1--3.
                if not opened and {n-3, n-2} <= I:
                    skipped += 1
                    continue
                name, words = bank.words(n, I, opened, chord)
                source = bank.symbolic_source(
                    n, I, [(1, n-2, int(opened)), (2, n-2, int(chord))])
                # symbolic_source explicitly rejects an inconsistent source.
                feasible += 1
                current = source
                edge_xor = set()
                black = []
                blues = []
                for index, q in enumerate(words):
                    need(sorted(q) == list(range(n)),
                         ('word permutation', n, ports, name, index))
                    need(all(current[u][v] == 1 for u, v in zip(q, q[1:])),
                         ('affine directedness', n, ports, name, index))
                    equations = bank.blue_equations(current, q)
                    black.append(bank.eliminate(equations) is None)
                    blues.append(equations)
                    current = bank.k_move(current, q)
                    edge_xor ^= edge_set(q)
                need(edge_xor == edge_set(list(range(n))), 'affine path XOR')
                need(current == bank.k_move(source, list(range(n))),
                     'affine literal target')
                tested += 1
                cases[(name, tuple(black))] += 1
                if not all(black):
                    bad = black.index(False)
                    rows = affine_model(source, blues[bad])
                    literal_check(rows, ports, opened, chord, words,
                                  [not value for value in black])
                    hostiles.append({
                        'n': n, 'ports': ports, 'middle_open': opened,
                        'x_to_c': chord, 'branch': name, 'source_rows': rows,
                        'blue_step': bad+1, 'words': words,
                        'directed_steps': [True, True, True],
                        'literal_blue_flags': [not value for value in black],
                        'literal_target_and_XOR': True,
                    })
        need(feasible == tested, 'no vacuous source cell')
        totals[n] = {
            'all_source_branch_cells': tested+skipped,
            'base_source_cells': tested,
            'closed_endpair_cells_deferred_to_proved_repairs': skipped,
            'branches': [{'branch': name, 'all_completion_black': black,
                          'cells': count}
                         for (name, black), count in sorted(cases.items())],
        }
    need([totals[n]['base_source_cells'] for n in (12, 13, 14)] ==
         [322, 464, 642], 'complete base-cell counts')
    expected = {
        (12, (2, 7, 9), True, True, 1),
        (12, (6, 7, 9), True, True, 1),
        (14, (3, 5, 10), False, True, 3),
        (14, (3, 7, 10), False, True, 3),
        (14, (5, 9, 10), False, False, 3),
        (14, (7, 9, 10), False, False, 3),
    }
    need({(h['n'], h['ports'], h['middle_open'], h['x_to_c'], h['blue_step'])
          for h in hostiles} == expected, 'exact failed-ledger cells')
    boundary = [bank.boundary(n, branch)
                for n in (12, 13, 14) for branch in ('open', 'closed')]
    need(all(r['inconsistent'] or r['minimum_middle_inports'] >= 4
             for r in boundary), 'second-word port budget')
    need(not any(h['n'] == 13 for h in hostiles), 'order13 full base closure')
    print(json.dumps({
        'status': 'PASS',
        'method': 'complete affine source cells; pinned independent GF2 solver; local literal hostile replay',
        'dependency_sha256': PINS,
        'orders': totals, 'second_word_boundary': boundary,
        'failed_ledger_hostiles': hostiles,
        'scope': 'order13 complete base closure FINITE-EXACT; order12/14 failures are ledger failures, not connectivity obstructions',
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

