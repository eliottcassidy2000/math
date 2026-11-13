"""Literal controls for the sharp common-inflation gap and ready root."""
from pathlib import Path
from itertools import combinations, permutations
import json
from wrapper_repair import birth, legal, require, state, strong

HERE = Path(__file__).resolve().parent
T = (22, 8, 2, 5, 14)
U = (22, 8, 2, 21, 6)


def inflate(rows, sizes):
    groups = []
    offset = 0
    for size in sizes:
        groups.append(tuple(range(offset, offset + size)))
        offset += size
    out = [0] * offset
    for i, group in enumerate(groups):
        for p, vertex in enumerate(group):
            out[vertex] |= sum(1 << other for other in group[p + 1:])
            for j, other_group in enumerate(groups):
                if rows[i] >> j & 1:
                    out[vertex] |= sum(1 << other for other in other_group)
    return tuple(out), groups


def direct_ready(a, r):
    # Parent has quotient A -> u -> D -> A.  Append a child at u and
    # retain that newborn, the first vertex of the new block (newborn,D).
    require(a >= 1 and r >= 2, 'positive two-width construction required')
    parent, groups = inflate((2, 4, 1), (a, 1, r - 1))
    return birth(parent, groups[1][0], True)


def substitute_one(rows, vertex, module):
    sizes = [1] * len(rows)
    sizes[vertex] = len(module)
    out, groups = inflate(rows, sizes)
    out = list(out)
    group = groups[vertex]
    mask = sum(1 << i for i in group)
    for i, old in enumerate(group):
        out[old] &= ~mask
        out[old] |= sum(1 << new for j, new in enumerate(group) if module[i] >> j & 1)
    return tuple(out), groups


def run_coefficients(rows, block):
    # Every module run is replaced by one occurrence of the same marker.
    # Nonadjacent marker occurrences are the ordered distinct module runs.
    one = sum(legal(rows, p) for p in permutations(range(len(rows))))
    outside = [i for i in range(len(rows)) if i != block]
    two = 0
    for p in permutations(outside):
        for slots in combinations(range(len(rows) + 1), 2):
            iterator = iter(p)
            word = tuple(block if i in slots else next(iterator)
                         for i in range(len(rows) + 1))
            two += legal(rows, word)
    return one, two


def main():
    controls = []
    selected = []
    for a in range(1, 5):
        for b in range(1, 5):
            for c in range(1, 5):
                A, B, C = 2**a, 2**b, 2**c
                left, groups = inflate(T, (1, a, b, 1, c))
                right, _ = inflate(U, (1, a, b, 1, c))
                p, q = state(left, 0), state(right, 0)
                qT = (A - 1) * (B - 1) + 2
                hT = (A - 1) * (B * (C + 1) - 1) + C + 2
                qU = (A - 1) * (B * C - 1) + 2
                hU = (A - 1) * (2 * B * C - 1) + 2
                require(p == (hT, qT, qT, hT), 'old inflated profile failed')
                require(q == (hU, qU - 1, qU, hU + 1), 'ready inflated profile failed')
                gap = B * (A - 2) * (C - 1) + (B - 2) * (C - 1) + C - 2
                require(hU - hT == gap and ((gap == 0) == (a == b == c == 1)),
                        'sharp common-inflation gap failed')
                incoming = [i for i, row in enumerate(left) if row & 1]
                require(incoming == [groups[3][0]] and
                        [i for i, row in enumerate(right) if row & 1] == incoming,
                        'physical incoming port changed')
                require(strong(left) and strong(right), 'inflated source is not strong')
                control = {'a': a, 'b': b, 'c': c, 'old_profile': p,
                           'new_profile': q, 'H_gap': gap, 'incoming_port': incoming[0]}
                controls.append(control)
                if (a, b, c) in ((1, 1, 1), (2, 1, 1), (2, 3, 4)):
                    selected.append({**control, 'old_rows': left, 'new_rows': right})
    ready, root = direct_ready(4, 2)
    require(state(ready, root) == (47, 16, 17, 48), 'H47 ready chart failed')
    exit_graph, exit_root = birth(ready, root, True)
    exit_graph, exit_root = birth(exit_graph, exit_root)
    require(state(exit_graph, exit_root) == (143, 47, 47, 143), 'H47 NR exit failed')
    repaired29 = (24, 1, 11, 34, 46, 7)
    parent29 = (8, 1, 11, 18, 7)
    h_coeff = run_coefficients(repaired29, 1)
    q_coeff = run_coefficients(parent29, 1)
    raw_source_coeff = run_coefficients(U, 2)
    require(h_coeff == (29, 15) and q_coeff == (9, 3), 'H29 skeleton coefficients failed')
    require(raw_source_coeff == (9, 4), 'raw U9 skeleton coefficients failed')
    boundary_controls = []
    for b in (1, 2, 3):
        graph, groups = inflate(repaired29, (1, b, 1, 1, 1, 1))
        root29 = groups[4][0]
        z29 = state(graph, root29)
        h29, q29 = 15 * 2**b - 1, 3 * 2**b + 3
        require(z29 == (h29, q29 - 1, q29, h29 + 1), 'H29 inflated ready state failed')
        boundary_controls.append({'b': b, 'rows': graph, 'root': root29, 'profile': z29,
                                  'H_drift_from_raw_NR': 3 * 2**b - 6})
    triangle = (2, 4, 1)
    triangle_modified, triangle_groups = substitute_one(repaired29, 1, triangle)
    modified_state = state(triangle_modified, triangle_groups[4][0])
    require(modified_state == (177, 42, 45, 180), 'arbitrary-module modified profile failed')
    triangle_raw, raw_groups = substitute_one(U, 2, triangle)
    triangle_raw, raw_root = birth(triangle_raw, raw_groups[0][0], True)
    triangle_raw, raw_root = birth(triangle_raw, raw_root)
    raw_state = state(triangle_raw, raw_root)
    require(raw_state == (159, 51, 51, 159), 'arbitrary-module raw profile failed')
    out = {'status': 'FINITE-EXACT controls; universal proof in H9_INFLATION.md',
           'universe': 'Common TT_a,TT_b,TT_c inflation at1,2,4; all positive sizes<=4',
           'controls': controls, 'literal_selected_sources': selected,
           'direct_H47': {'rows': ready, 'root': root, 'profile': [47, 16, 17, 48]},
           'H47_NR_exit': {'rows': exit_graph, 'root': exit_root,
                          'profile': [143, 47, 47, 143]},
           'H29_handoff_boundary': {'h_run_coefficients': h_coeff,
                                   'q_run_coefficients': q_coeff,
                                   'raw_U9_run_coefficients': raw_source_coeff,
                                   'literal_controls': boundary_controls,
                                   'nontransitive_module_hostile': {
                                       'module_rows': triangle, 'H_module': 3,
                                       'ordered_nonempty_two_cover_count': 6,
                                       'raw_NR_profile': raw_state,
                                       'modified_profile': modified_state,
                                       'H_drift': 18, 'repaired_readiness_defect': 3}}}
    (HERE / 'h9_inflation_receipt.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'passed': True, 'literal_graphs': 128,
                      'neutral_positive_triples': [[1, 1, 1]],
                      'H47_ready_profile': [47, 16, 17, 48]}))


if __name__ == '__main__':
    main()
