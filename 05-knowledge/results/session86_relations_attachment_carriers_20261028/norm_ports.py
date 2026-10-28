#!/usr/bin/env python3
"""Exact normalized port and restricted partner-probe audit; stdlib only."""
import itertools
import json
from functools import lru_cache
from check_field_forms import Field, rank, nullspace, need


def audit(n, polynomial):
    f = Field(n, polynomial)
    m, q = n // 2, 1 << (n // 2)
    need(f.power(2, f.size - 1) == 1 and
         {f.power(2, i) for i in range(f.size - 1)} == set(range(1, f.size)),
         'alpha generates all nonzero residues, certifying primitive field model')
    U = [u for u in range(1, f.size) if f.power(u, q + 1) == 1]
    need(len(U) == q + 1 and rank(U) == n, 'norm-one size and full span')
    Kstar = [k for k in range(1, f.size) if f.power(k, q) == k]
    norms = {x: f.power(x, q + 1) for x in range(1, f.size)}
    normalized = {x: f.mul(x, f.power(f.power(norms[x], 1 << (m - 1)), f.size - 2))
                  for x in range(1, f.size)}
    need(all(normalized[x] in U for x in normalized), 'normalization has norm one')
    need(all(normalized[f.mul(x, k)] == normalized[x] for x in normalized for k in Kstar),
         'K-scaling invariance')
    need(all(normalized[f.mul(x, y)] == f.mul(normalized[x], normalized[y])
             for x in normalized for y in normalized), 'multiplicative normalization')
    need(all(normalized[u] == u for u in U), 'retraction onto U')
    need(all(f.mul(u, v) in U for u in U for v in U), 'closed under product')
    need(all(f.power(u, f.size - 2) in U and f.mul(u, u) in U for u in U),
         'closed under inverse and Frobenius')
    lines = sorted({tuple(sorted(f.mul(x, k) for k in Kstar)) for x in normalized})
    line_ports = [normalized[L[0]] for L in lines]
    need(sorted(line_ports) == U, 'one normalized port per line')
    rows = [[f.trace(f.mul(u, f.power(v, q))) for v in U] for u in U]
    need(all(rows[i][j] == f.trace(f.mul(u, f.power(v, f.size - 2)))
             for i, u in enumerate(U) for j, v in enumerate(U)), 'pairing is trace of ratio')
    need(len(set(map(tuple, rows))) == len(U), 'full pairing response separates U')
    tested_counts, successful = [], []
    for k in range(n + 1):
        successful = []
        tested = 0
        for inds in itertools.combinations(range(len(U)), k):
            tested += 1
            if len({tuple(row[i] for i in inds) for row in rows}) == len(U):
                successful.append(list(inds))
        tested_counts.append({'size': k, 'tested': tested, 'injective': len(successful)})
        if successful:
            break
    need(k == {4: 4, 6: 6, 8: 5}[n], 'minimum restricted partner probes')
    need(len(successful) == {4: 5, 6: 27, 8: 714}[n], 'all minimum probe-set counts')
    adaptive = {'minimum_worst_case_depth': {4: 4, 6: 5, 8: 5}[n]}
    if n in (4, 6):
        # Exact deterministic adaptive minimax recurrence on every candidate subset.
        one_masks = [sum(rows[i][j] << i for i in range(len(U))) for j in range(len(U))]
        decisions = {}

        @lru_cache(None)
        def depth(states):
            if states.bit_count() <= 1:
                return 0
            best = len(U) + 1
            for j, ones in enumerate(one_masks):
                yes, no = states & ones, states & ~ones
                if not yes or not no:
                    continue
                value = 1 + max(depth(yes), depth(no))
                if value < best:
                    best, decisions[states] = value, j
            need(best <= len(U), 'adaptive subsets separable')
            return best

        all_depths = [depth(states) for states in range(1 << len(U))]
        full = (1 << len(U)) - 1
        need(depth(full) == adaptive['minimum_worst_case_depth'], 'adaptive minimax cost')

        def tree(states):
            if states.bit_count() == 1:
                return {'port': U[states.bit_length() - 1]}
            j = decisions[states]
            return {'partner': U[j], 'zero': tree(states & ~one_masks[j]),
                    'one': tree(states & one_masks[j])}

        decision_tree = tree(full)
        path_lengths = []
        for i, u in enumerate(U):
            node, used = decision_tree, 0
            while 'partner' in node:
                response = rows[i][U.index(node['partner'])]
                node = node['one' if response else 'zero']
                used += 1
            need(node['port'] == u, 'adaptive witness identifies every port')
            path_lengths.append(used)
        need(max(path_lengths) == depth(full), 'adaptive witness attains optimum')
        adaptive.update({'method': 'complete candidate-subset minimax DP',
                         'subset_count': len(all_depths), 'all_subset_depths': all_depths,
                         'decision_tree': decision_tree, 'path_lengths_in_port_order': path_lengths})
        if n == 6:
            need(all((full & ~a).bit_count() == 7 for a in one_masks), 'first zero leaves seven')
            min_after_two_zeros = min((full & ~a & ~b).bit_count() for a in one_masks for b in one_masks)
            need(min_after_two_zeros == 5, 'two-zero adaptive lower bound')
            adaptive['minimum_states_after_any_two_zero_answers'] = min_after_two_zeros
    else:
        need((len(U) - 1).bit_length() == k == adaptive['minimum_worst_case_depth'],
             'adaptive counting bound meets fixed witness')
        adaptive['method'] = 'binary counting lower bound equals fixed five-probe witness'
    triple_blocks = []
    if n == 6:
        mu3 = [u for u in U if f.power(u, 3) == 1]
        triple_blocks = sorted({tuple(sorted(f.mul(u, z) for z in mu3)) for u in U})
        need(len(triple_blocks) == 3 and all(rank(T) == 2 for T in triple_blocks), 'three planes')
        block_of = {u: i for i, T in enumerate(triple_blocks) for u in T}
        need(all(rows[i][j] == int(u != v and block_of[u] == block_of[v])
                 for i, u in enumerate(U) for j, v in enumerate(U)), 'orthogonal triple blocks')
        need(all(all(sum(U[i] in T for i in inds) == 2 for T in triple_blocks)
                 for inds in successful), 'every minimum takes two probes per triple')
    pairing_graph = {}
    if n == 4:
        need(all(rows[i][j] == int(i != j) for i in range(len(U)) for j in range(len(U))),
             'F16 normalized pairing is K5')
        pairing_graph = {'type': 'complete_graph', 'vertices': 5}
    if n == 8:
        zeta = 36
        exponent_vertices = [f.power(zeta, i) for i in range(17)]
        squares = {i * i % 17 for i in range(1, 17)}
        need(set(exponent_vertices) == set(U), 'Paley exponent vertex marking')
        need(all(f.trace(f.mul(u, f.power(v, q))) == int((i - j) % 17 in squares)
                 for i, u in enumerate(exponent_vertices) for j, v in enumerate(exponent_vertices)),
             'normalized pairing graph is exactly Paley17')
        pairing_graph = {'type': 'Paley_graph', 'vertices': 17, 'generator': zeta,
                         'exponent_vertices': exponent_vertices,
                         'quadratic_residues_mod_17': sorted(squares)}
    unrestricted = {'minimum_linear_probes': (len(U) - 1).bit_length()}
    if n == 4:
        differences = {u ^ v for u, v in itertools.combinations(U, 2)}
        avoiding = [(a,) for a in range(1, f.size) if a not in differences]
        need(len(differences) == 10 and len(avoiding) == 5, 'F16 one-plane kernels')
        kernel = avoiding[0]
        constraints = [sum(f.trace(f.mul(a, f.power(1 << i, q))) << i for i in range(n))
                       for a in kernel]
        probes = nullspace(constraints, n)
        probe_rows = [[f.trace(f.mul(u, f.power(w, q))) for w in probes] for u in U]
        need(len(probes) == 3 and len(set(map(tuple, probe_rows))) == len(U),
             'F16 unrestricted three-probe witness')
        unrestricted.update({'difference_set_size': len(differences), 'one_planes_tested': f.size - 1,
                             'avoiding_kernel_count': len(avoiding), 'all_avoiding_kernels': avoiding,
                             'example_kernel_nonzero': kernel, 'example_probe_values': probes,
                             'example_response_rows': probe_rows})
    if n == 6:
        differences = {u ^ v for u, v in itertools.combinations(U, 2)}
        planes = {tuple(sorted((a, b, a ^ b)))
                  for a, b in itertools.combinations(range(1, f.size), 2)}
        avoiding = [P for P in sorted(planes) if not (set(P) & differences)]
        need(len(planes) == 651 and len(avoiding) == 36 and len(differences) == 36,
             'all binary two-plane kernels')
        kernel = avoiding[0]
        constraints = [sum(f.trace(f.mul(a, f.power(1 << i, q))) << i for i in range(n))
                       for a in kernel]
        probes = nullspace(constraints, n)
        probe_rows = [[f.trace(f.mul(u, f.power(w, q))) for w in probes] for u in U]
        need(len(probes) == 4 and len(set(map(tuple, probe_rows))) == len(U),
             'unrestricted four-probe witness')
        unrestricted.update({'difference_set_size': len(differences), 'two_planes_tested': len(planes),
                             'avoiding_kernel_count': len(avoiding), 'all_avoiding_kernels': avoiding,
                             'example_kernel_nonzero': kernel, 'example_probe_values': probes,
                             'example_response_rows': probe_rows})
    witness = successful[0]
    if n == 8:
        unrestricted.update({'example_probe_values': [U[i] for i in witness],
                             'example_response_rows': [[row[i] for i in witness] for row in rows]})
    return {'degree': n, 'q': q, 'polynomial_mask': polynomial,
            'norm_one_ports': U, 'pairing_graph': pairing_graph, 'span_dimension': rank(U),
            'additive_closure_size': 1 << rank(U),
            'unrestricted_linear_probes': unrestricted,
            'adaptive_partner_queries': adaptive,
            'line_ports_in_sorted_spread_order': line_ports,
            'pairing_rows': rows, 'information_lower_bound': (len(U) - 1).bit_length(),
            'probe_subset_search': tested_counts, 'minimum_partner_probes': k,
            'minimum_probe_set_count': len(successful),
            'all_minimum_probe_index_sets': successful,
            'example_probe_values': [U[i] for i in witness],
            'example_response_rows': [[row[i] for i in witness] for row in rows],
            'f64_orthogonal_triple_blocks': triple_blocks}


if __name__ == '__main__':
    print(json.dumps({'scope': 'Norm-one representatives, probes restricted to norm-one partners',
                      'cases': [audit(6, 0x43), audit(8, 0x11d)],
                      'controls': [audit(4, 0x13)]}, indent=2, sort_keys=True))
