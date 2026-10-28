#!/usr/bin/env python3
"""Independent stdlib-only exact K4 / F64 connection checks. No removable asserts."""
import itertools as it
import json
from collections import Counter
from functools import lru_cache


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def mul(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 64:
            a ^= 0x43
    return value


def power(a, n):
    result = 1
    while n:
        if n & 1:
            result = mul(result, a)
        a = mul(a, a)
        n >>= 1
    return result


def trace(a):
    result = 0
    for _ in range(6):
        result ^= a
        a = mul(a, a)
    need(result in (0, 1), 'absolute trace is binary')
    return result


def pairing(a, b):
    return trace(mul(a, power(b, 8)))


def span(generators):
    result = {0}
    for g in generators:
        result |= {x ^ g for x in result}
    return result


EDGES = list(it.combinations(range(4), 2))


def boundary(mask):
    value = 0
    for i, (u, v) in enumerate(EDGES):
        if mask >> i & 1:
            value ^= (1 << u) | (1 << v)
    return value


def output(mask):
    return (boundary(mask) & 7) | ((mask.bit_count() % 2) << 3)


def partitions(values):
    if not values:
        yield ()
    else:
        a = values[0]
        for b in values[1:]:
            for rest in partitions([x for x in values[1:] if x != b]):
                yield ((a, b),) + rest


def port_masks(parts):
    return [
        x for a, b in parts for x in (1 << a, 1 << b, (1 << a) | (1 << b))]


def edge_permute(mask, p):
    return sum(1 << EDGES.index(tuple(sorted((p[u], p[v]))))
               for i, (u, v) in enumerate(EDGES) if mask >> i & 1)


def orientation_permute(mask, p):
    result = 0
    for i, (u, v) in enumerate(EDGES):
        source, target = (v, u) if mask >> i & 1 else (u, v)
        a, b = p[source], p[target]
        if a > b:
            result |= 1 << EDGES.index(tuple(sorted((a, b))))
    return result


def scores(mask):
    out = [0] * 4
    for i, (u, v) in enumerate(EDGES):
        out[v if mask >> i & 1 else u] += 1
    return tuple(sorted(out))


def outdegree_parities(mask):
    value = 0
    for i, (u, v) in enumerate(EDGES):
        value ^= 1 << (v if mask >> i & 1 else u)
    return value


def arc(mask, u, v):
    return bool(mask >> EDGES.index(tuple(sorted((u, v)))) & 1) == (u > v)


def main():
    need({power(2, i) for i in range(63)} == set(range(1, 64)), 'primitive field')
    U = [x for x in range(1, 64) if power(x, 9) == 1]
    block_bases = [(1, 58), (6, 26), (11, 20)]
    blocks = [sorted(span(b) - {0}) for b in block_bases]
    need(sorted(sum(blocks, [])) == U, 'three norm-one planes')
    need(len(span(sum((list(b) for b in block_bases), []))) == 64, 'direct sum')
    need(all(pairing(x, y) == int(i == j and x != y)
             for i, A in enumerate(blocks) for j, B in enumerate(blocks)
             for x in A for y in B), 'pairing graph is three triangles')
    D = {x ^ y for x, y in it.combinations(U, 2)}
    planes = {tuple(sorted((a, b, a ^ b))) for a, b in it.combinations(range(1, 64), 2)}
    avoiding = sorted(P for P in planes if not D.intersection(P))
    need(len(D) == 36 and len(planes) == 651 and len(avoiding) == 36, 'kernel census')
    fixed_counts = []
    for k in range(7):
        fixed_counts.append(sum(len({tuple(pairing(u, w) for w in W) for u in U}) == 9
                                for W in it.combinations(U, k)))
    need(fixed_counts == [0, 0, 0, 0, 0, 0, 27], 'restricted fixed costs')
    ones = [sum(pairing(u, w) << i for i, u in enumerate(U)) for w in U]

    @lru_cache(None)
    def depth(S):
        if S.bit_count() <= 1:
            return 0
        return 1 + min(max(depth(S & Q), depth(S & ~Q))
                       for Q in ones if S & Q and S & ~Q)

    depths = [depth(S) for S in range(512)]
    need(depths[-1] == 5, 'adaptive restricted cost')
    cuts = {sum(1 << i for i, (u, v) in enumerate(EDGES)
                if ((s >> u) ^ (s >> v)) & 1) for s in range(16)}
    cycles = {x for x in range(64) if boundary(x) == 0}
    bicycles = sorted(cuts & cycles)
    need(len(cuts) == len(cycles) == 8 and len(bicycles) == 4, 'cycle/cut dimensions')
    need(bicycles == [x for x in range(64) if output(x) == 0], 'bicycle kernel')
    need({x.bit_count() for x in bicycles[1:]} == {4}, 'bicycle four-cycles')
    parts = list(partitions(list(range(6))))
    successful_parts = [P for P in parts if len({output(x) for x in port_masks(P)}) == 9]
    all_adjacent = [P for P in parts if all(set(EDGES[a]) & set(EDGES[b]) for a, b in P)]
    need(len(parts) == 15 and len(successful_parts) == 8 and successful_parts == all_adjacent,
         'exact graph-scaffold criterion')
    scaffold_set = {frozenset(frozenset(pair) for pair in P) for P in successful_parts}
    rooted_scaffolds = set()
    for root in range(4):
        others = [v for v in range(4) if v != root]
        for cycle in (others, [others[0], others[2], others[1]]):
            P = []
            for i, v in enumerate(cycle):
                P.append(frozenset((EDGES.index(tuple(sorted((v, root)))),
                                    EDGES.index(tuple(sorted((v, cycle[(i + 1) % 3])))))))
            rooted_scaffolds.add(frozenset(P))
    need(rooted_scaffolds == scaffold_set, 'root plus directed triangle parametrization')
    opposite = [(0, 5), (1, 4), (2, 3)]
    adjacent_edges = [((0, 3), (0, 1)), ((1, 3), (1, 2)), ((2, 3), (0, 2))]
    adjacent = [(EDGES.index(a), EDGES.index(b)) for a, b in adjacent_edges]
    field_basis = sum((list(b) for b in block_bases), [])

    def field_to_edges(parts, bases=block_bases):
        table = {0: 0}
        for (a, b), (e, f) in zip(bases, parts):
            for v, bit in ((a, e), (b, f)):
                table.update({x ^ v: y ^ (1 << bit) for x, y in list(table.items())})
        need(len(table) == len(set(table.values())) == 64, 'linear coordinate isomorphism')
        return table

    phi = field_to_edges(adjacent)
    opp_phi = field_to_edges(opposite)
    need(all(pairing(x, y) == sum(((phi[x] >> a & 1) * (phi[y] >> b & 1)
                                  + (phi[x] >> b & 1) * (phi[y] >> a & 1))
                                 for a, b in adjacent) % 2
             for x in range(64) for y in range(64)), 'transported symplectic pairing')
    kernel = sorted(x for x in range(64) if output(phi[x]) == 0)
    probes = []
    for j in range(4):
        possible = [w for w in range(64) if all(pairing(x, w) == ((output(phi[x]) >> j) & 1)
                                                for x in range(64))]
        need(len(possible) == 1, 'nondegenerate trace probe representation')
        probes.append(possible[0])
    rows = [{'field': u, 'edge_mask': phi[u], 'response': output(phi[u]),
             'bits_b0_b1_b2_p': [(output(phi[u]) >> j) & 1 for j in range(4)],
             'tournament_scores': scores(phi[u])} for u in U]
    need(len({r['response'] for r in rows}) == 9, 'adjacent realization separates')
    need(output(phi[9]) == 0 and output(phi[mul(2, 9)]) != 0, 'multiplication does not descend')
    hostile = next((a, b, output(opp_phi[a])) for a, b in it.combinations(U, 2)
                   if output(opp_phi[a]) == output(opp_phi[b]))
    recovered_kernels = set()
    for bases in it.product(*(list(it.permutations(B, 2)) for B in blocks)):
        table = field_to_edges(adjacent, bases)
        recovered_kernels.add(tuple(sorted(x for x in range(1, 64) if output(table[x]) == 0)))
    need(recovered_kernels == set(avoiding), 'all 36 kernels realized by block bases')
    census = Counter(scores(x) for x in range(64))
    need(sorted(census.values()) == [8, 8, 24, 24], 'four labelled orbit sizes')
    permutations = list(it.permutations(range(4)))
    hamilton_counts = {x: sum(all(arc(x, p[i], p[i + 1]) for i in range(3)) for p in permutations)
                       for x in range(64)}
    expected_H = {(0, 1, 2, 3): 1, (0, 2, 2, 2): 3, (1, 1, 1, 3): 3, (1, 1, 2, 2): 5}
    need(all(hamilton_counts[x] == expected_H[scores(x)] for x in range(64)), 'actual Hamilton counts')
    permutation_orders = set()
    for p in permutations:
        q = tuple(range(4))
        for order in range(1, 25):
            q = tuple(p[q[i]] for i in range(4))
            if q == tuple(range(4)):
                permutation_orders.add(order)
                break
    need(permutation_orders == {1, 2, 3, 4}, 'no faithful cyclic order-nine action through S4')
    representatives = {min(orientation_permute(x, p) for p in permutations) for x in range(64)}
    need(len(representatives) == 4, 'four unlabelled tournament classes')
    for p in permutations:
        cocycle = orientation_permute(0, p)
        permutation_parity = sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4)) % 2
        need(cocycle.bit_count() % 2 == permutation_parity, 'orientation parity sign cocycle')
        need(all(orientation_permute(x, p) == (edge_permute(x, p) ^ cocycle)
                 for x in range(64)), 'affine orientation action')
    need(all(outdegree_parities(x) == (5 ^ boundary(x)) for x in range(64)), 'degree parity bridge')
    fibres = []
    for y in range(16):
        masks = [x for x in range(64) if output(x) == y]
        hist = Counter(scores(x) for x in masks)
        parities = {outdegree_parities(x) for x in masks}
        need(len(masks) == 4 and len(parities) == 1, 'four-element fibres and degree parity')
        p = parities.pop()
        expected = ({(0, 1, 2, 3): 2, (1, 1, 2, 2): 2} if p.bit_count() == 2 else
                    {(0, 2, 2, 2): 4} if p == 0 else {(1, 1, 1, 3): 4})
        need(hist == expected, 'entire decoder fibre classification')
        fibres.append({'response': y, 'masks': masks, 'outdegree_parity_mask': p,
                       'score_counts': {str(k): v for k, v in sorted(hist.items())}})
    output_not_class = next((a, b) for a, b in it.combinations(range(64), 2)
                           if output(a) == output(b) and scores(a) != scores(b))
    fano_columns = [boundary(1 << i) & 7 for i in range(6)]
    need(set(fano_columns) == set(range(1, 7)), 'Fano deletion six columns')
    need(all(sum(1 for c in fano_columns if c == a) == 1 for a in range(1, 7)), 'distinct columns')
    columns7 = fano_columns + [7]
    hamming = set()
    for word in range(128):
        syndrome = 0
        for i, col in enumerate(columns7):
            if word >> i & 1:
                syndrome ^= col
        if syndrome == 0:
            hamming.add(word)
    simplex = {sum(((a & col).bit_count() % 2) << i for i, col in enumerate(columns7))
               for a in range(8)}
    need({w & 63 for w in hamming if w < 64} == cycles, 'shortened Hamming is cycle code')
    need({w & 63 for w in simplex} == cuts, 'punctured simplex is cut code')
    need({w & 63 for w in simplex if w < 64} == set(bicycles), 'shortened simplex is bicycle')
    row_space = {w for w in range(64) if all((w & b).bit_count() % 2 == 0 for b in bicycles)}
    enumerators = {name: dict(sorted(Counter(w.bit_count() for w in code).items()))
                   for name, code in [('cycles', cycles), ('cuts', cuts),
                                      ('bicycles', bicycles), ('observations', row_space)]}
    need(enumerators['cycles'] == enumerators['cuts'] == {0: 1, 3: 4, 4: 3}, 'cycle/cut enumerators')
    need(enumerators['bicycles'] == {0: 1, 4: 3}, 'bicycle enumerator')
    need(enumerators['observations'] == {0: 1, 2: 3, 3: 8, 4: 3, 6: 1}, 'observation enumerator')
    print(json.dumps({
        'status': 'ALL CHECKS PASSED', 'source_pin': 'b8672af00bcc6d440a9a00a6915325faf39fcdff',
        'field_polynomial_mask': 67, 'norm_one_ports': U, 'block_bases': block_bases,
        'restricted_fixed_success_counts_sizes_0_to_6': fixed_counts,
        'adaptive_depth': depths[-1], 'difference_count': len(D), 'all_plane_count': len(planes),
        'avoiding_kernel_count': len(avoiding), 'edge_order': EDGES,
        'cycle_masks': sorted(cycles), 'cut_masks': sorted(cuts), 'bicycle_masks': bicycles,
        'all_pair_partitions': len(parts), 'successful_adjacent_pair_partitions': len(successful_parts),
        'rooted_cycle_scaffold_count': len(rooted_scaffolds), 'adjacent_scaffold_edges': adjacent_edges,
        'field_to_edge_basis_pairs': list(zip(field_basis, [e for P in adjacent for e in P])),
        'adjacent_field_kernel': kernel, 'adjacent_field_probes': probes,
        'adjacent_port_rows': rows, 'opposite_pair_hostile_ports_and_common_response': hostile,
        'opposite_pair_hostile_edge_masks': [opp_phi[hostile[0]], opp_phi[hostile[1]]],
        'all_36_kernels_recovered_by_basis_changes': len(recovered_kernels),
        'tournament_class_counts': {str(k): v for k, v in sorted(census.items())},
        'tournament_H_by_scores': {str(k): v for k, v in expected_H.items()},
        'S4_element_orders': sorted(permutation_orders),
        'four_bit_output_not_tournament_class_witness': output_not_class,
        'four_bit_fibres_all64': fibres,
        'quotient_multiplication_hostile': {'kernel_element': 9, 'multiplier': 2,
                                           'product': mul(2, 9), 'product_response': output(phi[mul(2, 9)])},
        'fano_deletion_columns': fano_columns,
        'binary_code_weight_enumerators': enumerators,
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
