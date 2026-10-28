#!/usr/bin/env python3
"""Exact controls for one-active-block linear query compression.

No external dependencies. Assertions are not used: python -O runs all checks.
Run: python3 check_block_compression.py > block-compression.json
"""
from functools import lru_cache, reduce
from itertools import combinations, product
import json
from operator import xor


def need(condition, label):
    if not condition:
        raise RuntimeError(label)


def span(columns):
    out = {0}
    for c in columns:
        out |= {x ^ c for x in out}
    return out


def image(columns, x):
    return reduce(xor, (c for i, c in enumerate(columns) if x >> i & 1), 0)


def active_states(r, k):
    return [a << (i * k) for i in range(r) for a in range(1, 1 << k)]


def block_weight(x, r, k):
    return sum(bool(x >> (i * k) & ((1 << k) - 1)) for i in range(r))


def partial_spread(blocks, k):
    seen = set()
    for cols in blocks:
        nonzero = span(cols) - {0}
        if len(nonzero) != (1 << k) - 1 or seen & nonzero:
            return False
        seen |= nonzero
    return True


def all_maps_control(r, k, s):
    """Every binary s by rk matrix, with no rank/symmetry filter."""
    states = active_states(r, k)
    total = good = 0
    for cols in product(range(1 << s), repeat=r * k):
        injective = len({image(cols, x) for x in states}) == len(states)
        spread = partial_spread([cols[i*k:(i+1)*k] for i in range(r)], k)
        kernel_ok = all(image(cols, x) != 0
                        for x in range(1, 1 << (r*k))
                        if block_weight(x, r, k) <= 2)
        need(injective == spread == kernel_ok, 'three equivalent criteria')
        total += 1
        good += injective
    return {'r': r, 'k': k, 's': s, 'matrices': total, 'separating': good}


def kernel_control():
    """All subspaces of F2^6 of dimensions 0,1,2,3; no orbit filter."""
    levels = [{(0,)}]
    for d in range(3):
        nxt = set()
        for space in levels[-1]:
            for x in range(1, 64):
                if x not in space:
                    nxt.add(tuple(sorted(set(space) | {a ^ x for a in space})))
        levels.append(nxt)
    good = [{space for space in level
             if all(block_weight(x, 3, 2) >= 3 for x in space if x)}
            for level in levels]
    need([len(z) for z in levels] == [1, 63, 651, 1395], 'subspace counts')
    need([len(z) for z in good] == [1, 27, 36, 0], 'avoiding kernel counts')
    gl2 = [(a, b) for a in range(1, 4) for b in range(1, 4) if a != b]
    diagonal = {tuple(sorted(x | (image(a, x) << 2) | (image(b, x) << 4)
                             for x in range(4)))
                for a in gl2 for b in gl2}
    need(diagonal == good[2], '36 kernels exactly twisted diagonals')
    return {'ambient_dimension': 6, 'block_dimensions': [2, 2, 2],
            'subspaces_by_dimension_0_to_3': [len(z) for z in levels],
            'valid_by_dimension_0_to_3': [len(z) for z in good],
            'all_valid_dimension_2_kernels': sorted(map(list, good[2]))}


def fixed_point_free_columns(m):
    """T and T+I invertible, from companion blocks of sizes 2 and 3."""
    need(m >= 2, 'T dimension')
    sizes = ([3] if m % 2 else []) + [2] * ((m - (3 if m % 2 else 0)) // 2)
    out = []
    offset = 0
    for size in sizes:
        out.extend(v << offset for v in ([2, 3] if size == 2 else [2, 4, 3]))
        offset += size
    need(offset == m, 'T full dimension')
    need(len(span(out)) == 1 << m, 'T invertible')
    need(len(span([v ^ (1 << i) for i, v in enumerate(out)])) == 1 << m,
         'T plus I invertible')
    return out


def line_spread(s):
    """Optimal binary partial line spread, returned as pairs of basis columns."""
    need(s >= 2, 'line ambient dimension')
    if s in (2, 3):
        return [(1, 2)]
    m = s - 2
    t = fixed_point_free_columns(m)
    graphs = [(a | (1 << m), image(t, a) | (1 << (m + 1)))
              for a in range(1 << m)]
    return graphs + line_spread(m)


def maximum_lines(s):
    return ((1 << s) - (1 if s % 2 == 0 else 5)) // 3


def line_construction_control():
    rows = []
    for s in range(2, 11):
        blocks = line_spread(s)
        need(partial_spread(blocks, 2), 'constructed lines pairwise disjoint')
        need(len(blocks) == maximum_lines(s), 'optimal line construction size')
        used = set().union(*(span(b) - {0} for b in blocks))
        holes = sorted(set(range(1, 1 << s)) - used)
        need(len(holes) == (0 if s % 2 == 0 else 4), 'hole count')
        need(reduce(xor, range(1, 1 << s), 0) == 0, 'ambient XOR')
        need(all(reduce(xor, span(b) - {0}, 0) == 0 for b in blocks), 'line XOR')
        need(reduce(xor, holes, 0) == 0, 'hole XOR')
        rows.append({'s': s, 'maximum_blocks': len(blocks), 'holes': holes})
    minima = []
    for r in range(1, 11):
        s = next(s for s in range(2, 7) if maximum_lines(s) >= r)
        blocks = line_spread(s)[:r]
        cols = sum((list(b) for b in blocks), [])
        states = active_states(r, 2)
        need(len({image(cols, x) for x in states}) == 3*r, 'all active states separated')
        minima.append({'r': r, 'states': 3*r, 'minimum_fixed_linear_queries': s,
                       'columns': cols})
    need([z['minimum_fixed_linear_queries'] for z in minima]
         == [2, 4, 4, 4, 4, 5, 5, 5, 5, 6], 'r 1 to 10 exact minima')
    return {'optimal_constructions': rows, 'r_1_to_10': minima,
            'nine_lines_in_F2_5': [sorted(span(b)-{0}) for b in line_spread(5)]}


def local_adaptive_control(r, k):
    """Exact minimax over candidate subsets; all block-local nonzero forms."""
    states = active_states(r, k)
    probes = active_states(r, k)
    ones = [sum((1 << i) for i, x in enumerate(states) if (x & w).bit_count() % 2)
            for w in probes]

    @lru_cache(None)
    def depth(mask):
        if mask.bit_count() <= 1:
            return 0
        best = len(states)
        for q in ones:
            a, b = mask & q, mask & ~q
            if a and b:
                best = min(best, 1 + max(depth(a), depth(b)))
        return best

    ans = depth((1 << len(states)) - 1)
    expected = k if r == 1 else r*k - 1
    need(ans == expected, 'adaptive block-local exact minimax')
    return {'r': r, 'k': k, 'states': len(states), 'minimum_worst_case_depth': ans,
            'candidate_subsets_evaluated': depth.cache_info().currsize}


def exceptional_controls():
    # k=1, nozero state universe: a zero column is legal, and r=1 takes zero bits.
    rows = []
    for r in range(1, 5):
        s = (r - 1).bit_length()
        cols = list(range(r))
        need(all(c < 1 << s for c in cols), 'k1 counting witness')
        need(len({image(cols, x) for x in active_states(r, 1)}) == r,
             'k1 permits a zero column')
        with_zero = r.bit_length()
        nzcols = list(range(1, r + 1))
        need(all(c < 1 << with_zero for c in nzcols), 'k1 zero-included witness')
        need(len({image(nzcols, x) for x in [0] + active_states(r, 1)}) == r+1,
             'k1 zero-included states')
        rows.append({'r': r, 'without_zero_bits': s, 'with_zero_bits': with_zero})
    # Distinguishing just kernel nonzeros of weight1 is not enough: two-block collision.
    cols = [1, 2, 1, 2]
    need(len(span(cols[:2])) == len(span(cols[2:])) == 4, 'restrictions injective')
    need(image(cols, 1) == image(cols, 4), 'hostile cross-block collision')
    # An adaptive local tree can beat every fixed linear map (r=2,k=2: 3 vs4).
    need(local_adaptive_control(2, 2)['minimum_worst_case_depth'] == 3,
         'fixed versus adaptive distinction')
    diagonal_map = [5, 10, 1, 2, 4, 8]  # (x,y,z) -> (x+y,x+z)
    for a in range(1, 4):
        need(image(diagonal_map, a | (a << 2)) == image(diagonal_map, a << 4),
             'one-block/two-block additive closure collision')
    c2 = [x for x in range(64) if block_weight(x, 3, 2) <= 2]
    need(len(c2) == 37 and len({image(diagonal_map, x) for x in c2}) == 16,
         'four bits cannot distinguish full two-block closure')
    for t in range(1, 4):
        ct = [x for x in range(64) if block_weight(x, 3, 2) <= t]
        differences = {x ^ y for x in ct for y in ct}
        expected = {x for x in range(64) if block_weight(x, 3, 2) <= min(2*t, 3)}
        need(differences == expected, 'sparse closure difference identity')
    gl4 = 1
    for i in range(4):
        gl4 *= 16 - (1 << i)
    need(gl4 == 20160 and 36 * gl4 == 725760 and 36 * gl4 // 24 == 30240,
         'ordered and unordered optimal probe counts')
    return {'k1': rows, 'hostile_two_identical_image_blocks': cols,
            'nozero_does_not_mean_no_zero_column_when_k_is_1': True,
            'closure_hostile': {'map_columns': diagonal_map,
                                'two_block_states_including_zero': len(c2),
                                'distinct_outputs': 16},
            'optimal_three_two_bit_block_probes': {
                'ordered_four_query_lists': 36 * gl4,
                'unordered_four_query_sets': 36 * gl4 // 24}}


def main():
    report = {
        'schema': 'one-active-block-compression-v1',
        'universe': 'V=(F2^k)^r, U consists of exactly one nonzero block',
        'all_maps_controls': [all_maps_control(r, k, s)
                              for r, k, s in [(1, 2, 2), (2, 2, 2), (2, 2, 3),
                                              (2, 2, 4), (3, 2, 2)]],
        'kernel_census': kernel_control(),
        'partial_line_spreads': line_construction_control(),
        'adaptive_local_minimax': [local_adaptive_control(r, k)
                                   for k, rs in [(2, range(1, 7)),
                                                 (3, range(1, 4)),
                                                 (4, range(1, 3))] for r in rs],
        'exceptional_and_hostile_controls': exceptional_controls(),
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
