#!/usr/bin/env python3
"""Exact carousel/block taxonomy and a labelled-C3 automorphism hostile.

Standard library only. The graph-isomorphism oracle sees adjacency matrices,
not modules, bit words, or necklaces. All gates survive Python -O.
"""
from collections import Counter
from itertools import combinations, permutations, product
from functools import lru_cache
from pathlib import Path
import argparse
import json

GATES = 0


def check(value, message):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(message)


def carousel(n):
    return tuple(sum(1 << j for j in range(n) if 1 <= (j-i) % n <= n//2) for i in range(n))


def opposite(rows):
    full = (1 << len(rows)) - 1
    return tuple(full ^ row ^ (1 << i) for i, row in enumerate(rows))


def inflate(base, blocks):
    offsets = [0]
    for block in blocks:
        offsets.append(offsets[-1] + len(block))
    output = []
    for i, block in enumerate(blocks):
        external = sum(((1 << len(blocks[j])) - 1) << offsets[j]
                       for j in range(len(base)) if base[i] >> j & 1)
        output.extend(external | (row << offsets[i]) for row in block)
    return tuple(output)


def cyclic_triples(rows):
    return frozenset(t for t in combinations(range(len(rows)), 3)
                     if all(sum(rows[i] >> j & 1 for j in t if i != j) == 1 for i in t))


def modules(rows):
    full = (1 << len(rows)) - 1
    return [s for s in range(1, full+1)
            if all(rows[v] & s in (0, s) for v in range(len(rows)) if not s >> v & 1)]


def iso_count(source, target, stop_after_one=False):
    """Full adjacency backtracking, pruned only by necessary graph conditions."""
    n = len(source)
    if n != len(target):
        return 0
    ds = [r.bit_count() for r in source]
    dt = [r.bit_count() for r in target]
    if sorted(ds) != sorted(dt):
        return 0
    classes = {d: tuple(i for i in range(n) if dt[i] == d) for d in set(dt)}
    mapping, used = {}, set()

    def recurse():
        if len(mapping) == n:
            return 1
        chosen, candidates = None, None
        for v in range(n):
            if v in mapping:
                continue
            cv = [w for w in classes[ds[v]] if w not in used
                  and all((source[v] >> u & 1) == (target[w] >> z & 1)
                          for u, z in mapping.items())]
            if not cv:
                return 0
            if candidates is None or len(cv) < len(candidates):
                chosen, candidates = v, cv
        total = 0
        for w in candidates:
            mapping[chosen] = w
            used.add(w)
            total += recurse()
            used.remove(w)
            del mapping[chosen]
            if total and stop_after_one:
                return total
        return total
    return recurse()


def hamilton_paths(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            count = dp[mask][last]
            if count:
                for nxt in range(n):
                    if not mask >> nxt & 1 and rows[last] >> nxt & 1:
                        dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def rotations(word):
    return {word[i:] + word[:i] for i in range(len(word))}


def period(word):
    return next(d for d in range(1, len(word)+1)
                if all(word[i] == word[(i+d) % len(word)] for i in range(len(word))))


def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]


def mobius(n):
    result, p = 1, 2
    while p*p <= n:
        if n % p == 0:
            n //= p
            result = -result
            if n % p == 0:
                return 0
        p += 1
    return -result if n > 1 else result


def primitive_words(d):
    return sum(mobius(e) * 2**(d//e) for e in divisors(d))


def induced(rows, vertices):
    return tuple(sum((rows[v] >> w & 1) << j for j, w in enumerate(vertices)) for v in vertices)


def order_bits(rows, ordering):
    return sum((rows[ordering[i]] >> ordering[j] & 1) << bit
               for bit, (i, j) in enumerate(combinations(range(len(rows)), 2)))


@lru_cache(None)
def tree_code(rows):
    """Brute-force strong modules and coloured quotient canonization."""
    n = len(rows)
    if n <= 1:
        return (n,)
    ms = modules(rows)
    strong = [s for s in ms if all(not s & t or s & t in (s, t) for t in ms)]
    proper = [s for s in strong if s != (1 << n)-1]
    children = [s for s in proper if not any(s != t and s & t == s for t in proper)]
    vertices = [tuple(i for i in range(n) if child >> i & 1) for child in children]
    codes = [tree_code(induced(rows, vs)) for vs in vertices]
    q = induced(rows, tuple(vs[0] for vs in vertices))
    degrees = [row.bit_count() for row in q]
    if sorted(degrees) == list(range(len(q))):
        ordering = sorted(range(len(q)), key=lambda i: -degrees[i])
        return (2, tuple(codes[i] for i in ordering))
    decorated = min((order_bits(q, p), tuple(codes[i] for i in p))
                    for p in permutations(range(len(q))))
    return (3, len(q), decorated)


def small_tree_oracle():
    result = []
    expected = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12}
    for n in range(1, 6):
        pairs = tuple(combinations(range(n), 2))
        orders = tuple(permutations(range(n)))
        tree_to_full, full_to_tree = {}, {}
        for mask in range(1 << len(pairs)):
            rows = [0] * n
            for bit, (i, j) in enumerate(pairs):
                if mask >> bit & 1:
                    rows[i] |= 1 << j
                else:
                    rows[j] |= 1 << i
            rows = tuple(rows)
            code = tree_code(rows)
            full = min(order_bits(rows, p) for p in orders)
            check(code not in tree_to_full or tree_to_full[code] == full, 'tree code implies full isomorphism')
            check(full not in full_to_tree or full_to_tree[full] == code, 'full isomorphism implies tree code')
            tree_to_full[code], full_to_tree[full] = full, code
        check(len(full_to_tree) == expected[n], 'small complete class census')
        result.append({'order': n, 'all_labelled_tournaments': 1 << len(pairs),
                       'tree_classes': len(tree_to_full), 'full_permutation_classes': len(full_to_tree)})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    # A=C3 -> singleton, with C3 arcs 0->1->2->0.
    a = (10, 12, 9, 0)
    a_op = opposite(a)
    check(iso_count(a, a) == iso_count(a_op, a_op) == 3, 'block automorphism orders')
    check(iso_count(a, a_op) == 0, 'opposite blocks are nonisomorphic')
    check(cyclic_triples(a) == cyclic_triples(a_op), 'opposite block cyclic triples')
    check(all(iso_count(carousel(n), carousel(n)) == n for n in (3, 5, 7, 9)), 'carousel full automorphisms')
    carousel_modules = {}
    for n in (3, 5, 7, 9, 11, 13):
        ms = modules(carousel(n))
        check(len(ms) == n+1, 'carousel is indecomposable in finite control')
        carousel_modules[n] = len(ms)

    t0 = inflate(carousel(3), [a, a, a])
    t1 = inflate(carousel(3), [a, a, a_op])
    c0, c1 = cyclic_triples(t0), cyclic_triples(t1)
    aut0, aut1 = iso_count(t0, t0), iso_count(t1, t1)
    check(c0 == c1, 'literal labelled C3 equality for 81/27 witness')
    check((aut0, aut1) == (81, 27), 'full graph oracle 81/27')
    check(iso_count(t0, t1) == 0, '81/27 witness nonisomorphism')
    h0, h1 = hamilton_paths(t0), hamilton_paths(t1)
    check(h0 == h1, 'Hamilton count survives the hostile')
    block_masks = {15, 15 << 4, 15 << 8}
    for graph in (t0, t1):
        proper = [m for m in modules(graph) if m != (1 << 12)-1]
        maximal = {m for m in proper if not any(m != other and m & other == m for other in proper)}
        check(maximal == block_masks, 'full subset module oracle recovers blocks')

    graph_oracle = []
    for n in (3, 5):
        words = list(product((0, 1), repeat=n))
        graphs = {w: inflate(carousel(n), [a_op if bit else a for bit in w]) for w in words}
        representatives = {w: min(rotations(w)) for w in words}
        auts = {}
        for w in words:
            auts[w] = iso_count(graphs[w], graphs[w])
            check(auts[w] == 3**n * (n//period(w)), 'full adjacency automorphism formula')
            check(iso_count(graphs[w], opposite(graphs[w])) == 0, 'no self-converse member')
        pair_checks = 0
        for i, w in enumerate(words):
            for v in words[i:]:
                isomorphic = iso_count(graphs[w], graphs[v], stop_after_one=True) != 0
                check(isomorphic == (representatives[w] == representatives[v]), 'full graph necklace iff')
                pair_checks += 1
        graph_oracle.append({'base_order': n, 'expanded_order': 4*n,
                             'words': len(words), 'all_unordered_pairs_with_diagonal': pair_checks,
                             'necklace_classes': len(set(representatives.values())),
                             'automorphism_histogram_on_words': dict(sorted(Counter(auts.values()).items()))})

    period_census = []
    for n in (3, 5, 7, 9, 11, 13, 15):
        reps = {min(rotations(w)) for w in product((0, 1), repeat=n)}
        strata = Counter(period(w) for w in reps)
        predicted = {d: primitive_words(d)//d for d in divisors(n)}
        check(dict(strata) == predicted, 'Mobius least-period strata')
        conjugate = lambda w: min(rotations(tuple(1-bit for bit in reversed(w))))
        check(all(conjugate(w) != w and conjugate(conjugate(w)) == w for w in reps),
              'complemented-reflection converse pairing')
        period_census.append({'n': n, 'classes': len(reps), 'least_period_counts': dict(sorted(strata.items())),
                              'converse_merged_classes': len(reps)//2})

    # Both a missing block-type hypothesis and a larger circulant group are tested.
    b = (2, 0)
    selfdual_graphs = [inflate(carousel(3), [opposite(b) if bit else b for bit in w])
                      for w in product((0, 1), repeat=3)]
    check(all(iso_count(selfdual_graphs[0], g, stop_after_one=True) for g in selfdual_graphs),
          'self-converse block destroys binary faithfulness')
    paley = tuple(sum(1 << j for j in range(7) if (j-i) % 7 in (1, 2, 4)) for i in range(7))
    check(iso_count(paley, paley) == 21, 'Paley7 multiplier hostile')
    check(all((paley[i] >> j & 1) == (paley[(2*i) % 7] >> ((2*j) % 7) & 1)
              for i in range(7) for j in range(7)), 'Paley multiplier by two is an automorphism')
    first = (1, 1, 0, 0, 0, 0, 0)
    second = tuple(first[(4*i) % 7] for i in range(7))  # image under multiplication by two
    check(second not in rotations(first), 'multiplier identifies distinct necklaces')
    check(iso_count(inflate(paley, [a_op if bit else a for bit in first]),
                    inflate(paley, [a_op if bit else a for bit in second]), stop_after_one=True) != 0,
          'full graph oracle confirms Paley non-necklace identification')
    linear3, singleton = (6, 4, 0), (0,)
    cut_a = inflate(linear3, [singleton, b, singleton])
    cut_b = inflate(linear3, [b, singleton, singleton])
    check(iso_count(cut_a, cut_b) == 1, 'nonintrinsic block boundaries can disappear')
    check(tree_code(cut_a) == tree_code(cut_b), 'canonical tree removes arbitrary cuts')
    strong6 = inflate(carousel(3), [a, singleton, singleton])
    check(cyclic_triples(strong6) == cyclic_triples(opposite(strong6))
          and iso_count(strong6, opposite(strong6)) == 0,
          'strong-component isomorphism types are not C3-invariant')
    tree_oracle = small_tree_oracle()
    result = {'status': 'FINITE-EXACT controls for all-order proved taxonomy', 'gates': GATES,
              'block_A_rows': a, 'block_Aop_rows': a_op,
              'counterexample': {'order': 12, 'uniform_rows': t0, 'mixed_rows': t1,
                                'identical_labelled_cyclic_triples': sorted(c0), 'c3_count': len(c0),
                                'automorphism_orders': [aut0, aut1], 'hamilton_paths': [h0, h1],
                                'maximal_proper_module_masks': sorted(block_masks)},
              'carousel_module_counts': carousel_modules, 'full_graph_oracle': graph_oracle,
              'period_census': period_census, 'full_strong_module_tree_oracle': tree_oracle,
              'hostiles': {'self_converse_block_all_eight_words_isomorphic': True,
                           'Paley7_automorphisms': 21, 'Paley_nonrotation_words': [first, second],
                           'arbitrary_linear_block_sizes': [[1, 2, 1], [2, 1, 1]],
                           'nonselfconverse_strong6_rows': strong6}}
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(payload)
    print(payload, end='')


if __name__ == '__main__':
    main()
