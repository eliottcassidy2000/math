#!/usr/bin/env python3
"""Exact named-port controls. Standard library; checks survive python -O."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def rows(mask, n):
    result = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if mask >> k & 1:
            result[j] |= 1 << i
        else:
            result[i] |= 1 << j
    return tuple(result)


def add(parent, mask):
    # bit i=1 means old i -> new vertex.
    n = len(parent)
    return tuple(row | (((mask >> i) & 1) << n) for i, row in enumerate(parent)) + (((1 << n) - 1) ^ mask,)


def transport(parent, order):
    return tuple(sum(((parent[v] >> w) & 1) << j for j, w in enumerate(order)) for v in order)


def hamilton(parent):
    n = len(parent)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            nxt = parent[last] & ~mask
            while nxt:
                bit = nxt & -nxt
                dp[mask | bit][bit.bit_length() - 1] += value
                nxt -= bit
    return sum(dp[-1])


def literal_h(parent):
    n = len(parent)
    return sum(all(parent[p[i]] >> p[i + 1] & 1 for i in range(n - 1)) for p in permutations(range(n)))


def endpoint_cover_state(parent):
    n = len(parent)
    s, e = [0] * n, [0] * n
    Q = [[0] * n for _ in range(n)]
    for p in permutations(range(n)):
        valid = [bool(parent[p[i]] >> p[i + 1] & 1) for i in range(n - 1)]
        if all(valid):
            s[p[0]] += 1
            e[p[-1]] += 1
        for k in range(n - 1):
            if all(valid[:k]) and all(valid[k + 1:]):
                Q[p[k]][p[k + 1]] += 1
    return s, e, Q


def canonical(parent):
    return min(transport(parent, p) for p in permutations(range(len(parent))))


def determinant(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for k in range(len(a)):
        j = next((i for i in range(k, len(a)) if a[i][k]), None)
        if j is None:
            return Fraction(0)
        if j != k:
            a[k], a[j] = a[j], a[k]
            result = -result
        pivot = a[k][k]
        result *= pivot
        for i in range(k + 1, len(a)):
            factor = a[i][k] / pivot
            a[i] = [x - factor * y for x, y in zip(a[i], a[k])]
    return result


def rank(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    k = 0
    for col in range(len(a[0])):
        j = next((i for i in range(k, len(a)) if a[i][col]), None)
        if j is None:
            continue
        a[k], a[j] = a[j], a[k]
        pivot = a[k][col]
        a[k] = [x / pivot for x in a[k]]
        for i in range(k + 1, len(a)):
            factor = a[i][col]
            a[i] = [x - factor * y for x, y in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def run():
    parent_count = cut_count = 0
    for n in range(1, 6):
        for code in range(1 << (n * (n - 1) // 2)):
            parent = rows(code, n)
            H = hamilton(parent)
            s, e, Q = endpoint_cover_state(parent)
            require(sum(s) == sum(e) == H, 'literal endpoint/DP disagreement')
            w = [[(Q[i][j] + Q[j][i]) // 2 for j in range(n)] for i in range(n)]
            require(all((Q[i][j] + Q[j][i]) % 2 == 0 and w[i][j] > 0
                        for i in range(n) for j in range(i + 1, n)), 'integer positive pair weights')
            for i in range(n):
                require(s[i] - e[i] == sum(Q[i][j] for j in range(n) if parent[i] >> j & 1)
                        - sum(Q[j][i] for j in range(n) if parent[j] >> i & 1), 'endpoint divergence')
            for cut in range(1 << n):
                bits = [(cut >> i) & 1 for i in range(n)]
                endpoint = H + sum((e[i] - s[i]) * bits[i] for i in range(n))
                endpoint += sum(Q[i][j] * bits[i] * (1 - bits[j]) for i in range(n) for j in range(n))
                network = H + 2 * sum(w[i][j] * (1 - bits[i]) * bits[j]
                                      for i in range(n) for j in range(n) if parent[i] >> j & 1)
                require(endpoint == network == hamilton(add(parent, cut)), 'all-cut formula')
                cut_count += 1
            parent_count += 1

    positive, negative = (2, 4, 1), (4, 1, 2)
    first = [hamilton(add(positive, a)) for a in range(8)]
    require(first == [hamilton(add(negative, a)) for a in range(8)] == [3, 5, 5, 5, 5, 5, 5, 3], 'one-port equality')
    difference_counts, classification = Counter(), Counter()
    for a, b, c in product(range(8), range(8), range(2)):
        aa, bb = [(a >> i) & 1 for i in range(3)], [(b >> i) & 1 for i in range(3)]
        delta = (aa[1] - aa[0]) * (bb[2] - bb[0]) - (aa[2] - aa[0]) * (bb[1] - bb[0])
        plus, minus = add(add(positive, a), b | (c << 3)), add(add(negative, a), b | (c << 3))
        hp, hm = hamilton(plus), hamilton(minus)
        require(hp == literal_h(plus) and hm == literal_h(minus), 'two-port literal/DP control')
        require(hp - hm == 2 * delta * (1 - 2 * c + sum(bb) - sum(aa)), 'chirality determinant formula')
        repeated = len(set(zip(aa, bb))) < 3
        same_class = canonical(plus) == canonical(minus)
        require(same_class == repeated, 'exact two-port isomorphism criterion')
        difference_counts[hp - hm] += 1
        classification[(hp == hm, same_class)] += 1

    canonical_cases = []
    for missing in product(range(2), repeat=2):
        signatures = sorted(set(product(range(2), repeat=2)) - {missing})
        a = sum(u << i for i, (u, v) in enumerate(signatures))
        b = sum(v << i for i, (u, v) in enumerate(signatures))
        for c in range(2):
            children = [add(add(t, a), b | (c << 3)) for t in (positive, negative)]
            canonical_cases.append({'missing_signature': missing, 'new_arc_x_to_y': c,
                'a': a, 'b': b, 'H': [hamilton(t) for t in children],
                'deletion_H_decks': [sorted(hamilton(transport(t, tuple(j for j in range(5) if j != i)))
                                           for i in range(5)) for t in children]})

    sources = [rows(code, 3) for code in range(8)]
    one = [[hamilton(add(t, a)) for t in sources] for a in range(8)]
    contexts = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 9), (0, 10), (1, 2)]
    bank = [[hamilton(add(add(t, a), b)) for t in sources] for a, b in contexts]
    require(rank(one) == 4 and rank(bank) == 7, 'one/two-birth response rank')
    minor = determinant([row[:7] for row in bank])
    require(minor == -128, 'exact seven-query minor')
    require(all(sum(v * (-1 if code.bit_count() % 2 else 1) for code, v in enumerate(row)) == 0
                for row in bank), 'triangle Fourier kernel')

    module_controls = 0
    for quotient_order in range(1, 4):
        for code in range(1 << (quotient_order * (quotient_order - 1) // 2)):
            quotient = rows(code, quotient_order)
            n = quotient_order + 2
            def expand(cycle):
                r = list(cycle) + [0] * (quotient_order - 1)
                for old in range(1, quotient_order):
                    v = old + 2
                    for u in range(3):
                        if quotient[0] >> old & 1:
                            r[u] |= 1 << v
                        else:
                            r[v] |= 1 << u
                    for other in range(1, quotient_order):
                        if quotient[old] >> other & 1:
                            r[v] |= 1 << (other + 2)
                return tuple(r)
            p, m = expand(positive), expand(negative)
            for cut in range(1 << n):
                u, v = next((u, v) for u, v in combinations(range(3), 2) if ((cut >> u) & 1) == ((cut >> v) & 1))
                order = list(range(n + 1)); order[u], order[v] = order[v], order[u]
                require(transport(add(p, cut), order) == add(m, cut), 'module physical-port transport')
                module_controls += 1

    witness_plus, witness_minus = add(add(positive, 1), 2), add(add(negative, 1), 2)
    require((hamilton(witness_plus), hamilton(witness_minus)) == (13, 11), 'prime-target hostile')
    return {'status': 'FINITE-EXACT', 'universe': 'all labelled parents orders1..5; all128 two-new attachments of each named cyclic triangle',
        'all_cut_controls': {'parents': parent_count, 'cuts': cut_count},
        'one_new_triangle_responses': first,
        'difference_histogram': dict(sorted(difference_counts.items())),
        'two_new_classification': [{'same_H': h, 'same_class': c, 'count': count} for (h, c), count in sorted(classification.items())],
        'canonical_eight_cases': canonical_cases,
        'prime_target_witness': {'first_old_to_new_mask': 1, 'second_old_to_new_mask': 2,
            'positive_child_rows': witness_plus, 'negative_child_rows': witness_minus, 'H': [13, 11]},
        'linear_response_span': {'one_new_dimension': 4, 'two_new_dimension': 7, 'all_fixed_dimension': '7 (inherited; elementary no-old-triangle upper bound)',
            'seven_contexts': contexts, 'evaluation_matrix': bank, 'first_seven_columns_determinant': int(minor)},
        'triangle_module_controls': module_controls,
        'check_semantics': 'explicit exceptions active under python -O'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
