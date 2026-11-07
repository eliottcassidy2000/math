#!/usr/bin/env python3
"""Exact checks for the binary-comb path-count barrier; standard library only."""
from functools import lru_cache
from itertools import permutations, product
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def tournament(m, q, tau, plus):
    check(m >= 1 and 0 <= q < 2**m and tau in (0, 1), 'invalid parameters')
    n = m + 4
    a, b = n - 2, n - 1
    out = [0] * n
    for i in range(m + 2):
        for j in range(i + 1, m + 2):
            out[i] |= 1 << j
        u, v = (a, i) if i == 1 else (i, a)
        out[u] |= 1 << v
        toward_b = (i == 1 and tau) or (i >= 2 and q >> (i - 2) & 1)
        u, v = (i, b) if toward_b else (b, i)
        out[u] |= 1 << v
    u, v = (a, b) if plus else (b, a)
    out[u] |= 1 << v
    for i in range(n):
        check(not (out[i] >> i & 1), 'loop')
        for j in range(i + 1, n):
            check(((out[i] >> j & 1) + (out[j] >> i & 1)) == 1, 'not tournament')
    return out


def path_dp(out, start=None):
    n = len(out)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        if start is None or v == start:
            dp[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            choices = out[v] & (full ^ mask)
            while choices:
                bit = choices & -choices
                choices -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return dp[full]


def hamilton_count(out):
    return sum(path_dp(out))


def cycle_count(out):
    # Fix vertex zero as the cycle start; no rotational overcount.
    return sum(c for v, c in enumerate(path_dp(out, start=0)) if out[v] & 1)


def permutation_count(out):
    return sum(all(out[u] >> v & 1 for u, v in zip(word, word[1:]))
               for word in permutations(range(len(out))))


def odd_cycle_formula(out):
    # Independently list each odd directed cycle, anchored at its least vertex.
    # Since n <= 8 in these checks, complete permutations are intentionally used.
    n = len(out)
    cycles = [0] * (1 << n)
    for mask in range(1 << n):
        if mask.bit_count() < 3 or mask.bit_count() % 2 == 0:
            continue
        vertices = [v for v in range(n) if mask >> v & 1]
        first = vertices[0]
        for tail in permutations(vertices[1:]):
            word = (first,) + tail
            if all(out[u] >> v & 1 for u, v in zip(word, word[1:] + word[:1])):
                cycles[mask] += 1

    @lru_cache(None)
    def collections(mask):
        if not mask:
            return 1
        first = mask & -mask
        value = collections(mask ^ first)  # First vertex belongs to no cycle.
        sub = mask
        while sub:
            if sub & first and cycles[sub]:
                value += 2 * cycles[sub] * collections(mask ^ sub)
            sub = (sub - 1) & mask
        return value

    return collections((1 << n) - 1)


def pair_weight(j, t):
    return 2**t * (3**j + 1) - (1 + (-1)**j) * int(t == 0)


def certificates(j, t):
    for prefix in product(range(3), repeat=j):
        # 0 unused, 1 A, 2 B
        for suffix in product(range(2), repeat=t):
            A = [i + 2 for i, color in enumerate(prefix) if color == 1]
            A += [j + 3 + i for i, color in enumerate(suffix) if color == 1]
            B = [i + 2 for i, color in enumerate(prefix) if color == 2]
            if len(A) % 2 == 1 and len(B) % 2 == 0:
                yield A, B


def valid_cycle(out, word):
    return len(word) >= 3 and len(word) % 2 == 1 and len(set(word)) == len(word) and all(
        out[u] >> v & 1 for u, v in zip(word, word[1:] + word[:1]))


def check_growth(m, q, tau, plus, epsilon, old):
    new = tournament(m + 1, 2*q + epsilon, tau, plus)
    # New bit vertex has label 2; every old label >=2 shifts up by one.
    image = [v if v < 2 else v + 1 for v in range(len(old))]
    for u in range(len(old)):
        for v in range(len(old)):
            check((old[u] >> v & 1) == (new[image[u]] >> image[v] & 1), 'growth changes old arc')


def main():
    cases = perm_cases = ocf_cases = growth_cases = nonzero = 0
    order_data = []
    for m in range(1, 9):
        values = set()
        minimum_slack = None
        for q in range(2**m):
            for tau, plus in product((0, 1), repeat=2):
                out = tournament(m, q, tau, plus)
                h = hamilton_count(out)
                c = cycle_count(out)
                check(c == q + tau + plus, 'inherited cycle formula')
                check(h >= 1 + 2**m, 'baseline bound')
                if q:
                    j = q.bit_length() - 1
                    t = m - j - 1
                    bound = 1 + 2**m + pair_weight(j, t)
                    check(h >= bound, 'pair lower bound')
                    check(bound >= 2**t * 3**j, 'coarse lower bound')
                    slack = h - bound
                    minimum_slack = slack if minimum_slack is None else min(slack, minimum_slack)
                    nonzero += 1
                if m <= 3:
                    check(permutation_count(out) == h, 'permutation count')
                    perm_cases += 1
                if m <= 4:
                    check(odd_cycle_formula(out) == h, 'odd-cycle formula')
                    ocf_cases += 1
                for epsilon in (0, 1):
                    check_growth(m, q, tau, plus, epsilon, out)
                    growth_cases += 1
                cases += 1
                values.add(h)
        order_data.append({'m': m, 'order': m+4, 'distinct_total_path_counts': len(values),
                           'minimum_H': min(values), 'maximum_H': max(values),
                           'minimum_pair_bound_slack_over_nonzero_q': minimum_slack})
    cert_cases = cert_total = 0
    for j in range(10):
        for t in range(10 - j):
            m = j + t + 1
            out = tournament(m, 1 << j, 0, False)
            a, b = m + 2, m + 3
            seen = set()
            for A, B in certificates(j, t):
                C1 = [a, 1] + A
                C2 = [b, 0] + B + [j + 2]
                check(valid_cycle(out, C1) and valid_cycle(out, C2), 'invalid cycle certificate')
                check(not(set(C1) & set(C2)), 'certificate cycles intersect')
                key = (tuple(C1), tuple(C2))
                check(key not in seen, 'duplicate cycle certificate')
                seen.add(key)
            check(4 * len(seen) == pair_weight(j, t), 'pair parity formula')
            cert_cases += 1
            cert_total += len(seen)
    collision = []
    for m, q, expected in ((3, 6, 81), (3, 7, 81), (4, 12, 225), (4, 14, 245)):
        out = tournament(m, q, 0, False)
        h = hamilton_count(out)
        ocf = odd_cycle_formula(out)
        check(h == expected and ocf == expected, 'growth collision count')
        collision.append({'m':m, 'q':q, 'tau':0, 'apex_arc':'b->a',
                          'out_rows':out, 'H_subset_dp':h, 'H_odd_cycle_formula':ocf})
    print(json.dumps({'status':'PASS', 'universe':'all B+/B- with m=1..8, every q, both tau values',
                      'family_cases':cases, 'nonzero_code_bound_checks':nonzero,
                      'complete_permutation_cases_orders_5_to_7':perm_cases,
                      'odd_cycle_collection_cases_orders_5_to_8':ocf_cases,
                      'old_arc_preservation_checks':growth_cases,
                      'pair_certificate_parameter_cases_j_plus_t_at_most_9':cert_cases,
                      'explicit_pair_certificates':cert_total,
                      'order_data':order_data, 'scalar_growth_collision':collision,
                      'scope':'Finite implementation checks only; all-order density-zero theorem is in PROOF.md.'},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
