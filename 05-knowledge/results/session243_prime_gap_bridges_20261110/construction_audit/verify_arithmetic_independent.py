#!/usr/bin/env python3
"""Independent supplementary audit of final wheel, arithmetic and root witnesses."""
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
from math import gcd, factorial, lcm
import json
import sys
from verify_independent import require, decode, profile, ab, gadget, ordinal, transitive, vertices, validate
from verify_gaps_independent import compact_words, word_state


def prime(n):
    if n < 2: return False
    p = 2
    while p * p <= n:
        if n % p == 0: return False
        p += 1 if p == 2 else 2
    return True


def hfull(rows):
    @lru_cache(None)
    def path(mask, v):
        rest = mask ^ (1 << v)
        return 1 if not rest else sum(path(rest, u) for u in vertices(rest & rows[v]))
    full = (1 << len(rows)) - 1
    return sum(path(full, v) for v in range(len(rows)))


def main():
    if '--fail' in sys.argv: require(False, 'intentional failure')
    primes = [p for p in range(3, 32) if prime(p)]
    expected = {3: (2, []), 5: (4, [0, 2]), 7: (3, [0]), 11: (10, [1, 6]),
                13: (12, [1, 7]), 17: (8, []), 19: (18, [4, 13]),
                23: (11, [2]), 29: (28, [8, 22]), 31: (5, [])}
    wheel = []
    period = 1
    for p in primes:
        power, order = 2 % p, 1
        while power != 1:
            power = power * 2 % p; order += 1
        zeros = [t for t in range(order) if any((6 * pow(2, t, p) + a) % p == 0 for a in [-1, 1])]
        require((order, zeros) == expected[p], ('wheel row', p, order, zeros))
        period = lcm(period, order)
        wheel.append({'prime': p, 'order': order, 'forbidden': zeros})
    survivors = []
    for t in range(period):
        direct = all((6 * pow(2, t, p) + a) % p for p in primes for a in [-1, 1])
        compact = t % 6 == 5 and t % 5 != 1 and t % 11 != 2
        require(direct == compact, 'reduced wheel iff')
        if t < 330 and direct: survivors.append(t)
    require(len(survivors) == 40, 'density')
    twins, four, six, consecutive_six = [], [], [], []
    samples = []
    for m in range(1, 31):
        vals = [6 * 2 ** m + a for a in [-1, 1, 3, 5]]
        if prime(vals[0]) and prime(vals[1]): twins.append(m)
        if prime(vals[1]) and prime(vals[3]): four.append(m)
        if prime(vals[0]) and prime(vals[3]):
            six.append(m)
            if not any(prime(x) for x in range(vals[0] + 1, vals[3])): consecutive_six.append(m)
        samples.append({'m': m, 'values': vals, 'prime_flags': [prime(x) for x in vals]})
    require(twins == [1, 5, 17] and four == [1, 4, 5, 7] and six == [1, 2, 3, 5, 6], 'prime lists')
    require(consecutive_six == [2, 3, 6], 'consecutive gap six')
    require(6563 * 7669 == 6 * 2 ** 23 - 1 and 61 * 825109 == 6 * 2 ** 23 + 1, 'surviving composite factorization')
    E6 = list(transitive(6)); E6[0] ^= 1 << 5; E6[5] |= 1
    require(gadget(E6) == (17, 224), 'E6 pair')
    require([6 * 224 + a * 17 for a in [-1, 1, 3, 5]] == [1327, 1361, 1395, 1429], 'E6 window')
    brown = []
    for rows, H, n in [((6, 12, 8, 1), 5, 4), ((24, 17, 19, 6, 8), 11, 5),
                       ((206, 252, 248, 240, 225, 65, 128, 32), 71, 7)]:
        validate(rows)
        direct = sum(all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))
                     for p in permutations(range(len(rows))))
        square = ordinal(rows, rows)
        got_square = hfull(square)
        require(direct == H == hfull(rows), 'Brown root')
        require(got_square == H * H == factorial(n) + 1, 'Brown square')
        brown.append({'rows': rows, 'H': H, 'joined_H': got_square, 'joined_order': len(square)})
    states, pairs, roots = set(), set(), 0
    for n in range(1, 6):
        for k in range(1 << (n * (n - 1) // 2)):
            rows = decode(n, k)
            for root in range(n):
                if not rows[root]: continue
                roots += 1
                s = profile(rows, root)
                states.add(s); pairs.add(ab(s))
    require((roots, len(states), len(pairs)) == (5044, 43, 39), 'full source census')
    for N in range(1, 10001):
        _, upper = compact_words(N)
        got = word_state((1, 0, 2, 1), upper)[0]
        equality = N == 1 or (N >= 3 and ((N - 1) & (N - 2)) == 0)
        require((got == 5 * N) == equality, 'compact exact equality classification')
    rows4 = {}
    basis = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
    for n in range(5):
        for letters in product('OZJ', repeat=n):
            word = ''.join(letters)
            rows4[word] = tuple(word_state(x, word)[0] for x in basis)
    rowset = set(rows4.values())
    require(not any((r[0], r[1]+2, r[2], r[3]) in rowset for r in rowset), 'short endpoint suffix absence')
    print(json.dumps({'status': 'PASS', 'all_words_length_at_most4': len(rows4), 'short_endpoint_gap_pairs': 0, 'wheel': wheel, 'wheel_full_lcm': period,
          'compact_period': 330, 'compact_survivors': survivors, 'prime_examples': samples,
          'prime_exponents': {'twins': twins, 'gap4': four, 'gap6': six, 'consecutive_gap6': consecutive_six},
          'E6': {'rows': E6, 'pair': [17, 224]}, 'Brown': brown,
          'source_census': {'nonsink_roots': roots, 'full_states': len(states), 'ab_pairs': len(pairs)},
          'compact_equality_targets': 10000}, sort_keys=True, indent=2))

if __name__ == '__main__': main()
