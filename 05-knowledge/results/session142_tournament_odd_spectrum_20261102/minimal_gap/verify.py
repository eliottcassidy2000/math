#!/usr/bin/env python3
"""Portable exact checks for prime-lane asymmetry and the sharp 7p threshold.
No imports from a source repository, random input, numerical tolerances,
or assertions that disappear under python -O.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate(rows):
    n = len(rows)
    for i, row in enumerate(rows):
        check(0 <= row < 1 << n and not (row >> i & 1), 'invalid row')
        for j in range(i):
            check(((row >> j) & 1) + ((rows[j] >> i) & 1) == 1,
                  'not a tournament')


def count_end_dp(rows):
    n = len(rows)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            nxt = rows[v] & ~mask
            while nxt:
                b = nxt & -nxt
                nxt -= b
                w = b.bit_length() - 1
                dp[mask | b][w] += count
    return sum(dp[-1])


def count_recursive(rows):
    n = len(rows)
    @lru_cache(None)
    def follow(last, remaining):
        if not remaining:
            return 1
        return sum(follow(w, remaining ^ (1 << w))
                   for w in range(n)
                   if (remaining >> w & 1) and (rows[last] >> w & 1))
    return sum(follow(v, ((1 << n) - 1) ^ (1 << v)) for v in range(n)) if n else 1


def count_literal(rows):
    return sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(rows)-1))
               for p in permutations(range(len(rows))))


def decode(n, code):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if code >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            k += 1
    return rows


def strong(rows):
    n = len(rows)
    for v in range(n):
        seen = todo = 1 << v
        while todo:
            bit = todo & -todo
            todo -= bit
            w = bit.bit_length() - 1
            new = rows[w] & ~seen
            seen |= new
            todo |= new
        if seen != (1 << n) - 1:
            return False
    return True


def automorphism(rows, p):
    return all(((rows[i] >> j) & 1) == ((rows[p[i]] >> p[j]) & 1)
               for i in range(len(rows)) for j in range(len(rows)))


def circulant(n, bits):
    steps = [d if bit else n-d for d, bit in enumerate(bits, 1)]
    return [sum(1 << ((v+d) % n) for d in steps) for v in range(n)]


def substitute_marked(outer, inner):
    """Replace outer vertex 0 by the homogeneous module inner."""
    n, k = len(inner), len(outer)-1
    rows = list(inner) + [0] * k
    for v in range(n):
        for j in range(1, k+1):
            if outer[0] >> j & 1:
                rows[v] |= 1 << (n+j-1)
            else:
                rows[n+j-1] |= 1 << v
    for i in range(1, k+1):
        for j in range(1, k+1):
            if outer[i] >> j & 1:
                rows[n+i-1] |= 1 << (n+j-1)
    validate(rows)
    return rows


def induced(rows, vertices):
    return [sum(1 << j for j, b in enumerate(vertices) if rows[a] >> b & 1)
            for a in vertices]


def prime(n):
    return n > 1 and all(n % d for d in range(2, 1 + __import__('math').isqrt(n)))


def floor(n):
    if n % 3 == 0:
        return 3 * 5 ** ((n-3)//3)
    if n % 3 == 1:
        return 5 ** ((n-1)//3)
    return 9 * 5 ** ((n-5)//3)


def main():
    receipt = {'status': 'PASS', 'claims': 'finite checks; theorems proved in PROOF.md'}
    # These three bases plus f(n+3)=5f(n) prove the infinite inequalities.
    check(all(floor(n) > n for n in (4,5,6)), 'f(n)>n bases')
    check(all(floor(n) > 7*n for n in (11,12,13)), 'f(n)>7n bases')
    check(all(floor(n+3) == 5*floor(n) for n in range(3,300)), 'floor recurrence')
    check(all(5*n > n+3 for n in range(4,300)), 'linear propagation')
    receipt['floor_bases'] = {
        'strict_order': [[n, floor(n)] for n in (4,5,6)],
        'seven_times_order': [[n, floor(n)] for n in (11,12,13)]}
    # Independent literal counts for every labelled tournament through order5.
    small = []
    for n in range(1,6):
        prime_count = 0
        symmetric_prime = 0
        marked_4_checks = 0
        for code in range(1 << (n*(n-1)//2)):
            rows = decode(n, code)
            h = count_end_dp(rows)
            check(h == count_literal(rows), 'literal DP mismatch')
            if h > 3 and prime(h):
                prime_count += 1
                aut = sum(automorphism(rows,p) for p in permutations(range(n)))
                symmetric_prime += (aut != 1)
                check(aut == 1, 'small prime asymmetry')
            if n >= 4 and strong(rows):
                for v in range(n):
                    from itertools import combinations
                    check(any(strong(induced(rows,[v]+list(rest)))
                              for rest in combinations([w for w in range(n) if w != v],3)),
                          'marked strong4 test')
                    marked_4_checks += 1
        small.append({'order':n,'tournaments':1 << (n*(n-1)//2),
                      'prime_gt3_count':prime_count,
                      'nontrivial_prime_automorphisms':symmetric_prime,
                      'marked_strong4_checks':marked_4_checks})
    receipt['small_literal_controls'] = small
    banks = {}
    for n in (5,7):
        bank = []
        for bits in product((0,1),repeat=(n-1)//2):
            rows = circulant(n,bits)
            validate(rows)
            h = count_end_dp(rows)
            check(h == count_literal(rows), 'circulant literal mismatch')
            check(strong(rows), 'circulant strong')
            check(automorphism(rows,tuple((v+1)%n for v in range(n))), 'rotation')
            bank.append((bits,rows,h))
        expected = {15} if n == 5 else {175,189}
        check({h for _,_,h in bank} == expected, 'circulant spectrum')
        banks[n] = bank
    receipt['circulants'] = {
        str(n): [{'bits':bits,'H':h} for bits,_,h in banks[n]] for n in (5,7)}
    tables = {}
    for n, m in ((5,3),(7,3),(7,4)):
        histogram = Counter()
        strong_outer = 0
        for code in range(1 << (m*(m-1)//2)):
            outer = decode(m,code)
            if not strong(outer):
                continue
            strong_outer += 1
            for bits,inner,_ in banks[n]:
                rows = substitute_marked(outer,inner)
                check(strong(rows),'strong substitution')
                h = count_end_dp(rows)
                check(h == count_recursive(tuple(rows)), 'independent recursion mismatch')
                histogram[h] += 1
        tables[f'{n}_orbit_{m}_quotient'] = {
            'strong_quotients':strong_outer,
            'substitutions':sum(histogram.values()),
            'H_histogram':dict(sorted(histogram.items()))}
        if (n,m) == (5,3):
            check(dict(histogram) == {135:8}, '5 triangle blowup')
        elif (n,m) == (7,3):
            check(dict(histogram) == {2163:12,2205:4}, '7 triangle blowup')
        else:
            check(dict(histogram) == {4151:72,4221:24,5789:72,5859:24}, 'marked4 bank')
    receipt['module_tables'] = tables
    check(175**2 > 4151, 'two seven-orbit barrier')
    check(not prime(175//7) and not prime(189//7), 'orbit quotient composite')
    check(not prime(2163//7) and not prime(2205//7), 'triangle quotient composite')
    witness = [142,156,184,240,225,195,135,512,255,383]
    validate(witness)
    check(strong(witness),'threshold witness strong')
    check(count_end_dp(witness) == count_recursive(tuple(witness)) == 4151, 'threshold H')
    check(prime(593) and 7*593 == 4151,'threshold prime')
    rotations = [tuple((v+k)%7 if v<7 else v for v in range(10)) for k in range(7)]
    check(all(automorphism(witness,p) for p in rotations), 'seven automorphisms')
    # H=7*593 and order10<593 exclude all other automorphism factors,
    # so the witnessed seven elements are the entire automorphism group.
    receipt['sharp_witness'] = {'rows':witness,'H':4151,'prime_quotient':593,
                                'visible_automorphisms':7,'exact_group_order_by_proof':7}
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
