#!/usr/bin/env python3
"""Finite controls of the proved sharp boundary and the two unbounded signs.
The full <=11 theorem is a proof, not inferred from this q<=4 test universe.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json

CHECKS = 0

def check(ok, text):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(text)

def compositions(n, q):
    if q == 1:
        if n > 0:
            yield (n,)
        return
    for first in range(1, n-q+2):
        for tail in compositions(n-first, q-1):
            yield (first,) + tail

def calc(E, c):
    q = len(c)
    L = [p for p in permutations(range(q))
         if all(c[i] == c[p[i]] for i in range(q))
         and {tuple(sorted((p[i], p[j]))) for i, j in E} == E]
    signs = [(-1) ** sum(c[i]*c[j] for i, j in E if p[i] > p[j]) for p in L]
    total = 0
    for p, sign in zip(L, signs):
        value, used = sign, set()
        for i in range(q):
            if i in used:
                continue
            j, length = i, 0
            while j not in used:
                used.add(j)
                length += 1
                j = p[j]
            value *= sum((-1) ** (length*a) for a in range(c[i]//3+1))
        total += value
    D = Fraction(total, len(L))
    check(D.denominator == 1, 'nonintegral D')
    even = all(s == 1 for s in signs)
    delta = -int(D) if max(c) >= 3 else 0
    if even:
        check(delta <= 0, 'even graph with positive defect')
    counts = Counter(c)
    if any((size//3) % 2 and count % 2 for size, count in counts.items()):
        check(D == 0, 'odd multiplicity vanishing test')
    return dict(D=int(D), delta=delta, globally_even=even)

def main():
    counts = Counter()
    exceptions = []
    for q in range(1, 5):
        pairs = list(combinations(range(q), 2))
        for mask in range(1 << len(pairs)):
            E = {e for i, e in enumerate(pairs) if mask >> i & 1}
            nb = [frozenset(v if u == i else u for u, v in E if i in (u, v))
                  for i in range(q)]
            if len(set(nb)) != q:
                continue
            degrees = sorted(map(len, nb))
            for n in range(q, 13):
                for c in compositions(n, q):
                    result = calc(E, c)
                    counts[n] += 1
                    d, even = result['delta'], result['globally_even']
                    if n <= 11:
                        check(abs(d) <= 1, 'sharp <=11 magnitude bound')
                        check(even or d >= 0, 'sharp <=11 odd-sign bound')
                    else:
                        large_exception = (q == 2 and c == (6, 6)) or (
                            q == 4 and c == (3, 3, 3, 3) and degrees == [1, 1, 2, 2])
                        odd_negative = q == 4 and c == (3, 3, 3, 3) and degrees == [1]*4
                        check((abs(d) > 1) == large_exception, 'order12 magnitude classification')
                        check((not even and d < 0) == odd_negative, 'order12 odd negative classification')
                        if large_exception or odd_negative:
                            exceptions.append(dict(q=q, mask=mask, weights=c, **result))
    family = []
    for r in range(1, 41):
        c = 6*r+3
        positive = calc({(0, 1)}, (c, c))
        negative = calc({(0, 1), (2, 3)}, (3, 3, c, c))
        check(positive['delta'] == r+1, 'positive odd family')
        check(negative['delta'] == -(r+1), 'negative odd family')
        check(not positive['globally_even'] and not negative['globally_even'], 'odd family sign')
        family.append(dict(r=r, positive=positive, negative=negative))
    print(json.dumps(dict(status='PASS', checks=CHECKS,
                          finite_scope='All false-twin-free labelled quotients of orders 1..4, all positive integer weights with total at most 12; two odd-sign families for r=1..40.',
                          weighted_quotient_cases=sum(counts.values()),
                          cases_by_total_order=dict(sorted(counts.items())),
                          order12_labelled_exceptions=exceptions,
                          opposite_sign_unbounded_family=family), indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
