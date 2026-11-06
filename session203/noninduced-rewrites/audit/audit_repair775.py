#!/usr/bin/env python3
"""Independent literal construction and path DFS for the E129 repair."""
import itertools
import json
import math
from collections import Counter

def need(b, s):
    if not b:
        raise RuntimeError(s)

def construct(k):
    # Labels: r,c,p0,...,pk,x,y. Build by pair choices, not copied rows.
    n = k + 5
    x, y = n - 2, n - 1
    a = [0] * n
    def arc(i, j):
        a[i] |= 1 << j
    arc(0, 1)
    for p in range(2, k + 3):
        arc(p, 0)
        # The original first bit is 1, all following bits 0.
        # Reverse c -> p0, so now every private vertex beats c.
        arc(p, 1)
        for q in range(p + 1, k + 3):
            arc(p, q)
    for old in range(x):
        arc(x, old) if old in (0, 1, 3) else arc(old, x)
    for old in range(y):
        arc(y, old) if old in (0, 1, 4, 6) else arc(old, y)
    for i, row in enumerate(a):
        need(not row & (1 << i), 'loop')
        for j in range(i):
            need(((row >> j) & 1) + ((a[j] >> i) & 1) == 1, 'pair')
    return a

def literal_count(a):
    full = (1 << len(a)) - 1
    calls = 0
    def visit(last, mask):
        nonlocal calls
        calls += 1
        if mask == full:
            return 1
        total = 0
        choices = a[last] & ~mask
        while choices:
            bit = choices & -choices
            choices -= bit
            total += visit(bit.bit_length() - 1, mask | bit)
        return total
    result = sum(visit(v, 1 << v) for v in range(len(a)))
    return result, calls

def main():
    core = construct(4)
    outgoing = {0, 1, 7, 8}
    # Check literal clone relation and the core embedding over extra lengths.
    controls = []
    for k in range(4, 11):
        a = construct(k)
        emb = list(range(7)) + [k + 3, k + 4]
        need([[bool(a[i] & (1 << j)) for j in emb] for i in emb] ==
             [[bool(core[i] & (1 << j)) for j in range(9)] for i in range(9)], 'core embedding')
        for p in range(7, k + 3):
            need({i for i, q in enumerate(emb) if a[p] & (1 << q)} == outgoing, 'clone type')
            for q in range(p + 1, k + 3):
                need(a[p] & (1 << q), 'clone chain')
        h, calls = literal_count(a)
        need(h == 26 * 3 ** (k - 4) + 10 * 2 ** (k - 4) - 7, 'formula control')
        controls.append({'k': k, 'H': h, 'target': 6 * ((1 << k) + 1) + 1,
                         'literal_DFS_calls': calls})
    # Count core orders by (number of allowed clone gaps, mandatory gaps).
    histogram = Counter()
    examined = 0
    for p in itertools.permutations(range(9)):
        examined += 1
        allowed = int(p[0] in outgoing) + int(p[-1] not in outgoing)
        mandatory = 0
        good = True
        for left, right in zip(p, p[1:]):
            possible = left not in outgoing and right in outgoing
            allowed += int(possible)
            if not core[left] & (1 << right):
                if not possible:
                    good = False
                    break
                mandatory += 1
        if good:
            histogram[allowed, mandatory] += 1
    coefficients = [0] * 11
    for (allowed, mandatory), multiplicity in histogram.items():
        for j in range(mandatory + 1):
            coefficients[allowed - j] += multiplicity * (-1) ** j * math.comb(mandatory, j)
    need(coefficients == [0, -7, 10, 26] + [0] * 7, 'exact coefficients')
    need(construct(7) == [2, 0, 4091, 3059, 2019, 4035, 1923, 3843, 3587, 3075, 2059, 83], 'published rows')
    print(json.dumps({'status': 'PASS', 'core_rows': core, 'repair775_rows': construct(7),
                      'clone_outgoing_core_indices': sorted(outgoing),
                      'controls': controls, 'core_orders_checked': examined,
                      'allowed_mandatory_histogram': {str(k): v for k, v in sorted(histogram.items())},
                      'exact_exponential_coefficients': coefficients}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
