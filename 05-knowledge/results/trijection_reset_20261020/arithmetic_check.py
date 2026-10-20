"""Independent finite controls for the divisor decomposition note."""
from collections import Counter
from itertools import combinations_with_replacement, permutations, product
from math import isqrt
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def divisors(n):
    out = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.update((d, n // d))
    return out - {1, n}


def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def fsu(n):
    profiles = [factor(d) for d in divisors(n)]
    return (len(profiles), sum(max(f.values()) == 1 for f in profiles),
            sum(sum(f.values()) == 1 for f in profiles))


def sign(x):
    return (x > 0) - (x < 0)


rows = []
numeric_checks = 0
for k in range(2, 41):
    states = set(product(range(k + 1), range(2), range(2))) - {
        (0, 0, 0), (k, 1, 1)}
    C = lambda x: (k - x[0], 1 - x[1], 1 - x[2])
    Q = lambda x: (x[0], x[2], x[1])
    f = lambda x: 2 * x[0] - k * (x[1] + x[2])
    g = lambda x: x[1] - x[2]
    counts_f = Counter(sign(f(x)) for x in states)
    counts_g = Counter(sign(g(x)) for x in states)
    fixed = [sum(op(x) == x for x in states)
             for op in (C, Q, lambda x: C(Q(x)))]
    expected_f = (2*k, 2, 2*k) if k % 2 == 0 else (2*k+1, 0, 2*k+1)
    require(tuple(counts_f[s] for s in (-1, 0, 1)) == expected_f, (k, "f"))
    require(tuple(counts_g[s] for s in (-1, 0, 1)) == (k+1, 2*k, k+1), (k, "g"))
    require(fixed == [0, 2*k, 2 if k % 2 == 0 else 0], (k, "fixed"))
    chart = {(f(x), g(x)): x for x in states}
    require(len(chart) == len(states) and (0, 0) not in chart, (k, "faithful chart"))
    for (height, row), x in chart.items():
        if row:
            recovered = ((height + k)//2, int(row == 1), int(row == -1))
        elif height > 0:
            recovered = (height//2, 0, 0)
        else:
            recovered = ((height+2*k)//2, 1, 1)
        require(recovered == x, (k, x, "chart inverse"))
    full = set(product(range(k + 1), range(2), range(2)))
    fibres = Counter((f(x), g(x)) for x in full)
    require(fibres[(0, 0)] == 2 and all(n == 1 for y, n in fibres.items() if y != (0, 0)), (k, "endpoint collision"))
    for x in states:
        require(C(x) in states and Q(x) in states, (k, x, "domain"))
        require(C(C(x)) == Q(Q(x)) == x and C(Q(x)) == Q(C(x)), (k, x, "group"))
        require(f(C(x)) == -f(x) and f(Q(x)) == f(x), (k, x, "f law"))
        require(g(C(x)) == g(Q(x)) == -g(x), (k, x, "g law"))
        require((f(x) == 0) == (C(Q(x)) == x), (k, x, "zero"))
    if k <= 8:
        for p, q, r in ((2, 3, 5), (3, 2, 5), (5, 3, 2)):
            n = p**k * q * r
            values = {x: p**x[0] * q**x[1] * r**x[2] for x in states}
            require(set(values.values()) == divisors(n), (n, "numeric divisor set"))
            for x, d in values.items():
                require(values[C(x)] == n // d, (n, d, "complement"))
                for y, e in values.items():
                    require((e % d == 0) == all(a <= b for a, b in zip(x, y)), (n, d, e, "order"))
                    require((e % d == 0) == (values[C(x)] % values[C(y)] == 0), (n, d, e, "reverse"))
            require(fsu(n) == (4*k+2, 7, 3), (n, "fsu"))
            numeric_checks += 1
    if k <= 6:
        rows.append({"k": k, "sign_f": expected_f, "sign_g": (k+1, 2*k, k+1), "fixed_C_Q_CQ": fixed})

require(fsu(60) == (10, 7, 3) and fsu(8) == (2, 1, 1), "positive controls")
require(fsu(4) == (1, 1, 1) and fsu(240) == (18, 7, 3), "hostile controls")
for exponent in range(2, 13):
    require(fsu(2**exponent * 3) == (2*exponent, 3, 2), (exponent, "literal shape"))

pairs = tuple(combinations_with_replacement(range(4), 2))
fixed_pair_histogram = Counter()
for permutation in permutations(range(4)):
    if all(permutation[permutation[i]] == i for i in range(4)):
        fixed_count = sum(tuple(sorted((permutation[a], permutation[b]))) == (a, b) for a, b in pairs)
        fixed_pair_histogram[fixed_count] += 1
require(dict(fixed_pair_histogram) == {10: 1, 4: 6, 2: 3}, "geometric fixed-point obstruction")

print(json.dumps({"status": "FINITE-EXACT", "numeric_prime_assignments": numeric_checks,
                  "exponent_universe": [2, 40], "first_rows": rows,
                  "tetrahedral_involutions_by_fixed_pairs": dict(sorted(fixed_pair_histogram.items()))},
                 sort_keys=True, indent=2))
