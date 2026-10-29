#!/usr/bin/env python3
"""Finite-exact independent controls; no floating point counting."""
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial, gcd
import json
from pathlib import Path

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def oddparts(n, least=1):
    if not n:
        yield ()
    for k in range(least, n + 1, 2):
        for tail in oddparts(n - k, k):
            yield ((k,) + tail)

def count(n):
    acc = F(0)
    for parts in oddparts(n):
        mult = Counter(parts)
        z = 1
        for k, m in mult.items():
            z *= k ** m * factorial(m)
        e = sum(((k - 1) // 2 for k in parts))
        e += sum((gcd(x, y) for i, x in enumerate(parts) for y in parts[i + 1:]))
        support = sum((k for k in parts if k > 1))
        loss = comb(n, 2) - e
        require(3 * loss >= support * (2 * n - support - 1), '3 * loss >= support * (2 * n - support - 1)')
        if support == 3:
            require(loss == 2 * n - 4, 'loss == 2 * n - 4')
        require(support not in (1, 2, 4), 'support not in (1, 2, 4)')
        acc += F(2 ** e, z)
    require(acc.denominator == 1, 'acc.denominator == 1')
    return int(acc)
N = 30
a = [count(n) for n in range(N + 1)]
require(a[:16] == [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248, 154108311168, 48542114686912, 28401423719122304, 31021002160355166848], 'a[:16] == [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248, 154108311168, 48542114686912, 28401423719122304, 31021002160355166848]')
s = [0] * (N + 1)
for n in range(1, N + 1):
    s[n] = a[n] - sum((s[k] * a[n - k] for k in range(1, n)))
require(s[:8] == [0, 1, 0, 1, 1, 6, 35, 353], 's[:8] == [0, 1, 0, 1, 1, 6, 35, 353]')
b = [F(factorial(n) * a[n], 2 ** comb(n, 2)) for n in range(N + 1)]
d = [F(factorial(n) * s[n], 2 ** comb(n, 2)) for n in range(N + 1)]

def kernel(n, k):
    return F(comb(n, k), 2 ** (k * (n - k)))
for n in range(1, N + 1):
    require(b[n] == sum((kernel(n, k) * d[k] * b[n - k] for k in range(1, n + 1))), 'b[n] == sum((kernel(n, k) * d[k] * b[n - k] for k in range(1, n + 1)))')
for n in range(N + 1):
    for i in range(n + 1):
        for j in range(n - i + 1):
            k = n - i - j
            require(kernel(n, i + j) * kernel(i + j, i) == kernel(n, i) * kernel(j + k, j), 'kernel(n, i + j) * kernel(i + j, i) == kernel(n, i) * kernel(j + k, j)')
    leading = F(32 * comb(n, 3), 4 ** n) if n >= 3 else F(0)
    require(b[n] - 1 - leading >= 0, 'b[n] - 1 - leading >= 0')
require(b[3] != sum((d[k] * b[3 - k] for k in range(1, 4))), 'b[3] != sum((d[k] * b[3 - k] for k in range(1, 4)))')
report = {'status': 'PASS', 'n_max': N, 'a_initial': a[:16], 's_initial': s[:11], 'checks': ['Burnside integers', 'known initial counts', 'support loss bound', 'SCC convolution', 'twisted normalization', 'associative kernel', 'nonnegative correction remainder'], 'normalized_diagnostics': [{'n': n, 'b_exact': str(b[n]), 'b_minus_one': str(b[n] - 1)} for n in [3, 5, 10, 20, 30]]}
Path(__file__).with_name('verification.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
