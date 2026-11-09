import sys

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

if "--false-gate" in sys.argv:
    check(False, "deliberate-false gate sentinel")

from functools import lru_cache
from itertools import product
import json
from pathlib import Path
C = [4, 1, 18, 7, 75, 31, 175, 63]
cut = 124

def hc(rows):
    n = len(rows)
    full = (1 << n) - 1

    @lru_cache(None)
    def paths(v, remaining):
        if not remaining:
            return 1
        candidates = remaining & rows[v]
        total = 0
        while candidates:
            bit = candidates & -candidates
            candidates -= bit
            total += paths(bit.bit_length() - 1, remaining ^ bit)
        return total
    return sum((paths(v, full ^ 1 << v) for v in range(n)))

def literal(prefix, suffix):
    a = len(prefix)
    bits = list(prefix) + [cut >> i & 1 for i in range(8)] + list(suffix)
    n = len(bits)
    rows = []
    for i in range(n):
        row = sum((1 << j for j in range(i + 1, n)))
        if a <= i < a + 8:
            row = row & ~(255 << a) | C[i - a] << a
        rows.append(row | 1 - bits[i] << n)
    rows.append(sum((b << i for i, b in enumerate(bits))))
    return hc(rows)

def transfer(prefix, suffix):
    a = t = e = 0
    for c in prefix:
        a, t, e = (2 * a + c, t + c * (a + e), c)
    a, t, e = (1440 * a + 878, 577 * a + 45 * t + 332 + 25 * e, 14)
    for c in suffix:
        a, t, e = (2 * a + 45 * c, t + c * (a + e), 45 * c)
    return a - t + 45 - e
cases = 0
for total in range(5):
    for a in range(total + 1):
        for bits in product((0, 1), repeat=total):
            p, s = (bits[:a], bits[a:])
            check(literal(p, s) == transfer(p, s), (p, s))
            pn = p
            while pn and pn[0] == 0:
                pn = pn[1:]
            sn = s
            while sn and sn[-1] == 1:
                sn = sn[:-1]
            check(transfer(p, s) == transfer(pn, sn), 'correctness gate failed')
            cases += 1

def pair(w):
    return sum((2 ** (j - i - 1) + (j == i + 1) for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] and w[j]))

def integer(w):
    return sum((v * 2 ** (len(w) - 1 - i) for i, v in enumerate(w)))
for total in range(9):
    for a in range(total + 1):
        for bits in product((0, 1), repeat=total):
            p, s = (bits[:a], bits[a:])
            P = integer(p)
            R = integer(s)
            b = len(s)
            f = (1440 * P + 878) * (2 ** b - integer(s[::-1])) + 45 * (R + 1 - pair(s) - (s[-1] if s else 0)) - 577 * P - 45 * pair(p) - 332 - 25 * (p[-1] if p else 0) - 14 * (s[0] if s else 1)
            check(f == transfer(p, s), 'correctness gate failed')
if __name__ == '__main__':
 print(json.dumps({'literal_cases': cases, 'scalar_max_context_length': 8, 'normalization': 'PASS', 'literal': 'PASS', 'scalar': 'PASS'}))
