"""Literal cover/color controls for the slot-polynomial precision theorem."""
from pathlib import Path
from itertools import product
from math import comb, factorial, gcd
from functools import reduce
import hashlib
import json

HERE = Path(__file__).resolve().parent
CORE = (56, 117, 129, 22, 100, 76, 141, 59)
EXPECTED = (513, 6102, 32508, 99432, 187200, 216000, 141120, 40320)


def path_counts(rows):
    n = len(rows)
    ends = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        ends[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            previous = mask ^ (1 << v)
            if previous:
                ends[mask][v] = sum(ends[previous][u] for u in range(n)
                                     if previous >> u & 1 and rows[u] >> v & 1)
    return [1] + [sum(e) for e in ends[1:]]


def cover_counts(rows):
    """Ordered subset convolution; independent of the slot coefficients."""
    n = len(rows)
    h = path_counts(rows)
    covers = [[0] * (n + 1) for _ in range(1 << n)]
    covers[0][0] = 1
    for mask in range(1, 1 << n):
        sub = mask
        while sub:
            rest = mask ^ sub
            for k in range(1, mask.bit_count() + 1):
                covers[mask][k] += h[sub] * covers[rest][k - 1]
            sub = (sub - 1) & mask
    return h, tuple(covers[-1][1:])


def coefficients(covers):
    n = len(covers)
    out = [0] * (n + 1)
    falling = [1]
    for k, c in enumerate(covers, 1):
        if c % factorial(k):
            raise ValueError('unordered cover integrality')
        falling = [0] + falling
        for j in range(k):
            falling[j] -= (k - 1) * falling[j + 1]
        for j, value in enumerate(falling):
            out[j] += (c // factorial(k)) * value
    return out


def eval_poly(coeff, x):
    out = 0
    for a in reversed(coeff):
        out = out * x + a
    return out


def color_count(h, n, q):
    result = 0
    for colors in product(range(q), repeat=n):
        masks = [0] * q
        for v, c in enumerate(colors):
            masks[c] |= 1 << v
        term = 1
        for mask in masks:
            term *= h[mask]
        result += term
    return result


def valuation(n, p):
    if n == 0:
        return None
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def main():
    h, covers = cover_counts(CORE)
    if covers != EXPECTED:
        raise ValueError('literal cover profile')
    coeff = coefficients(covers)
    if coeff != [0, 0, 246, 0, 232, 0, 34, 0, 1]:
        raise ValueError('power basis')
    colors = []
    for q in range(1, 5):
        actual = color_count(h, len(CORE), q)
        if actual != eval_poly(coeff, q):
            raise ValueError(('color evaluation', q))
        colors.append({'slots': q, 'literal_colorings': q ** 8, 'value': actual})
    for x in range(-200, 201):
        value = eval_poly(coeff, x)
        if value % 27 != (18 if x % 9 in (4, 5) else 0):
            raise ValueError(('residue', x))
        if x and x % 3 == 0 and valuation(value, 3) != 2 * valuation(x, 3) + 1:
            raise ValueError(('nonunit valuation', x))
    padding = []
    for m in (0, 1, 2):
        rows = tuple(row | sum(1 << j for j in range(8, 8 + m)) for row in CORE)
        rows += tuple(sum(1 << j for j in range(i + 1, 8 + m)) for i in range(8, 8 + m))
        _, c = cover_counts(rows)
        power = coefficients(c)
        if power != [0] * m + coeff:
            raise ValueError(('padding polynomial', m))
        raw = reduce(gcd, c)
        context = reduce(gcd, (factorial(k - 1) * a for k, a in enumerate(c, 1)))
        if (raw, context) != (9, 27):
            raise ValueError(('padding content', m))
        padding.append({'suffix': m, 'order': len(rows), 'h': raw, 'g': context})
    precision_controls = 0
    for p in (2, 3, 5):
        for s in (1, 2):
            for k in (1, 2, 3, 4, 8, 9, 12):
                e, power = 0, p
                while power <= k:
                    e += 1
                    power *= p
                period = p ** (s + e)
                for x in range(40):
                    if (comb(x + period, k) - comb(x, k)) % (p ** s):
                        raise ValueError(('period', p, s, k, x))
                    precision_controls += 1
    receipt = {
        'status': 'FINITE-EXACT controls; universal conclusions proved in PROOF.md',
        'core': CORE, 'ordered_covers': covers, 'power_coefficients': coeff,
        'independent_color_controls': colors,
        'residue_integer_interval': [-200, 200],
        'normalized_mod3_period': [eval_poly(coeff, x) // 9 % 3 for x in range(9)],
        'literal_padding_controls': padding,
        'binomial_period_controls': precision_controls,
        'hostile_nonunit_only_explanation': {'unit_input': 4, 'value': eval_poly(coeff, 4), 'v3': 2},
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'color_slots': [1, 2, 3, 4],
                      'residue_inputs': 401, 'padding_orders': [8, 9, 10],
                      'period_controls': precision_controls}))


if __name__ == '__main__':
    main()
