"""Exact controls for funnel determinant/growth and the 462 rewrite debt."""
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def coefficients(word, m):
    """Inherited six-state tail convention, stated explicitly in the proof."""
    length = len(word)
    left = (m - 1) ** length + sum(m ** i * (m - 1) ** (length - i - 1)
                                     for i, bit in enumerate(word) if bit)
    down = (m - 1) ** length + sum((m - 1) ** i * m ** (length - i - 1)
                                     for i, bit in enumerate(word) if not bit)
    states = {(0, 0): 1}
    for bit in word:
        nxt = {}
        for (last, started), weight in states.items():
            for key, factor in [((last, started), m - 1),
                                ((1 if bit else 2, started), 1)] + (
                                    [((last, 1), 1)] if started or not bit else []):
                nxt[key] = nxt.get(key, 0) + factor * weight
        states = nxt
    return left, down, sum(v for (last, _), v in states.items() if last != 2)


def matrix(word):
    l2, d2, b2 = coefficients(word, 2)
    l3, d3, _ = coefficients(word, 3)
    return ((2 * b2, -2 * d2, 2 * l3), (l2, -1, 0),
            (d3, 0, 3 ** len(word)))


def det(m):
    a, b, c = m
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def val3(n):
    need(n != 0, 'zero determinant')
    out = 0
    while n % 3 == 0:
        n //= 3
        out += 1
    return out


def leading_ones(word):
    return next((i for i, bit in enumerate(word) if not bit), len(word))


def advance(z, word):
    h, e, s, p, f, t = z
    l1, d1, b1 = coefficients(word, 1)
    l2, d2, b2 = coefficients(word, 2)
    l3, d3, _ = coefficients(word, 3)
    return (l2 * s + b1 * h - d1 * e, l1 * h - int(not word) * e,
            2 ** len(word) * s + d2 * h,
            2 * (l3 * t + b2 * p - d2 * f), l2 * p - f,
            3 ** len(word) * t + d3 * p)


def funnel(rows, root, word):
    n, length = len(rows), len(word)
    newborn = n + length
    tail = ((1 << length) - 1) << n
    out = [row | tail | (0 if i == root else 1 << newborn)
           for i, row in enumerate(rows)]
    for i, bit in enumerate(word):
        out.append(sum(1 << (n + j) for j in range(i + 1, length))
                   | (1 << newborn if bit else 0))
    out.append((1 << root) | sum(1 << (n + i) for i, bit in enumerate(word) if not bit))
    return tuple(out), newborn


def profile(rows, root):
    n = len(rows)
    full = (1 << n) - 1

    @lru_cache(None)
    def end(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(end(rest, u) for u in range(n)
                   if rest >> u & 1 and rows[u] >> v & 1)

    @lru_cache(None)
    def start(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(start(rest, u) for u in range(n)
                   if rest >> u & 1 and rows[v] >> u & 1)

    @lru_cache(None)
    def h(mask):
        return sum(end(mask, v) for v in range(n) if mask >> v & 1) if mask else 1

    @lru_cache(None)
    def two(mask):
        sub, out = mask, 0
        while True:
            out += h(sub) * h(mask ^ sub)
            if not sub:
                return out
            sub = (sub - 1) & mask

    e = end(full, root)
    s2 = f2 = s3 = 0
    for mask in range(1 << n):
        if mask >> root & 1:
            rest = full ^ mask
            s2 += start(mask, root) * h(rest)
            f2 += end(mask, root) * h(rest)
            s3 += start(mask, root) * two(rest)
    return h(full), e, s2, two(full), f2, s3


def main():
    words = [()] + [w + (1,) for length in range(1, 11)
                    for w in product((0, 1), repeat=length - 1)]
    matrix_controls = []
    for word in words:
        r = leading_ones(word)
        l2, d2, b2 = coefficients(word, 2)
        l3, d3, _ = coefficients(word, 3)
        determinant = det(matrix(word))
        need(determinant == 2 * ((d2 * l2 - b2) * 3 ** len(word) + l3 * d3),
             'determinant formula')
        need(val3(determinant) == r, 'initial-one valuation')
        need(2 * b2 - d2 >= 1 and l3 >= 3 ** r, 'positive growth constants')
        if len(word) <= 3:
            matrix_controls.append({'word': ''.join(map(str, word)), 'matrix': matrix(word),
                                    'determinant': determinant, 'v3': r,
                                    'growth_lower_multiplier': l3 + 2 * b2 - d2})

    source_checks = child_checks = 0
    for n in range(1, 5):
        pairs = list(combinations(range(n), 2))
        for bits in range(1 << len(pairs)):
            rows = [0] * n
            for i, (u, v) in enumerate(pairs):
                a, b = (u, v) if bits >> i & 1 else (v, u)
                rows[a] |= 1 << b
            rows = tuple(rows)
            for root in range(n):
                z = profile(rows, root)
                h, _, s, p, f, t = z
                need(p % 2 == 0 and 2 * f <= p <= 2 * t and h <= s <= p,
                     'cover inequalities')
                d, j, k = s - h, t - p // 2, p // 2 - f
                need(j >= d >= 0 and k >= 0, 'bad-prefix injections')
                debt = 18 * h - 18 * s - 14 * p - 2 * f + 30 * t
                need(debt == 2 * (15 * j - 9 * d + k), '462 exact offset debt')
                need((debt == 0) == (rows[root] == 0), '462 sink equality boundary')
                source_checks += 1
                for word in words[:4]:  # empty,1,01,11; complete legal words of length<=2
                    child, marked = funnel(rows, root, word)
                    actual = profile(child, marked)
                    need(actual == advance(z, word), 'literal six-state funnel map')
                    l3 = coefficients(word, 3)[0]
                    need(actual[3] >= (l3 + 1) * p, 'literal pruning bound')
                    child_checks += 1

    chains = []
    for address in [('1',), ('11', '1'), ('1', '11', '')]:
        z = (1, 1, 1, 2, 1, 1)
        for word in address:
            z = advance(z, tuple(map(int, word)))
        chains.append({'address': address, 'six_state': z})
    need([x['six_state'] for x in chains] ==
         [(3, 1, 3, 12, 3, 7), (15, 5, 15, 216, 57, 115),
          (27, 12, 30, 564, 171, 327)], 'inherited 462 core states')
    slope_row, offset_row = [], []
    for i in range(6):
        basis = tuple(int(j == i) for j in range(6))
        left = advance(advance(basis, (1, 1)), (1,))
        right = advance(advance(advance(basis, (1,)), (1, 1)), ())
        slope_row.append(7 * (right[2] + right[3]) - 18 * (left[2] + left[3]))
        offset_row.append(5 * right[0] - right[4] - 3 * (3 * left[0] - left[4]))
    need(slope_row == [16, 0, -76, 312, -102, -462], '462 slope row')
    need(offset_row == [18, 0, -18, -14, -2, 30], '462 offset row')
    out = {'status': 'FINITE-EXACT controls for the all-word proofs',
           'legal_words_through_length': 10, 'word_count': len(words),
           'source_universe': 'Every labelled tournament through order4, every root',
           'source_checks': source_checks, 'literal_funnel_checks': child_checks,
           'small_matrices': matrix_controls, 'inherited_462_core_states': chains,
           '462_inner_substitution_debt': {'slope_row': slope_row, 'offset_row': offset_row},
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(out, indent=2) + '\n'
    destination = HERE / 'receipt.json'
    if '--write' in sys.argv:
        destination.write_text(encoded)
    else:
        need(destination.read_text() == encoded, 'receipt mismatch')
    print(json.dumps({'passed': True, 'word_count': len(words), 'source_checks': source_checks,
                      'literal_funnel_checks': child_checks,
                      'receipt_sha256': hashlib.sha256(encoded.encode()).hexdigest()}))


if __name__ == '__main__':
    main()
