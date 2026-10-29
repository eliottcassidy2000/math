"""Positive fixed-point laws, numerical-semigroup gaps and minimal arity order."""
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import json
import charged

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def convolution(a, b):
    out = Counter()
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] += x * y
    return out


def power(a, n):
    out = Counter({0: 1})
    for _ in range(n):
        out = convolution(out, a)
    return out


def distribution(p, k):
    """Integer histogram in the full affine wreath action, independently of cycles."""
    a = Counter({1: 1})
    h = (p - 1) // 2
    for _ in range(k):
        size = sum(a.values())
        out = power(a, p)
        for x, m in a.items():
            out[x] += p * (h - 1) * size ** (p - 1) * m
        out[0] += (p - 1) * size ** p
        need(sum(out.values()) == p * h * size ** p, 'integer offspring law total')
        a = out
    return a


def subset_sum(a, b):
    """Boolean convolution; exact support, no coefficient cancellation assumed."""
    mask = sum(1 << v for v in b)
    out = 0
    for x in a:
        out |= mask << x
    return {i for i in range(out.bit_length()) if out >> i & 1}


def sum_power(a, n):
    out = {0}
    for _ in range(n):
        out = subset_sum(out, a)
    return out


def characteristic(support):
    out = [1]
    for a in sorted(support):
        nxt = [0] * (len(out) + 1)
        for i, x in enumerate(out):
            nxt[i] -= a * x
            nxt[i + 1] += x
        out = nxt
    return out


def rank_mod(a, prime):
    a = [row[:] for row in a]
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col] % prime), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][col], -1, prime)
        a[r] = [x * inv % prime for x in a[r]]
        for i in range(r + 1, len(a)):
            c = a[i][col]
            if c:
                a[i] = [(x - c * y) % prime for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def check_moments(hist):
    support = sorted(hist)
    m = len(support)
    coeff = characteristic(support)
    # Unnormalized moments have the same minimal recurrence and Hankel rank.
    moments = [sum(c * a ** r for a, c in hist.items()) for r in range(2 * m + 5)]
    for r in range(m + 5):
        need(sum(coeff[i] * moments[r + i] for i in range(m + 1)) == 0,
             'characteristic annihilates every checked arity window')
    modulus = 1000000007
    matrix = [[moments[i + j] % modulus for j in range(m)] for i in range(m)]
    need(rank_mod(matrix, modulus) == m, 'full modular Hankel minor certifies rational rank')
    return {'minimal_order_from_r0': m, 'minimal_order_from_r1': m - int(0 in hist),
            'characteristic_coefficients_ascending': coeff,
            'hankel_certificate_prime': modulus}


def main(write=True):
    global CHECKS
    CHECKS = 0
    rows = []
    for p in (7, 11, 19):
        support = {1}
        for k in range(1, 5):
            support = sum_power(support, p) | support | {0}
            N = p ** (k - 1)
            predicted = {a + p * b for b in range(N + 1) for a in range(N - b + 1)} if N < 200 else None
            # Independent coin reachability, followed by reflection at p^k.
            reachable = [False] * (p ** k + 1)
            reachable[0] = True
            for d in range(1, len(reachable)):
                reachable[d] = (d >= p - 1 and reachable[d - p + 1]) or (d >= p and reachable[d - p])
            reflected = {p ** k - d for d, ok in enumerate(reachable) if ok}
            need(support == reflected, 'recursive support equals reflected two-coin semigroup')
            if predicted is not None:
                need(support == predicted, 'bottom-level independent offspring support')
            expected = 3 if k == 1 else p ** k + 1 - (p - 1) * (p - 2) // 2
            need(len(support) == expected, 'exact all-depth Hankel order formula')
            row = {'prime': p, 'depth': k, 'fixed_support_size': len(support),
                   'missing_deficits': [i for i, ok in enumerate(reachable) if not ok]}
            if p in (7, 11) and k == 2:
                hist = distribution(p, k)
                need(set(hist) == support, 'positive coefficient support matches Boolean recurrence')
                row.update(check_moments(hist))
            rows.append(row)
    x, y = {1}, {1}
    ternary = []
    for k in range(1, 7):
        y = subset_sum(y, {2 * a for a in x})
        x = sum_power(x, 3) | {0}
        need(x == set(range(0, 3 ** k + 1, 3)), 'ternary full support arithmetic progression')
        need(y == set(range(3, 3 ** k + 1, 6)), 'ternary mirrored support arithmetic progression')
        ternary.append({'depth': k, 'automorphism_order_in_arity': len(x),
                        'centralizer_order_in_arity': len(y)})
    c = charged.inventory(3, 3, True)
    hist = Counter()
    for (_, cycles), multiplicity in c.items():
        hist[cycles.count(1)] += multiplicity
    ternary_example = check_moments(hist)
    need(ternary_example['minimal_order_from_r0'] == 5, '27-vertex five-state arity recurrence')
    out = {'status': 'FINITE-EXACT support and Hankel certificates for PROVED formulas',
           'checks': CHECKS, 'paley_rows': rows, 'ternary_rows': ternary,
           'ternary_depth3_recurrence': ternary_example}
    if write:
        (HERE / 'memory.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'checks': CHECKS, 'paley_ranks': [(r['prime'], r['depth'], r['fixed_support_size']) for r in rows],
                      'ternary': ternary}))
    return out


if __name__ == '__main__':
    main()
