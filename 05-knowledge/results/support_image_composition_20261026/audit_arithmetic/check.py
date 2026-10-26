#!/usr/bin/env python3
"""Independent literal audit of residue clocks; no producer imports."""
import itertools as it
import json
import math
from collections import Counter, defaultdict
from functools import lru_cache

GATES = 0


def check(condition, reason):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(reason)


def shortcut(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def native(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def trace(n, modulus, length):
    out = [n % modulus]
    for _ in range(length):
        out.append(n % 2)
        n = shortcut(n)
        out.append(n % modulus)
    return tuple(out)


@lru_cache(None)
def arrows(modulus, r):
    return tuple((e, s) for e in (0, 1) for s in range(modulus)
                 if (2 * s - 3 ** e * r - e) % modulus == 0)


@lru_cache(None)
def paths(modulus, r, length):
    result = [(r,)]
    for _ in range(length):
        result = [p + (e, s) for p in result for e, s in arrows(modulus, p[-1])]
    return tuple(result)


def project(path, modulus):
    return tuple(value if i % 2 else value % modulus for i, value in enumerate(path))


def affine(word):
    a, b, denominator = 1, 0, 1
    for e in word:
        a, b = 3 ** e * a, 3 ** e * b + e * denominator
        denominator *= 2
    return a, b, denominator


def multiply(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
                 for i in range(2))


def apply(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(2)) for i in range(2))


def determinant(a):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))


def guarded_wheel(n, bound):
    return all(x % p != 0 or x == p for p in range(2, bound + 1) if prime(p)
               for x in (n, n + 2))


def main():
    cylinders = 0
    for modulus in range(1, 13):
        for r in range(modulus):
            outgoing = arrows(modulus, r)
            check(len(outgoing) == 2, 'outgoing labelled degree')
            literal = {(n % 2, shortcut(n) % modulus)
                       for n in range(1, 2 * modulus + 1) if n % modulus == r}
            check(set(outgoing) == literal, 'one-step positive existential image')
            if modulus % 2 == 0:
                check({e for e, _ in outgoing} == {r % 2}, 'even-modulus label')
                check(abs(outgoing[0][1] - outgoing[1][1]) == modulus // 2,
                      'even-modulus dyadic choice')
            else:
                check({e for e, _ in outgoing} == {0, 1}, 'odd-modulus labels')
                check((outgoing[0][1] == outgoing[1][1]) == (r == (modulus - 1) // 2),
                      'odd-modulus coincident successor')

        for length in range(5):
            q = modulus * 2 ** length
            formal = {p for r in range(modulus) for p in paths(modulus, r, length)}
            raw = defaultdict(list)
            for a in range(q):
                raw[trace(a, modulus, length)].append(a)
            check(set(raw) == formal and len(formal) == q, 'complete finite path image')
            check(all(len(a) == 1 for a in raw.values()), 'one source residue per path')
            positive = defaultdict(list)
            for n in range(1, 2 * q + 2):
                positive[trace(n, modulus, length)].append(n)
            for p, aa in raw.items():
                cylinders += 1
                a = aa[0]
                for lift in (-2, -1, 0, 1, 2):
                    check(trace(a + lift * q, modulus, length) == p, 'full integer cylinder')
                n = a if a else q
                check(trace(n, modulus, length) == p and n > 0, 'positive witness')
                heights = {0, 1, 5, q - 1, q, q + 1, 2 * q, 2 * q + 1,
                           max(0, a - 1), a, a + 1}
                for height in heights:
                    # Height 5 can exceed the enumerated interval only for q=1.
                    literal_count = sum(trace(k, modulus, length) == p
                                        for k in range(1, height + 1))
                    floor_count = (height - a) // q - (-a) // q
                    check(literal_count == floor_count, 'exact positive interval count')
                    if q >= height:
                        check(literal_count <= 1, 'single bounded candidate')
                word = p[1::2]
                multiplier, carry, denominator = affine(word)
                value = n
                native_value = n
                elapsed = 0
                for e in word:
                    check(value % 2 == e, 'actual parity guard')
                    value = shortcut(value)
                    for _ in range(1 + e):
                        native_value = native(native_value)
                        elapsed += 1
                    check(native_value == value, 'native clock correspondence')
                check(multiplier * n + carry == denominator * value, 'affine evaluator')
                check(elapsed == length + sum(word), 'ordinary elapsed time')

    covering, noncovering = [], []
    for fine in range(1, 25):
        for coarse in range(1, fine + 1):
            if fine % coarse:
                continue
            odd = (fine // coarse) % 2 == 1
            one_step_is_cover = True
            for b in range(fine):
                for length in range(1, 4):
                    counts = Counter(project(p, coarse) for p in paths(fine, b, length))
                    coarse_paths = paths(coarse, b % coarse, length)
                    check(set(counts) <= set(coarse_paths), 'reduction preserves labelled paths')
                    if odd:
                        check(all(counts[p] == 1 for p in coarse_paths), 'odd-ratio unique path lift')
                    if length == 1:
                        fibre_sizes = sorted(counts[p] for p in coarse_paths)
                        one_step_is_cover &= fibre_sizes == [1, 1]
                        check(fibre_sizes == ([1, 1] if odd else [0, 2]),
                              'sharp one-arrow covering boundary')
            check(one_step_is_cover == odd, 'covering iff odd quotient')
            (covering if odd else noncovering).append((coarse, fine))

    for modulus in range(1, 13):
        for length in range(1, 5):
            q = modulus * 2 ** length
            for a in range(q):
                check(shortcut(a) % (q // 2) == shortcut(a + q) % (q // 2),
                      'shrinking-precision deterministic map')
            zero = (0,) + (0, 0) * length
            minus = ((modulus - 1) % modulus,) + (1, (modulus - 1) % modulus) * length
            check(trace(q, modulus, length) == zero, 'zero-loop least positive source')
            check(trace(q - 1, modulus, length) == minus, 'one-loop least positive source')
            check(all(trace(n, modulus, length) != zero for n in range(1, q)),
                  'zero-loop positive minimality')
            check(all(trace(n, modulus, length) != minus for n in range(1, q - 1)),
                  'one-loop positive minimality')
    check(shortcut(0) == 0 and shortcut(-1) == -1, 'nonpositive fixed-point boundary')
    check(trace(10, 7, 1) == (3, 0, 5) and trace(3, 7, 1) == (3, 1, 5),
          'parallel arrows retain native-time label')
    check(not any(trace(n, 3, 1) == (0, 0, 0) for n in range(1, 6)),
          'finite height hostile')
    for length in range(5):
        for word in it.product((0, 1), repeat=length):
            a, b, q = affine(word)
            for cut in range(length + 1):
                au, bu, qu = affine(word[:cut])
                av, bv, qv = affine(word[cut:])
                check((a, b, q) == (av * au, av * bu + qu * bv, qu * qv),
                      'ordered affine word composition')

    matrices = {'U': ((2, -1), (1, 0)), 'A': ((2, 1), (1, 0)), 'D': ((1, 2), (0, 1))}
    for modulus in range(1, 25):
        word = 'A' + 'D' * (modulus - 1) + 'UD'
        matrix = ((1, 0), (0, 1))
        vector = (2, 1)
        for letter in word:
            matrix = multiply(matrices[letter], matrix)
            vector = apply(matrices[letter], vector)
        check(vector == (16 * modulus + 2, 4 * modulus + 1), 'explicit Berggren family')
        check(tuple(v % modulus for v in vector) == (2 % modulus, 1 % modulus),
              'root-vector modular collision')
        check(determinant(matrix) == -1 and apply(matrix, (2, 1)) == vector,
              'retained full-matrix determinant')
        check(vector[0] > vector[1] > 0 and math.gcd(*vector) == 1 and sum(vector) % 2 == 1,
              'actual primitive Euclid guard')
        if modulus % 2 == 1 and modulus > 1:
            check(determinant(matrix) % modulus != 1, 'opposite determinants remain distinct')

    for bound, height in ((3, 22), (5, 46), (7, 118)):
        q = next(x for x in range(bound + 1, 2 * bound + 3) if prime(x))
        check(height + 2 < q * q, 'height-complete sieve hypothesis')
        wheel = math.prod(p for p in range(2, bound + 1) if prime(p))
        by_trace = defaultdict(list)
        for n in range(3, height + 1):
            check(guarded_wheel(n, bound) == (prime(n) and prime(n + 2)),
                  'corrected wheel / actual twin predicate')
            if guarded_wheel(n, bound):
                by_trace[trace(n, wheel, 3)].append(n)
        brute = [n for n in range(3, height + 1) if prime(n) and prime(n + 2)]
        check(sorted(x for values in by_trace.values() for x in values) == brute,
              'sieve and path certificate compose without losing sources')
    check(3 % 6 == 9 % 6 and guarded_wheel(3, 3) and not guarded_wheel(9, 3),
          'prime-equality exception is not residue invariant')
    print(json.dumps({'status': 'PASS', 'checks': GATES,
                      'clock_moduli': [1, 12], 'path_lengths': [0, 4],
                      'literal_source_cylinders': cylinders,
                      'refinement_fine_maximum': 24,
                      'covering_pairs': len(covering), 'noncovering_pairs': len(noncovering),
                      'cardinality_and_floors_use_positive_integer_sources': True,
                      'producer_imports': False}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
