#!/usr/bin/env python3
"""Exact signed-orbital, Paley-code and eta-series controls; stdlib only."""
from itertools import combinations, permutations
from math import comb, gcd
from collections import Counter
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def induced_pairs(p):
    pairs = list(combinations(range(len(p)), 2))
    lookup = {v: i for i, v in enumerate(pairs)}
    return tuple(lookup[tuple(sorted((p[a], p[b])))] for a, b in pairs)


def orbital_orientations(n, action):
    """Each action element is (permutation, global reversal bit)."""
    pairs = list(combinations(range(n), 2))
    lookup = {p: i for i, p in enumerate(pairs)}
    unseen = set(range(len(pairs)))
    basis = []
    while unseen:
        seed = min(unseen)
        u, v = pairs[seed]
        values = {}
        for p, reversal in action:
            a, b = p[u], p[v]
            j = lookup[tuple(sorted((a, b)))]
            bit = 1 ^ reversal ^ int(a > b)
            if j in values and values[j] != bit:
                return None
            values[j] = bit
        unseen.difference_update(values)
        basis.append(values)
    return basis


def linear_dimension(n, action):
    """Independent GF(2) elimination of all signed pair constraints."""
    pairs = list(combinations(range(n), 2))
    lookup = {p: i for i, p in enumerate(pairs)}
    pivots = {}
    for p, reversal in action:
        for i, (u, v) in enumerate(pairs):
            a, b = p[u], p[v]
            j = lookup[tuple(sorted((a, b)))]
            row = (1 << i) ^ (1 << j)
            rhs = reversal ^ int(a > b)
            while row:
                k = row.bit_length() - 1
                if k not in pivots:
                    pivots[k] = (row, rhs)
                    break
                old, value = pivots[k]
                row ^= old
                rhs ^= value
            if not row and rhs:
                return None
    return len(pairs) - len(pivots)


def Hamilton_count(n, mask):
    pairs = list(combinations(range(n), 2))
    incoming = [0] * n
    for k, (a, b) in enumerate(pairs):
        if (mask >> k) & 1:
            incoming[b] |= 1 << a
        else:
            incoming[a] |= 1 << b
    table = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        table[1 << v][v] = 1
    for subset in range(1, 1 << n):
        for v in range(n):
            if not subset & (1 << v):
                continue
            previous = subset ^ (1 << v)
            if not previous:
                continue
            allowed = previous & incoming[v]
            while allowed:
                bit = allowed & -allowed
                table[subset][v] += table[previous][bit.bit_length() - 1]
                allowed ^= bit
    return sum(table[-1])


def masks_from_basis(basis):
    for choices in range(1 << len(basis)):
        mask = 0
        for j, values in enumerate(basis):
            for k, bit in values.items():
                mask |= (bit ^ ((choices >> j) & 1)) << k
        yield mask


def check_action(name, n, action, expected, paths=False):
    basis = orbital_orientations(n, action)
    dim = None if basis is None else len(basis)
    independent = linear_dimension(n, action)
    require(dim == independent == expected, name + ': incompatible dimensions')
    record = {'name': name, 'vertices': n, 'group_order': len(action),
              'free_bits': dim, 'orientations': 0 if dim is None else 2 ** dim}
    if paths:
        histogram = Counter(Hamilton_count(n, mask) for mask in masks_from_basis(basis))
        record['Hamilton_histogram'] = dict(sorted(histogram.items()))
    return record


def action_checks():
    rows = []
    for m in (3, 5, 7, 9, 11, 13):
        rotations = [(tuple((i + a) % m for i in range(m)), 0) for a in range(m)]
        signed = rotations + [(tuple((a - i) % m for i in range(m)), 1) for a in range(m)]
        rows.append(check_action('D%d natural signed' % (2*m), m, signed,
                                 (m-1)//2, paths=m <= 5))
        pair_rotations = [(induced_pairs(p), s) for p, s in rotations]
        rows.append(check_action('C%d pair lift' % m, comb(m, 2), pair_rotations,
                                 (m-1)*(m-2)*(m+1)//8, paths=m == 5))
        pair_signed = [(induced_pairs(p), s) for p, s in signed]
        rows.append(check_action('D%d pair lift signed' % (2*m), comb(m, 2), pair_signed,
                                 1 if m == 3 else None))
    a5 = []
    for p in permutations(range(5)):
        inversions = sum(p[i] > p[j] for i in range(5) for j in range(i+1, 5))
        if inversions % 2 == 0:
            a5.append((p, 0))
    rows.append(check_action('A5 natural', 5, a5, None))
    rows.append(check_action('A5 Petersen pair vertices', 10,
                             [(induced_pairs(p), s) for p, s in a5], None))
    c3xc3 = [(tuple([(i+a)%3 for i in range(3)] +
                    [3+(i+b)%3 for i in range(3)]), 0)
              for a in range(3) for b in range(3)]
    rows.append(check_action('C3xC3 two triples', 6, c3xc3, 3, paths=True))
    return rows


def euler_product(N, exponent, degree):
    product = [1] + [0] * degree
    for d in range(1, degree+1):
        if d % N == 0:
            continue
        factor = [(j*d, (-1)**j*comb(exponent, j))
                  for j in range(min(exponent, degree//d)+1)]
        updated = [0] * (degree+1)
        for k, a in enumerate(product):
            for offset, b in factor:
                if k + offset <= degree:
                    updated[k+offset] += a*b
        product = updated
    return product


def multiply(a, b, degree):
    return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(degree+1)]


def v2(n):
    return (n & -n).bit_length()-1


def eta_coefficients(N, degree):
    exponent = 24 // (N-1)
    sigma = [0] * (degree+1)
    for d in range(1, degree+1):
        for k in range(d, degree+1, d):
            sigma[k] += d
    recurrence = [1]
    for n in range(1, degree+1):
        numerator = -exponent * sum(
            (sigma[k] - (N*sigma[k//N] if k % N == 0 else 0))*recurrence[n-k]
            for k in range(1, n+1))
        require(numerator % n == 0, 'nonintegral eta coefficient')
        recurrence.append(numerator // n)
    product = euler_product(N, exponent, degree)
    require(product == recurrence, 'eta independent paths disagree')
    depth = v2(exponent)
    for s in range(depth+1):
        root = euler_product(N, exponent//(2**s), degree)
        for _ in range(s):
            root = multiply(root, root, degree)
        require(root == recurrence, 'formal dyadic root failed')
    require(recurrence[1] % (2**(depth+1)) != 0, 'next root not obstructed')
    first_odd = next(n for n in range(1, degree+1) if recurrence[n] % 2)
    require(first_odd == 2**depth, 'wrong Frobenius support depth')
    for n in range(1, degree+1):
        require(recurrence[n] % (2**max(0, depth-v2(n))) == 0,
                'dyadic coefficient lower bound failed')
    return {'level': N, 'exponent': exponent, 'checked_degree': degree,
            'max_integral_dyadic_root_depth': depth,
            'first_odd_coefficient_after_constant': first_odd,
            'first_13_coefficients': recurrence[:13]}


def paley_code(p):
    squares = {i*i % p for i in range(1, p)}
    rows = [sum(1 << j for j in range(p) if (j-i) % p in squares) for i in range(p)]
    pivots = {}
    for row in rows:
        while row:
            k = row.bit_length()-1
            if k not in pivots:
                pivots[k] = row
                break
            row ^= pivots[k]
    code = [0]
    for row in pivots.values():
        code += [word ^ row for word in code]
    weights = Counter(x.bit_count() for x in code)
    extension = Counter(w + w % 2 for x in code for w in [x.bit_count()])
    return {'order': p, 'dimension': len(pivots), 'weights': dict(sorted(weights.items())),
            'parity_extension_weights': dict(sorted(extension.items()))}


def main():
    levels = [d+1 for d in range(1, 25) if 24 % d == 0]
    supplied = [4, 2, 5, 3, 7, 13]
    even_levels = [N for N in levels if (24//(N-1)) % 2 == 0]
    require(set(supplied) == set(even_levels), 'literal set changed')
    codes = [paley_code(p) for p in (7, 23)]
    require(codes[1]['parity_extension_weights'] == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1},
            'Golay control failed')
    print(json.dumps({'signed_actions': action_checks(), 'literal_set': supplied,
                      'integer_eta_levels': levels, 'even_exponent_levels': even_levels,
                      'eta': [eta_coefficients(N, 128) for N in levels],
                      'Paley_codes': codes}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
