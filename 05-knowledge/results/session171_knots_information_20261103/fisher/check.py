#!/usr/bin/env python3
"""Exact finite controls for information in linear holonomy probes."""
import json
from fractions import Fraction
from itertools import combinations


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def parity(x):
    return x.bit_count() & 1


def span(rows):
    out = {0}
    for row in rows:
        out |= {x ^ row for x in tuple(out)}
    return out


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def rm_rows(d, k):
    monomials = [mask for mask in range(1 << d) if mask.bit_count() <= k]
    rows = [sum(1 << j for j, mask in enumerate(monomials)
                if y & mask == mask) for y in range(1 << d)]
    return monomials, rows


def fibres(K, rows):
    out = {}
    for x in range(1 << K):
        t = sum(parity(x & row) << j for j, row in enumerate(rows))
        out.setdefault(t, []).append(x)
    return out


def variance(values):
    mean = Fraction(sum(values), len(values))
    return sum((v - mean) ** 2 for v in values) / len(values)


def projected_variance(values, fs):
    mean = Fraction(sum(values), len(values))
    return sum(len(xs) * (Fraction(sum(values[x] for x in xs), len(xs)) - mean) ** 2
               for xs in fs.values()) / len(values)


def rat(x):
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'


def case(d, k, selected, name, all_characters=False):
    monomials, all_rows = rm_rows(d, k)
    K = len(monomials)
    rows = [all_rows[y] for y in selected]
    a = rank(rows)
    fs = fibres(K, rows)
    rowspace = span(rows)
    require(len(fs) == 1 << a, 'image size')
    require({len(xs) for xs in fs.values()} == {1 << (K-a)}, 'equal affine fibre sizes')
    require(len(rowspace) == 1 << a, 'rowspace rank')
    probes = range(1, 1 << K) if all_characters else sorted({1, (1 << K)-1, *all_rows})
    for ell in probes:
        vals = [1 - 2 * parity(ell & x) for x in range(1 << K)]
        require(variance(vals) == 1, 'character full variance')
        require(projected_variance(vals, fs) == (1 if ell in rowspace else 0), 'character Fisher classification')
    vals = [int(x == 0) for x in range(1 << K)]
    full = variance(vals)
    observed = projected_variance(vals, fs)
    require(full == Fraction(1, 1 << K) - Fraction(1, 1 << (2*K)), 'identity full Fisher')
    require(observed == Fraction(1, 1 << (2*K-a)) - Fraction(1, 1 << (2*K)), 'identity observed Fisher')
    # Fourier expansion with known coefficients, including a constant term.
    terms = {ell: ((ell * 7) % 9) - 4 for ell in range(min(1 << K, 31))}
    vals = [sum(c * (1 - 2 * parity(ell & x)) for ell, c in terms.items()) for x in range(1 << K)]
    expected_full = sum(c*c for ell, c in terms.items() if ell)
    expected_observed = sum(c*c for ell, c in terms.items() if ell and ell in rowspace)
    require(variance(vals) == expected_full, 'Parseval full Fisher')
    require(projected_variance(vals, fs) == expected_observed, 'Fourier projected Fisher')
    # Direct S170 degree and minimum distance control, for these small codes.
    weights = [sum(parity(x & row) for row in all_rows) for x in range(1, 1 << K)]
    require(min(weights) == 1 << (d-k), 'RM minimum distance')
    return dict(name=name, d=d, degree=k, code_dimension=K, readout_count=len(rows), rank=a,
                residual_entropy_bits=K-a, fibre_count=len(fs), fibre_size=1 << (K-a),
                minimum_distance=min(weights), characters_checked=len(probes),
                identity_fisher_full=rat(full), identity_fisher_observed=rat(observed),
                identity_fisher_lost=rat(full-observed))


def adaptive_control():
    # First read coefficient 0; only on the 1-branch also read coefficient 1.
    K = 3
    fs = {}
    for x in range(1 << K):
        t = (0,) if not (x & 1) else (1, (x >> 1) & 1)
        fs.setdefault(t, []).append(x)
    require(sorted(map(len, fs.values())) == [2, 2, 4], 'adaptive fibres')
    entropy = sum(Fraction(len(xs), 1 << K) * (len(xs).bit_length()-1) for xs in fs.values())
    expected_rank = sum(Fraction(len(xs), 1 << K) * len(t) for t, xs in fs.items())
    require(entropy == Fraction(3, 2) and expected_rank == Fraction(3, 2), 'adaptive entropy')
    require(entropy == K - expected_rank, 'adaptive rank identity')
    return dict(dimension=K, leaf_sizes=[2, 2, 4], residual_entropy_bits=rat(entropy), expected_rank=rat(expected_rank))


def main():
    cases = [case(3, 1, [0], 'one sheet in level-four balanced family', True),
             case(3, 1, [0, 0, 0], 'duplicate reads carry no new noiseless information', True),
             case(3, 1, [0, 1, 2, 4], 'optimal level-four bank', True),
             case(4, 2, [0, 1, 2, 4, 8], 'degree-one bank misses degree-two coefficients'),
             case(4, 2, [y for y in range(16) if y.bit_count() <= 2], 'optimal level-five bank'),
             case(4, 2, list(range(16)), 'redundant full level-five cube')]
    print(json.dumps(dict(status='PASS', cases=cases, adaptive=adaptive_control()), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
