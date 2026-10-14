#!/usr/bin/env python3
"""Exact Pascal-word relative maxima and golden carry controls.

The bisection brackets are rational and verified by a max-plus word DAG.
Only the displayed normal-limit approximation uses floating point.
"""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path


def require(value, message):
    if not value:
        raise RuntimeError(message)


def join(a, b):
    la, sa, ma = a
    lb, sb, mb = b
    return la + lb, sa + sb, max(ma, sa + mb)


def row_support(k, slope_numerator, denominator):
    plus = (1, denominator - slope_numerator, max(0, denominator - slope_numerator))
    minus = (1, -denominator - slope_numerator, 0)
    row = [plus, minus]
    for _ in range(1, k):
        row = [plus] + [join(row[j-1], row[j]) for j in range(1, len(row))] + [minus]
    out = (0, 0, 0)
    for node in row:
        out = join(out, node)
    return out[2] - slope_numerator * (1 << k)


def max_ratio_bracket(k, bits=30):
    # lambda=max e(n)/n=2 A_k. Integer arithmetic at every DAG node.
    denominator = 1 << bits
    low, high = 0, denominator
    while high - low > 1:
        mid = (low + high) // 2
        if row_support(k, mid, denominator) > 0:
            low = mid
        else:
            high = mid
    require(row_support(k, low, denominator) > 0, 'strict lower support')
    require(row_support(k, high, denominator) <= 0, 'upper support')
    return {'k': k, 'A_lower': str(Fraction(low, 2*denominator)),
            'A_upper': str(Fraction(high, 2*denominator)),
            'scaled_lower_numeric': low/(2*denominator)*math.sqrt(k),
            'scaled_upper_numeric': high/(2*denominator)*math.sqrt(k)}


def golden_value(word):
    a, b = 0, 0
    for bit in word:
        a, b = b + int(bit), a + b
    return a, b


def normalize(word):
    word = '0' + word
    while '11' in word:
        j = word.index('11')
        require(j > 0 and word[j-1] == '0', 'leftmost golden carry context')
        word = '0' + word[:j-1] + '100' + word[j+2:]
    return word.lstrip('0') or '0'


def golden_audit():
    canonical = {}
    count = 0
    for length in range(1, 15):
        for bits in itertools.product('01', repeat=length):
            word = ''.join(bits)
            value = golden_value(word)
            norm = normalize(word)
            require(value == golden_value(norm), 'exact golden carry')
            require('11' not in norm, 'golden normal form')
            require(value not in canonical or canonical[value] == norm, 'finite canonical consistency')
            canonical[value] = norm
            count += 1
    require(golden_value('11') == golden_value('100') == (1, 1), '11=100')
    require(golden_value('111') == golden_value('1001') == (2, 2), 'overlapping carry hostile')
    return {'binary_words_through_length': 14, 'words': count,
            'exact_values': len(canonical), 'examples': {w: normalize(w) for w in ('11', '111', '1111', '1011', '11011')}}


def small_conway_control():
    kmax = 17
    a = [0, 1, 1]
    for n in range(3, (1 << (kmax + 1)) + 1):
        j = a[-1]
        a.append(a[j] + a[n-j])
    word_row = ['+', '-']
    bridge_checks = 0
    for k in range(1, kmax + 1):
        word = ''.join(word_row)
        require(all((2*(a[n]-a[n-1])-1 == (1 if bit == '+' else -1))
                    for n, bit in enumerate(word, (1 << k)+1)), 'classical Pascal encoding finite control')
        for j, cell in enumerate(word_row):
            b = len(cell)
            total = cell.count('+') - cell.count('-')
            require(total*k == (k-2*j)*b, 'Pascal signed total')
            prefix, deviation = 0, 0
            for i, bit in enumerate(cell, 1):
                prefix += 1 if bit == '+' else -1
                centered = b*prefix-total*i
                require(centered >= 0, 'nonnegative Pascal bridge')
                deviation = max(deviation, centered)
            if 0 < j < k:
                cut = math.comb(k-1, j-1)
                left_total = cell[:cut].count('+')-cell[:cut].count('-')
                require(Fraction(b*left_total-total*cut, b) == Fraction(2*j*(k-j)*b, k*k*(k-1)), 'bridge split height')
            if k >= 8 and 3*k <= 8*j <= 5*k:
                bound = Fraction(8, 3*k)+2*Fraction(3, 4)**(k//8)
                require(Fraction(deviation, b*b) <= bound, 'uniform within-cell bridge bound')
            bridge_checks += 1
        if k in (4, 10, 17):
            maximum = max(Fraction(2*a[n]-n, n) for n in range(1 << k, (1 << (k+1))+1))
            for delta in (-1, 0, 1):
                slope = maximum + Fraction(delta, 10**8)
                sign = row_support(k, slope.numerator, slope.denominator)
                require((sign > 0) == (delta < 0) and (sign == 0) == (delta == 0), 'DAG support against literal sequence')
        word_row = ['+'] + [word_row[j-1]+word_row[j] for j in range(1, len(word_row))] + ['-']
    return {'sequence_through': 1 << (kmax + 1), 'literal_support_levels': [4, 10, 17],
            'literal_Pascal_bridge_checks': bridge_checks}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-level', type=int, default=1024)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    levels = [k for k in (10, 20, 29, 64, 128, 256, 400, 512, 1024) if k <= args.max_level]
    rows = []
    for k in levels:
        rows.append(max_ratio_bracket(k))
        print(json.dumps(rows[-1]), flush=True)
    density = lambda z: math.exp(-z*z/2)/math.sqrt(2*math.pi)
    distribution = lambda z: (1+math.erf(z/math.sqrt(2)))/2
    low, high = -1., 0.
    for _ in range(70):
        z = (low+high)/2
        if z*(1+distribution(z))+density(z) < 0:
            low = z
        else:
            high = z
    z = (low+high)/2
    result = {'status': 'FINITE-EXACT rational Pascal maxima; NUMERICAL evaluation of proved limiting constant',
              'rows': rows, 'normal_limit': {'z_numeric': z, 'constant_numeric': -z/2,
                                            'peak_position_numeric': 1+distribution(z)},
              'golden_carries': golden_audit(), 'independent_sequence': small_conway_control()}
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
