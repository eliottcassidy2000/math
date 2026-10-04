"""Exact affine-axis carriers; finite curvature, word, and compiler controls.

No floating point chooses a sign. Checks remain active under python -O.
The infinite elementary proofs and inherited boundaries are in the note.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import gcd
from pathlib import Path
import importlib.util


def need(test, message):
    if not test:
        raise ValueError(message)


def data(word):
    need(bool(word) and all(type(a) is int and a >= 1 for a in word),
         'nonempty positive integer valuation word')
    p, q, carry = 1, 1, 0
    for a in word:
        p, q, carry = 3*p, q << a, 3*carry+q
    return p, q, carry


def independent_data(word):
    total = 0
    carry = 0
    for i, a in enumerate(word):
        carry += 3**(len(word)-1-i) * 2**total
        total += a
    return 3**len(word), 2**total, carry


def axis(word):
    p, q, carry = data(word)
    return F(carry, q-p)


def affine(word, point):
    need(type(point) is int or isinstance(point, F), 'exact rational point')
    p, q, carry = data(word)
    return F(p*point+carry, q)


def v2(n):
    need(type(n) is int and n != 0, 'nonzero integer valuation')
    n = abs(n)
    return (n & -n).bit_length()-1


def decode(p, q, carry):
    need(all(type(n) is int for n in (p, q, carry)), 'integer affine data')
    length, original_p = 0, p
    while p > 1 and p % 3 == 0:
        p //= 3
        length += 1
    need(p == 1 and length >= 1, 'positive power of three')
    need(q >= 2 and q & (q-1) == 0, 'positive power of two')
    total, original_carry = q.bit_length()-1, carry
    output = []
    for remaining in range(length, 1, -1):
        delta = carry-3**(remaining-1)
        need(delta > 0, 'positive remaining carry')
        a = v2(delta)
        need(a >= 1, 'positive next valuation')
        output.append(a)
        total -= a
        carry = delta >> a
    need(carry == 1 and total >= 1, 'valid last valuation and carry')
    output.append(total)
    word = tuple(output)
    need(data(word) == (original_p, q, original_carry), 'decoded map')
    return word


def compositions(total):
    if total == 0:
        yield ()
    else:
        for first in range(1, total+1):
            for tail in compositions(total-first):
                yield (first,)+tail


def primitive_root(word):
    for length in range(1, len(word)+1):
        if len(word) % length == 0 and word == word[:length]*(len(word)//length):
            return word[:length]
    raise ValueError('missing primitive root')


def replay_exact(source, word):
    current = source
    for a in word:
        numerator = 3*current+1
        need(v2(numerator) == a, 'literal exact valuation')
        current = numerator >> a
    return current


def earliest_strict_score_error(c1, c2, cap=5000):
    need(type(c1) is int and type(c2) is int and c1 > 0 > c2,
         'one-letter-correct integer score')
    ternary = 1
    for length in range(1, cap+1):
        ternary *= 3
        last_growing_twos = ternary.bit_length()-1-length
        # Drift is strictly decreasing with the number of 2-letters. These
        # two neighboring counts decide whether the two thresholds disagree.
        for twos in (last_growing_twos, last_growing_twos+1):
            if 0 <= twos <= length:
                score = c1*(length-twos)+c2*twos
                gap = ternary-(1 << (length+twos))
                if score*gap < 0:
                    return length, twos, score, gap
    raise ValueError('declared score-search cap reached')


def run():
    words = [w for total in range(1, 13) for w in compositions(total)]
    by_axis = defaultdict(list)
    replay_count = 0
    for word in words:
        p, q, carry = data(word)
        need((p, q, carry) == independent_data(word), 'independent carry sum')
        need(decode(p, q, carry) == word, 'lossless word decoder')
        a = axis(word)
        need(affine(word, a) == a, 'finite axis endpoint')
        need(F((p+q)**2, p*q) == 4+F((p-q)**2, p*q) > 4,
             'strict hyperbolic invariant')
        rotated = word[1:]+word[:1]
        need(affine(word[:1], a) == axis(rotated), 'marked cyclic rotation')
        by_axis[a].append(word)
        residue = (q-carry)*pow(p, -1, 2*q) % (2*q)
        need(residue > 0 and residue % 2 == 1, 'positive exact cylinder')
        for lift in (0, 1, 7):
            source = residue+2*q*lift
            end = replay_exact(source, word)
            need(F(end) == affine(word, F(source)), 'affine/literal endpoint')
            need(F(end-source) == (F(p, q)-1)*(F(source)-a),
                 'marked-source displacement')
            replay_count += 1
    for same_axis_words in by_axis.values():
        need(len({primitive_root(w) for w in same_axis_words}) == 1,
             'equal axes have one marked primitive root')

    small = [w for w in words if sum(w) <= 8]
    order_checks = 0
    for u, v in product(small, repeat=2):
        pu, qu, bu = data(u)
        pv, qv, bv = data(v)
        slope_u, slope_v = F(pu, qu), F(pv, qv)
        direct_difference = F(pv*bu+bv*qu-pu*bv-bu*qv, qu*qv)
        predicted = (slope_u-1)*(slope_v-1)*(axis(v)-axis(u))
        need(direct_difference == predicted, 'axis order defect')
        need((direct_difference == 0) == (u+v == v+u)
             == (primitive_root(u) == primitive_root(v)), 'centralizer criterion')
        order_checks += 1

    resonant = [w for w in words if len(w) == 7 and sum(w) == 11]
    integer_axes = sorted(int(axis(w)) for w in resonant if axis(w).denominator == 1)
    need(len(resonant) == 210 and len({axis(w) for w in resonant}) == 210,
         'same-spectrum distinct marked axes')
    need(integer_axes == [-91, -61, -55, -41, -37, -25, -17], 'inherited signed cycle')
    need(axis((1,2)) == -5 and axis((2,1)) == -7, 'rotation hostile')
    need(data((1,2))[:2] == data((2,1))[:2], 'same clock hostile')
    need(replay_exact(11, (1,2)) == 13, 'flat-score word actually grows')
    need(earliest_strict_score_error(1, -1) == (7, 4, -1, 139),
         'first strict equal-weight curvature error')
    shadow_word = (1,1,1,2,1,1,4)
    need(data(shadow_word) == (2187, 2048, 2363) and axis(shadow_word) == -17,
         'inherited minus17 anchor')
    need(sum(3-2*a for a in shadow_word) == -1, 'extended rational curvature score')
    need(all(data(shadow_word[:i])[0] > data(shadow_word[:i])[1]
             for i in range(1, len(shadow_word)+1)), 'all shadow prefixes grow')
    for m in range(1, 25):
        source = 2*2048**m-17
        current = source
        for a in shadow_word*m:
            numerator = 3*current+1
            need(v2(numerator) == a, 'independent exact minus17 shadow valuation')
            current = numerator >> a
            need(current > source, 'every shadow prefix exceeds immutable source')
        need(current == 2*2187**m-17, 'minus17 shadow endpoint')
    score_errors = []
    for c1 in range(1, 13):
        for c2 in range(-12, 0):
            g = gcd(c1, -c2)
            n1, n2 = -c2//g, c1//g
            need(c1*n1+c2*n2 == 0, 'formal curvature-flat word')
            need(3**(n1+n2) != 2**(n1+2*n2), 'flat is not neutral drift')
            score_errors.append((c1, c2, earliest_strict_score_error(c1, c2)))

    # Import the incoming compiler without running its larger census. This
    # checks an interface, not a second implementation of its all-length proof.
    path = Path(__file__).with_name('collatz_boundary_compiler_20261004.py')
    spec = importlib.util.spec_from_file_location('boundary_compiler', path)
    compiler = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(compiler)
    compiler_controls = 0
    compiler_status = defaultdict(int)
    for word in compiler.growth_words(6):
        for repetition in (1, 2, 3):
            repeated = word*repetition
            need(axis(repeated) == axis(word), 'repetition retains axis')
            row = compiler.compile_growth(word, repetition)
            coarse = tuple(row['coarse_word'])
            need(axis(repeated) < 0 < axis(coarse), 'guard changes axis side')
            pc, qc, bc = data(coarse)
            need(F(bc, qc*(qc-pc)) <= F(1121, 3328), 'inherited sharp carry bound')
            n = row['residue']+row['modulus']
            need(n > axis(coarse) and affine(coarse, F(n)) < n,
                 'positive lift lies beyond contracting axis')
            compiler_status[row['boundary_status']] += 1
            compiler_controls += 1
    hostile = (4,1,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3)
    need(replay_exact(165, hostile) == 167 and axis(hostile) > 165,
         'contracting word can grow below its own axis')
    need(affine((1,), F(27)) == 41 and affine((2,), F(41)) == 31 > 27,
         'moving checkpoint hostile')

    print('GEOMETRY / COLLATZ DRIFT CARRIERS: PROVED algebra; FINITE-EXACT controls')
    print('All positive words with total valuation<=12:', len(words))
    print('Distinct marked axes:', len(by_axis), '; exact integer replays:', replay_count)
    print('All ordered word pairs with both totals<=8:', order_checks)
    print('Same-clock j7,A11 words / distinct axes / integer axes:',
          len(resonant), len({axis(w) for w in resonant}), integer_axes)
    print('Letters1,2: axes -1,+1; words12,21: axes -5,-7; order defect -1/4')
    print('Equal-weight 5/7 formal curvature: flat word12 grows11->13; first strict error',
          earliest_strict_score_error(1, -1))
    print('Inherited minus17 shadow controls m1..24: score -m, all7m prefixes grow; m1:4079 ->4357')
    print('144 rational-score samples; latest first strict error:',
          max(score_errors, key=lambda item: item[2][0]))
    print('Incoming compiler growing words length<=6, repetitions1..3:',
          compiler_controls, dict(sorted(compiler_status.items())))
    print('Contracting-axis hostile source165 ->167; moving checkpoint27 ->41 ->31')
    print('ALL CHECKS PASSED; no new coverage or geometric Collatz conjugacy claimed')


if __name__ == '__main__':
    run()
