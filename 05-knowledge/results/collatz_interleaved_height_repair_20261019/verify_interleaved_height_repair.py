#!/usr/bin/env python3
"""Exact all-lift offset controls for the bounded-interleaving theorem.

The finite universe is 12,798 endpoint offsets, not Collatz starting integers.
No imports from earlier checkers, no network, and no disabled assertions.
"""
from fractions import Fraction as F
from collections import Counter
import hashlib
import json

checks = 0

def require(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)

def branch(word):
    a, b, modulus = 1, 0, 1
    for bit in word:
        if bit == '1':
            a, b = 3*a, 3*b + modulus
        elif bit != '0':
            raise ValueError(bit)
        modulus *= 2
    return a, b, modulus

def cylinder(word):
    a, b, modulus = branch(word)
    return (-b * pow(a, -1, modulus)) % modulus

def step(n):
    return (3*n+1)//2 if n % 2 else n//2

def parity(n, length):
    word = ''
    for _ in range(length):
        word += str(n % 2)
        n = step(n)
    return word, n

def rational(q):
    return f'{q.numerator}/{q.denominator}'

K = '110110110110110101011011101010'
Ku = '001111111111111101010110000100'
Fword, Fu = '10', '01'
H, Hu = '110101000', '001101001'
V, Vu = 4*K + Fword, 4*Ku + Fu
P, C0, Q = branch(K)
C = F(branch(V)[1], branch(V)[0]-branch(V)[2])
E = C + 1
require(C == F(83829391903560839836419813590128166431,
               157489106280555890461927692448373459), 'maximal center')
require(F(C0, P-Q) < C, 'K center below V center')
require(12798 < 24*E < 12799, '24-episode integer offset bound')
require(13332 < 25*E < 13333, '25-episode comparison')
for name, word, upper in [('K', K, Ku), ('F', Fword, Fu),
                         ('H', H, Hu), ('V', V, Vu)]:
    a, b, modulus = branch(word)
    au, bu, mu = branch(upper)
    require((au, bu, mu) == (a, b+modulus-a, modulus), name+' companion map')
    require(cylinder(upper) == cylinder(word)+1, name+' companion guard')
    if name in ('K', 'V'):
        for j in range(1, len(word)+1):
            aj, bj, mj = branch(word[:j])
            require(aj > mj, name+' every nonempty prefix expands')
    else:
        require(a < modulus and a+b <= modulus, name+' reset above one')
        for j in range(len(word)):
            aj, bj, mj = branch(word[:j])
            require(F(aj, mj) > F(a, modulus), name+' proper-prefix slope')
            require(F(aj+bj, mj) >= F(a+b, modulus), name+' proper prefix at one')

r = cylinder(K)
r4 = cylinder(K*4)
mod4 = 2**120
require(r4 == 588943162011782273059357097228697595, 'forced K^4 residue')
require(parity(r4, 120)[0] == K*4, 'direct forced-prefix replay')
require(r4 > 466400 and 0 < r4 < mod4, 'source threshold and positive lift')
for ell in (1, 2, 3):
    a, b, modulus = branch(K*ell + Fword)
    require(a < modulus and a*r+b < modulus*r, 'short expansion already descends')
ah, bh, mh = branch(H)
af, bf, mf = branch(Fword)
require(F(ah, mh) < F(af, mf) and F(bh, mh) < F(bf, mf), 'H below F')

rows = []
classes = {}
maximum_threshold = F(0)
maximizers = []
for d in range(1, 12799):
    n = r4 + d
    a, b, modulus = 1, 0, 1
    word = ''
    for j in range(1, 121):
        bit = n % 2
        word += str(bit)
        if bit:
            a, b = 3*a, 3*b + modulus
        modulus *= 2
        n = step(n)
        if a < modulus:
            break
    else:
        raise RuntimeError('Uncertified offset '+str(d))
    require(j <= 92, '92-step continuation certificate')
    require(branch(word) == (a, b, modulus), 'separate branch reconstruction')
    require((r4+d-cylinder(word)) % modulus == 0, 'complete lawful-word guard')
    require(mod4 % modulus == 0, 'all inherited lifts share the parity prefix')
    threshold = F(a*d+b, modulus-a)
    require(threshold < 466400, 'exact original-source descent threshold')
    if threshold > maximum_threshold:
        maximum_threshold = threshold
        maximizers = [d]
    elif threshold == maximum_threshold:
        maximizers.append(d)
    # A positive control at a second lift; the congruence and inequality above
    # are the proof for every lift, rather than this literal replay.
    x = r4+17*mod4
    w2, y2 = parity(x+d, j)
    require(w2 == word and y2 == (a*(x+d)+b)//modulus and y2 < x,
            'literal lifted positive control')
    rows.append((d, j, a, b, modulus, rational(threshold)))
    key = (j, d % modulus)
    if key not in classes:
        classes[key] = [j, d % modulus, modulus, word, a, b, 0, d, d]
    classes[key][6] += 1
    classes[key][8] = d
require(len(rows) == 12798, 'complete offset universe')
require(len(classes) == 718, 'compressed offset classes')
require(sum(row[6] for row in classes.values()) == 12798, 'class coverage count')
for d in range(1, 12799):
    hits = [row for row in classes.values() if d % row[2] == row[1]]
    require(len(hits) == 1, 'disjoint exhaustive congruence classes')
    require(hits[0][7] <= d <= hits[0][8], 'class endpoint range')
require(max(row[1] for row in rows) == 92, 'maximum repair depth')
require([row[0] for row in rows if row[1] == 92] == [4], 'unique deepest repair')
require(maximum_threshold == F(810692053755261919, 1738366812781),
        'largest exact threshold')
require(maximizers == [11808], 'largest-threshold offset')

# Hostile control: the next carry budget includes an offset for which this
# particular 120-bit inherited prefix does not supply a contracting prefix.
hostile = 12852
hostile_word, _ = parity(r4+hostile, 120)
for j in range(1, 121):
    a, b, modulus = branch(hostile_word[:j])
    require(a >= modulus, '12852 has no contraction in the inherited prefix')
require(hostile <= 25*E, 'hostile lies in next offset budget')

# A guarded near-resonant family rules out an episode-independent carry bound.
lam, mu, c = F(P, Q), F(3, 4), F(C0, P-Q)
s = lam**11 * mu**3
e = c*s + 1 - (c+1)*mu**3
require(1 < s and mu*s**24 < 1, '24 excursions followed by first crossing')
large_carry = mu*e*sum(s**i for i in range(24)) + F(1, 4)
require(large_carry > 534, 'single-episode carry bound fails for interleaving')
w = (K*11 + Fword*3)*24 + Fword
a, b, modulus = branch(w)
require((len(w), w.count('1')) == (8066, 5089), 'hostile word length and odd count')
require(F(a, modulus) == mu*s**24 and F(b, modulus) == large_carry,
        'independent full-word carry reconstruction')
# For one complete cycle of expansion/reset blocks, each reset endpoint is
# at least its final expanding endpoint, and all later copies multiply it.
for j in (1, 2, 3):
    require(lam**11 * mu**j > 1, 'no earlier reset crossing in a copy')

certificate = json.dumps(rows, separators=(',', ':')).encode()
print(json.dumps({
    'status': 'PASS',
    'explicit_checks': checks,
    'scope': '12798 offset certificates valid on every inherited lawful lift; no source census or arbitrary entry',
    'bank_center': rational(C),
    'per_episode_carry_bound': rational(E),
    'episode_limit': 24,
    'max_positive_offset': 12798,
    'forced_source_residue': r4,
    'forced_source_modulus': mod4,
    'repair_max_steps': 92,
    'unique_deepest_offset': 4,
    'max_affine_threshold': rational(maximum_threshold),
    'max_threshold_offset': 11808,
    'compressed_class_count': len(classes),
    'offset_rows_sha256': hashlib.sha256(certificate).hexdigest(),
    'repair_time_counts': dict(sorted(Counter(row[1] for row in rows).items())),
    'unsupported_25_episode_offset': hostile,
    'hostile_word': {'length': len(w), 'odd_bits': w.count('1'),
                     'normalized_carry_above': 534,
                     'normalized_carry_decimal_for_display': float(large_carry),
                     'total_slope_decimal_for_display': float(mu*s**24)},
}, indent=2, sort_keys=True))
