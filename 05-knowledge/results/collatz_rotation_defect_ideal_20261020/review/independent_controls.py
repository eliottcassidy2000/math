#!/usr/bin/env python3
"""Small, independent exact controls; imports no producer/repository code."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd
import json

counts = Counter()
digest = sha256()

def check(value, label):
    counts[label] += 1
    if not value:
        raise RuntimeError(label)

def carry(word, a):
    # Positional definition, deliberately not the producer's affine recurrence.
    indices = [i for i, bit in enumerate(word) if bit]
    return sum(2 ** i * a ** (len(indices) - k - 1)
               for k, i in enumerate(indices))

def rotations(word):
    return sorted({word[i:] + word[:i] for i in range(len(word))})

def common_context_options(u, v):
    # Includes unstripped, maximal, and every intermediate allowed context.
    left = 0
    while u[left] == v[left]:
        left += 1
    right = 0
    while u[len(u)-1-right] == v[len(v)-1-right]:
        right += 1
    return [(i, j) for i in range(left+1) for j in range(right+1)]

def decode(B, n, p, a):
    bits = []
    for _ in range(n):
        bit = B % 2
        bits.append(bit)
        B = (B - (a ** (p-1) if bit else 0)) // 2
        p -= bit
    return tuple(bits), B, p

def inspect(word, a):
    n = len(word)
    p = sum(word)
    D = 2**n - a**p
    B = carry(word, a)
    g = gcd(D, B)
    check(gcd(D, 2*a) == 1, 'clock_units')
    check(decode(B, n, p, a) == (word, 0, 0), 'carry_decoder')
    rs = rotations(word)
    bs = [carry(r, a) for r in rs]
    check(len(set(bs)) == len(rs), 'rotation_injectivity')

    # A mixed tree: vertex i attaches to floor((i-1)/2), neither star nor path
    # once there are sufficiently many rotations. Each edge uses a different
    # allowed amount of outside-context stripping.
    full_content = local_content = 0
    for i in range(1, len(rs)):
        j = (i-1)//2
        u, v = rs[i], rs[j]
        diff = bs[i] - bs[j]
        full_content = gcd(full_content, diff)
        options = common_context_options(u, v)
        selected = options[(i + a + sum(word)) % len(options)]
        for left, right in options:
            end = n-right
            core_u, core_v = u[left:end], v[left:end]
            ell, h = len(core_u), sum(core_u)
            delta = carry(core_u, a) - carry(core_v, a)
            check(sum(core_v) == h, 'equal_core_shape')
            check(diff == 2**left * a**sum(u[end:]) * delta,
                  'every_allowed_context')
            bound = (2**(ell-h)-1) * (a**h-2**h)//(a-2)
            check(0 < abs(delta) <= bound, 'short_core_bound')
            if (left, right) == selected:
                local_content = gcd(local_content, delta)
    check(full_content == local_content == g, 'mixed_tree_content')

    # Fraction reduction and actual parity-directed replay independently check
    # the integer gate, sign and the claimed itinerary for odd additive b.
    for b in (-5, -3, -1, 1, 3, 5):
        start = x = Fraction(b*B, D)
        nodes = []
        for bit in word:
            nodes.append(x)
            check(x.denominator % 2 == 1 and x.numerator % 2 == bit,
                  'odd_denominator_parity')
            x = (a*x+b)/2 if bit else x/2
        check(x == start, 'exact_rational_closure')
        integral = all(x.denominator == 1 for x in nodes)
        check(integral == (b % (abs(D)//g) == 0), 'integer_cycle_gate')
        check(all(x > 0 for x in nodes) == (b*D > 0), 'positive_sign_gate')
        check(len(set(nodes)) == len(rs), 'least_orbit_period')
        scaled_full_content = 0
        for i in range(1, len(bs)):
            scaled_full_content = gcd(scaled_full_content, b*(bs[i]-bs[0]))
        check(scaled_full_content == abs(b)*g, 'scaled_affine_carry_content')
    digest.update(repr((a, word, D, B, g, full_content, local_content)).encode())
    counts['word_base_cases'] += 1

shape_cases = 0
for a in (3, 5, 9):
    for n in range(2, 10):
        shape_values = {}
        for word in product((0, 1), repeat=n):
            p = sum(word)
            if not 0 < p < n:
                continue
            inspect(word, a)
            shape_values.setdefault(p, []).append((carry(word, a), word))
        for p, values in shape_values.items():
            check(len({B for B, _ in values}) == len(values), 'shape_injectivity')
            low = (a**p-2**p)//(a-2)
            high = 2**(n-p)*low
            check(min(values) == (low, (1,)*p+(0,)*(n-p)), 'unique_minimum')
            check(max(values) == (high, (0,)*(n-p)+(1,)*p), 'unique_maximum')
            shape_cases += 1

signed_word = tuple(map(int, '11110111000'))
inspect(signed_word, 3)
check(Fraction(carry(signed_word, 3), 2**11-3**7) == -17,
      'negative_seventeen_control')

powers = 0
for source in ('01', '000011', '0011011', '11110111000'):
    w = tuple(map(int, source))
    for a in (3, 5, 9):
        D = 2**len(w)-a**sum(w)
        B = carry(w, a)
        for k in (2, 3, 4):
            Dk = 2**(k*len(w))-a**(k*sum(w))
            H = Dk//D
            check(H > 0 and Dk == H*D, 'positive_power_factor')
            check(carry(w*k, a) == H*B, 'power_carry')
            check(gcd(Dk, carry(w*k, a)) == H*gcd(D, B), 'power_content')
            check(abs(Dk)//gcd(Dk, carry(w*k, a)) == abs(D)//gcd(D, B),
                  'power_reduced_denominator')
            powers += 1

out = {
    'status': 'PASS',
    'scope': 'Independent finite controls, not the all-length proof and not a convergence proof.',
    'universe': {
        'all_nonconstant_binary_words': 'lengths 2 through 9, bases 3, 5, 9',
        'additive_parameters': [-5, -3, -1, 1, 3, 5],
        'connected_graph': 'binary-parent tree with independently varied allowed context stripping',
        'extra_signed_word': '11110111000',
        'extremal_shape_base_cases': shape_cases,
        'word_power_cases': powers,
    },
    'counts': dict(sorted(counts.items())),
    'certificate_digest_sha256': digest.hexdigest(),
}
print(json.dumps(out, indent=2, sort_keys=True))
