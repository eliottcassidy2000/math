#!/usr/bin/env python3
"""Independent exact controls for the height-controlled Collatz entry proof.

These controls authenticate finite identities; the all-size theorem is proved
in the reviewed mathematical argument, not by a bounded numerical census.
No production checker is imported or executed.
"""
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json

checks = 0

def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(label)

def step(n):
    return (3*n+1)//2 if n % 2 else n//2

word = '11011011011011010101101110101010'
# Reconstruct the unique lawful residue by successive binary lifting, using
# only direct integer iteration; do not infer it from endpoint divisibility.
residue = 0
for j, bit in enumerate(word):
    options = []
    for x in (residue, residue + (1 << j)):
        y = x
        prefix = ''
        for _ in range(j + 1):
            prefix += str(y % 2)
            y = step(y)
        if prefix == word[:j + 1]:
            options.append(x)
    require(len(options) == 1, 'unique direct binary lift')
    residue = options[0]
require(residue == 2962948091, 'claimed word residue')
require(len(word) == 32 and word.count('1') == 20, 'word length/count')
modulus, multiplier = 1 << 32, 3**20
y = residue
for _ in word:
    y = step(y)
carry = modulus*y - multiplier*residue
require((multiplier, carry, y) == (3486784401, 15810005461, 2405410912),
        'affine branch from direct orbit')
require(multiplier < modulus and carry < 20*(modulus-multiplier),
        'fixed point below 20')
require(residue > 20, 'all positive lawful lifts descend')
for lift in (0, 1, 2, 17, 1234567):
    x = residue + lift*modulus
    y = x
    observed = ''
    for _ in word:
        observed += str(y % 2)
        y = step(y)
    require(observed == word, 'direct replay of positive lift')
    require(modulus*y == multiplier*x + carry and y < x,
            'lift affine identity and descent')

# Independent direct-orbit check of parity bijection and BOTH moments.
source_count = 0
for n in range(1, 15):
    size = 1 << n
    seen = set()
    sum_multiplier = [0]*(n+1)
    sum_carry = [0]*(n+1)
    avoidance = Counter()
    for x in range(size, 2*size):
        source_count += 1
        y, h, prefix = x, 0, ''
        sum_multiplier[0] += 1
        for j in range(1, n+1):
            bit = y % 2
            prefix += str(bit)
            h += bit
            y = step(y)
            numerator_a = 3**h
            numerator_b = (1 << j)*y - numerator_a*x
            require(numerator_b >= 0, 'nonnegative exact orbit carry')
            sum_multiplier[j] += numerator_a
            sum_carry[j] += numerator_b
        require(prefix not in seen, 'direct-source parity injection')
        seen.add(prefix)
        # Two tiny targets, embedded away from time zero when space permits.
        for target in ('0', '10', '110'):
            L = len(target)
            J = n//4
            k = (n-J)//L
            if all(prefix[J+t*L:J+(t+1)*L] != target for t in range(k)):
                avoidance[target] += 1
    require(len(seen) == size, 'parity bijection')
    for j in range(n+1):
        require(sum_multiplier[j] == size*(1 << j), 'mean multiplier is one')
        require(4*sum_carry[j] == j*size*(1 << j), 'mean carry equals j/4')
    for target in ('0', '10', '110'):
        L = len(target)
        J = n//4
        k = (n-J)//L
        expected = (1 << (n-L*k))*((1 << L)-1)**k
        require(avoidance[target] == expected, 'unconditional block avoidance')

q20 = Fraction(10**20, 9**9*11**11)
require(0 < q20 < Fraction(999, 1000)**20 < 1, 'exact Chernoff bound')
require(3**5 < 2**8, 'rational logarithm bound')
require(Fraction(11, 20)*Fraction(8, 5)-1 == -Fraction(3, 25),
        'multiplier exponent')
require(Fraction(3, 25)*Fraction(3, 4) == Fraction(9, 100),
        'entry-time multiplier contraction')
require(Fraction(15, 16)-Fraction(91, 100) == Fraction(11, 400),
        'height exponent margin')
require(73*Fraction(11, 400) >= 2 and 72*Fraction(11, 400) < 2,
        'exact height-comparison onset')
require(Fraction(1, 2**38) < Fraction(3, 4000), 'separated error exponents')

def parameters(m):
    N = 7*m//8
    J = (3*m+3)//4
    return N, J, (N-J)//32

# The rounding functions are affine-periodic with period 256. Checking a
# complete residue system plus the identities proves these rounding claims
# for every integer m>=512, not just the finite range checked.
minimum_k = 10**9
for m in range(512, 768):
    N, J, k = parameters(m)
    N1, J1, k1 = parameters(m+256)
    require((N1-N, J1-J, k1-k) == (224, 192, 1), 'rounding translation')
    require(0 < J < N <= m and k >= 1, 'nonempty inspected block set')
    require(J+32*k <= N, 'full final word inside parity horizon')
    require(Fraction(k) >= Fraction(m, 256)-2, 'global floor loss bound')
    minimum_k = min(minimum_k, k)

with localcontext() as context:
    context.prec = 65
    D = Decimal
    p = D(2) ** -32
    delta = -(1-p).ln()/D(2).ln()/D(256)
    q = (D(10**20)/D(9**9*11**11)).ln().exp() ** (D(1)/D(20))
    chernoff_exponent = -D(3)/D(4)*q.ln()/D(2).ln()
    require(0 < delta < D(2)**-38, 'display delta within proved rational bound')
    require(chernoff_exponent > D('0.0054'), 'display Chernoff exponent')
    display = {'delta': str(delta), 'q': str(q),
               'chernoff_exponent': str(chernoff_exponent)}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--proof', type=Path, help='Path to the audited HEIGHT_CONTROLLED_ENTRY.md')
args = parser.parse_args()
base = Path(__file__).resolve().parent.parent
candidates = [
    base / 'HEIGHT_CONTROLLED_ENTRY.md',
    base / 'collatz_height_controlled_entry_20261019/HEIGHT_CONTROLLED_ENTRY.md',
    base / 'checkpoint2/05-knowledge/results/collatz_height_controlled_entry_20261019/HEIGHT_CONTROLLED_ENTRY.md',
]
proof = args.proof or next((candidate for candidate in candidates if candidate.is_file()), None)
if proof is None:
    raise RuntimeError('Cannot locate the proof; pass --proof PATH explicitly.')
print(json.dumps({
    'status': 'PASS',
    'checked_proof_sha256': hashlib.sha256(proof.read_bytes()).hexdigest(),
    'checks': checks,
    'finite_source_controls': source_count,
    'finite_source_range': 'complete dyadic intervals of exponents 1 through 14',
    'rounding_control': 'complete period 512<=m<768, with exact 256-translation',
    'minimum_candidate_count_in_base_period': minimum_k,
    'word_residue': residue,
    'scope': 'Independent finite controls. All-size proof is the written audit.',
    **display,
}, indent=2, sort_keys=True))
