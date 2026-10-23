"""Exact controls for the recovered valuation chain and its decorated quotient.

No external packages, randomness, or optimization-sensitive assertions.
"""
import itertools
import json
from fractions import Fraction
from math import prod

checks = 0


def check(predicate, detail):
    global checks
    checks += 1
    if not predicate:
        raise RuntimeError(detail)


def valuation(n, p):
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count, n


def register(g, p, k):
    an, un = valuation(g.numerator, p)
    ad, ud = valuation(g.denominator, p)
    modulus = p**k
    return an-ad, (un * pow(ud, -1, modulus)) % modulus


def register_product(left, right, modulus):
    return left[0]+right[0], left[1]*right[1] % modulus


def quotient_period(edges, start):
    # The only endpoint identification is 4 ~ 13. Other states stay marked.
    representative = lambda v: 13 if v == 4 else v
    image = {}
    for a, b in edges.items():
        a, b = representative(a), representative(b)
        if a in image and image[a] != b:
            raise RuntimeError('This set quotient does not carry a function')
        image[a] = b
    visited = {}
    x = representative(start)
    while x not in visited:
        visited[x] = len(visited)
        x = image[x]
    return len(visited)-visited[x]


core = 2**5*3**7*5**3*7**4*11*13**2
edges = {p: valuation(core, p)[0] for p in (2, 3, 5, 7, 11, 13)}
check(core == 39046339332000, 'Core arithmetic')
check(edges == {2: 5, 3: 7, 5: 3, 7: 4, 11: 1, 13: 2}, 'Valuation arrows')
check(quotient_period(edges, 13) == 5, 'Five quotient states')
normalized = {p: valuation(core//2, p)[0] for p in edges}
check(quotient_period(normalized, 13) == 2, 'Normalized two-cycle')
for modulus, a, b in ((9, 11, 2), (3, 2, 5)):
    check((a-b) % modulus == 0, 'Input collision')
    check((edges[a]-edges[b]) % modulus != 0, 'Output disagreement')

chain = (13, 2, 5, 3, 7, 4)
gains = [Fraction(b, a) for a, b in zip(chain, chain[1:])]
debt = prod(gains)
check(debt == Fraction(4, 13), 'Exact five-edge holonomy')
check(debt*Fraction(13, 4) == 1, 'Explicit reset cost')
check(prod((Fraction(2, 13), Fraction(4, 2))) == debt, 'Period does not fix debt')
check(register(debt, 2, 3) == (2, 5), 'Mod-8 valuation/unit record')
check(register(debt, 3, 2) == (0, 1), 'Mod-9 invisible debt')
for p, k in ((2, 3), (3, 2)):
    check(register(Fraction(73), p, k) == (0, 1), 'Joint-register hostile')

word_count = 0
for length in range(2, 7):
    for labels in itertools.product(range(1, 6), repeat=length):
        word_count += 1
        local = [Fraction(b, a) for a, b in zip(labels, labels[1:])]
        h = [Fraction(labels[0], a) for a in labels[:-1]]
        transformed = [h[(i+1) % len(h)]*g/h[i] for i, g in enumerate(local)]
        endpoint = Fraction(labels[-1], labels[0])
        check(prod(local) == endpoint, 'Telescoping')
        check(transformed == [Fraction(1)]*(len(local)-1)+[endpoint], 'Rational seam gauge')
        check((prod(local) == 1) == (labels[-1] == labels[0]), 'Nonzero vertex potential gate')

fraction_list = [Fraction(a, b) for a in range(1, 9) for b in range(1, 9)]
pair_count = 0
for g, h in itertools.product(fraction_list, repeat=2):
    pair_count += 1
    for p, k in ((2, 3), (3, 2)):
        check(register(g*h, p, k) == register_product(register(g, p, k), register(h, p, k), p**k), 'Register homomorphism')

phase_counts = [0, 0, 0]
hamiltonian_count = 0
for order in itertools.permutations(range(6)):
    asc = sum(a < b for a, b in zip(order, order[1:]))
    phase_counts[asc % 3] += 1
    hamiltonian_count += asc == 5
all_coefficients = [phase_counts[0]-phase_counts[2], phase_counts[1]-phase_counts[2]]
check(sum(phase_counts) == 720, 'All permutation universe')
check(hamiltonian_count == 1, 'Transitive Hamiltonian subset')
check(all(c % 9 == 0 for c in all_coefficients), 'All-permutation phase cancellation')
check(hamiltonian_count % 9 != 0, 'Restricted phase sum is not divisible by nine')

def integral_clock(source, word):
    current = Fraction(source)
    for steps in range(1000):
        next_value = current*word[steps % len(word)]
        if next_value.denominator != 1:
            return steps
        current = next_value
    raise RuntimeError('Finite control did not exit within its independent bound')


integer_sources = set(range(1, 513)) | {
    13**k*b for k in range(13) for b in (1, 2, 3, 5, 7, 11, 17, 19)
}
for source in sorted(integer_sources):
    budget = valuation(source, 13)[0]
    check(integral_clock(source, gains) == 5*budget, 'Five-edge source-specific clock')
    check(integral_clock(source, [Fraction(2, 13), Fraction(2)]) == 2*budget,
          'Two-edge source-specific clock')

print(json.dumps({
    'status': 'PASS',
    'checks': checks,
    'universe': {
        'label_words': word_count,
        'labels': [1, 2, 3, 4, 5],
        'word_lengths': [2, 3, 4, 5, 6],
        'ordered_fraction_pairs': pair_count,
        'fraction_numerators_denominators': [1, 8],
        'registers': [[2, 3], [3, 2]],
    },
    'core': core,
    'quotient_periods': [5, 2],
    'chain': chain,
    'gains': [str(g) for g in gains],
    'holonomy': str(debt),
    'reset': '13/4',
    'holonomy_registers': {'2^3': register(debt, 2, 3), '3^2': register(debt, 3, 2)},
    'joint_register_hostile': 73,
    'transitive_order_six': {
        'phase_counts_mod_three': phase_counts,
        'all_permutation_coefficients_in_one_omega': all_coefficients,
        'hamiltonian_count': hamiltonian_count,
        'hamiltonian_coefficients_in_one_omega': [-1, -1],
    },
    'integer_lift_controls': {
        'distinct_sources': len(integer_sources),
        'consecutive_sources': [1, 512],
        'exponents_of_13': [0, 12],
        'coprime_bases': [1, 2, 3, 5, 7, 11, 17, 19],
        'successful_edge_clocks': ['5*v_13(s)', '2*v_13(s)'],
    },
}, indent=2, sort_keys=True))
