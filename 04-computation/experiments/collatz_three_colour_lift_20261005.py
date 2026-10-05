"""Exact four-state decorations, source likelihoods, and payment interfaces.

The independent colours are not the owner's ordered Zeckendorf stripe.
All proof checks use integers/Fractions and remain active under -O.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


def charge_counts(r):
    need(type(r) is int and r >= 0, 'nonnegative exact one-count')
    return ((3**r+3*(-1)**r)//4,)+((3**r-(-1)**r)//4,)*3


def rotate(q):
    a, b = q & 1, q >> 1
    return b+2*(a ^ b)


def shortcut(n):
    need(type(n) is int and n > 0, 'positive exact integer')
    return (3*n+1)//2 if n % 2 else n//2


def prefix(n, length):
    need(type(length) is int and length >= 0, 'nonnegative exact horizon')
    out = ''
    for _ in range(length):
        n = shortcut(n)
        out += str(n % 2)
    return out, n


def affine_prefix(word, initial=1):
    need(type(word) is str and set(word) <= {'0', '1'}, 'binary prefix')
    need(type(initial) is int and initial in (0, 1), 'initial parity')
    p, q, b, parity = 1, 1, 0, initial
    for next_bit in word:
        if parity:
            p, b = 3*p, 3*b+q
        q *= 2
        parity = int(next_bit)
    return p, q, b


def valuation_carrier(word):
    p, q, b = 1, 1, 0
    for a in word:
        need(type(a) is int and a >= 1, 'exact positive valuation')
        p, q, b = 3*p, q*(1 << a), 3*b+q
    return p, q, b


def replay(n, word):
    need(type(n) is int and n > 0 and n % 2, 'positive odd exact source')
    for expected in word:
        raw = 3*n+1
        actual = (raw & -raw).bit_length()-1
        need(actual == expected, 'actual source valuation guard')
        n = raw >> actual
    return n


def guarded_charge_partition(words):
    """Two scalar partition sums recover all four aggregate-charge counts."""
    s = sum(3**w.count('1') for w in words)
    t = sum((-1)**w.count('1') for w in words)
    return s, t, ((s+3*t)//4,)+((s-t)//4,)*3


def charge_dp(word):
    counts = [1, 0, 0, 0]
    for bit in word:
        choices = (0,) if bit == '0' else (1, 2, 3)
        new = [0]*4
        for old, multiplicity in enumerate(counts):
            for q in choices:
                new[old ^ q] += multiplicity
        counts = new
    return tuple(counts)


def zeck_support(n, weights):
    out = []
    for i in range(len(weights)-1, 1, -1):
        if weights[i] <= n:
            out.append(i)
            n -= weights[i]
    need(n == 0, 'adequate Fibonacci support')
    return tuple(sorted(out))


def owner_controls():
    weights = [1, 1, 1, 2]
    while weights[-1] <= 100:
        weights.append(weights[-1]+weights[-2])
    charges = [(-1, 1), (1, 0), (1, 0), (0, 1)]
    while len(charges) < len(weights):
        x, y = charges[-2:]
        charges.append((x[0]+y[0], x[1]+y[1]))

    def charge(indices):
        return tuple(sum(charges[i][j] for i in indices) for j in (0, 1))

    fibres = {n: set() for n in range(1, 101)}
    for mask in range(1 << len(weights)):
        if mask & (mask << 1):
            continue
        support = tuple(i for i in range(len(weights)) if mask >> i & 1)
        n = sum(weights[i] for i in support)
        if n in fibres:
            fibres[n].add(support)
    sizes = Counter()
    for n, fibre in fibres.items():
        ordinary = zeck_support(n, weights)
        previous = zeck_support(n-1, weights)
        expected = {ordinary, (0,)+previous}
        if 2 not in previous:
            expected.add((1,)+previous)
        need(fibre == expected, 'complete marked representation fibre')
        need(sum(charge(s) == charge(ordinary) for s in fibre) == 2,
             'exactly two full-charge-preserving representations')
        sizes[len(fibre)] += 1
    q2, q4 = (charge(zeck_support(n, weights)) for n in (2, 4))
    need(tuple(c % 2 for c in q2) == tuple(c % 2 for c in q4),
         'same auxiliary charge does not restore owner colours')
    need(charge((3,)) != (2, 0), '2 versus repeated1+1 charge hostile')
    words = ['K', 'B', 'R', 'RK', 'RKB']
    while len(words) <= 9:
        words.append(words[-1]+words[-2])
    supplied = 'RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB'
    need(words[9]+'B' == supplied, 'exact owner row35')
    need({''.join(words[i] for i in reversed(s)) for s in fibres[36]}
         == {words[9]+'RK'}, 'scoped row36 continuation obstruction')
    return dict(sorted(sizes.items()))


def main():
    report = {'owner_representation_fibre_sizes_n1_to_100': owner_controls()}
    for r in range(9):
        actual = Counter()
        for word in product((1, 2, 3), repeat=r):
            q = 0
            for colour in word:
                q ^= colour
            actual[q] += 1
        need(tuple(actual[q] for q in range(4)) == charge_counts(r),
             'independent exact colour enumeration')
    report['nonzero_colour_tuple_universe'] = sum(3**r for r in range(9))
    decoration_checks = 0
    for length in range(7):
        actual = Counter()
        for colours in product(range(4), repeat=length):
            bits = ''.join(str(int(q != 0)) for q in colours)
            charge = 0
            for q in colours:
                charge ^= q
                need(rotate(rotate(rotate(q))) == q and (rotate(q) == 0) == (q == 0),
                     'order-three action preserves the3+1 quotient')
            actual[bits, charge] += 1
            decoration_checks += 1
        for word in map(''.join, product('01', repeat=length)):
            need(tuple(actual[word, q] for q in range(4)) == charge_counts(word.count('1'))
                 == charge_dp(word), 'literal lift equals formula and transition matrix')
    report['all_four_state_words_through_length6'] = decoration_checks

    source_checks = 0
    for length in range(1, 11):
        seen = set()
        for n in range(1, 1 << (length+1), 2):
            word, target = prefix(n, length)
            seen.add(word)
            p, q, b = affine_prefix(word)
            need(F(p*n+b, q) == target, 'actual shortcut prefix equals carrier')
            likelihood = F(3**word.count('1'), 2**length)
            need(likelihood == F(p, q)*F(3)**(int(word[-1])-1),
                 'endpoint-parity likelihood correction')
            need(prefix(n+(1 << (length+1)), length)[0] == word,
                 'full binary source cylinder')
            source_checks += 1
        need(len(seen) == 1 << length, 'all output prefixes uniquely label odd residues')
    report['source_prefix_controls_all_odds_at_horizons1_to_10'] = source_checks
    composition_checks = 0
    for length in range(1, 9):
        for word in map(''.join, product('01', repeat=length)):
            p, q, b = affine_prefix(word)
            for cut in range(1, length):
                pu, qu, bu = affine_prefix(word[:cut])
                pv, qv, bv = affine_prefix(word[cut:], int(word[cut-1]))
                need((pv*pu, qv*qu, pv*bu+bv*qu) == (p, q, b),
                     'source and endpoint parity compose with the full carrier')
                composition_checks += 1
    report['every_binary_split_through_length8'] = composition_checks
    need(prefix(1, 1) == ('0', 2) and affine_prefix('0') == (3, 2, 1),
         'even-cut factor1/3 hostile')
    for r in range(1, 6):
        for word in product(range(1, 4), repeat=r):
            unary = ''.join('0'*(a-1)+'1' for a in word)
            need(affine_prefix(unary) == valuation_carrier(word),
                 'complete valuation boundaries erase parity correction')
    for length in range(1, 11):
        n = (1 << (length+1))-1
        word, target = prefix(n, length)
        need(word == '1'*length and target > n,
             'maximum decorated weight can accompany actual positive growth')
    report['three_bit_fibre_sizes_lex_order'] = [3**w.count('1')
        for w in map(''.join, product('01', repeat=3))]
    need(charge_dp('101') == charge_dp('011') == (3, 2, 2, 2),
         'same charge distribution can hide different ordered guards')
    need(valuation_carrier((1, 2)) == (9, 8, 5)
         and valuation_carrier((2, 1)) == (9, 8, 7)
         and prefix(11, 3)[0] == '101' and prefix(9, 3)[0] == '011',
         'word order and carry distinguish the source cells')

    # Existing paid primitives; coverage is only their exact source domains.
    bank = {
        'H': (729, 1024, 669, 155, 2048),
        'A': (81, 128, 85, 187, 256),
        'B': (81, 128, 73, 7, 256),
        'L': (9, 16, -3, 219, 256),
        'J': (243, 256, 147, 799, 1024),
    }
    selected = []
    for n in range(1, 2048, 2):
        hits = [letter for letter, (_, _, _, a, m) in bank.items() if n % m == a]
        need(len(hits) <= 1, 'five primitive paid guards are disjoint')
        if hits:
            selected.append(prefix(n, 10)[0])
    need(len(selected) == 27, 'exact existing paid guard source mass')
    s, t, charged = guarded_charge_partition(selected)
    summed = tuple(sum(charge_dp(w)[q] for w in selected) for q in range(4))
    need(charged == summed and sum(charged) == s, 'guarded two-scalar partition compiler')
    paid_checks = 0
    for p, q, b, a, m in bank.values():
        for k in range(32):
            n = a+k*m
            child = F(p*n+b, q)
            need(child.denominator == 1 and child.numerator % 2 and 0 < child < n,
                 'source-aware payment at inherited primitive guards')
            paid_checks += 1
    report['existing_paid_guard_bank_horizon10'] = {
        'letters': list(bank), 'odd_source_mass': '27/1024',
        'tilted_source_mass': str(F(s, 4**10)), 'S': s, 'T': t,
        'aggregate_charge_lift_counts': charged, 'literal_affine_payment_controls': paid_checks}
    extension = [w+bit for w in selected for bit in '01']
    s1, t1, charged1 = guarded_charge_partition(extension)
    need((s1, t1, charged1) == (4*s, 0, (s,)*4),
         'one unrestricted decoration erases every charge/finite-guard correlation')
    for incoming in range(4):
        need(Counter(incoming ^ q for q in range(4)) == Counter(range(4)),
             'unrestricted transfer matrix is rank-one all-ones')
    report['same_guard_bank_horizon11'] = {
        'S': s1, 'T': t1, 'aggregate_charge_lift_counts': charged1,
        'tilted_source_mass_unchanged': str(F(s1, 4**11))}
    for k in range(32):
        n = 799+1024*k
        child = (243*n+147)//256
        need(replay(n, (1, 1, 1, 1, 2, 3)) == replay(child, (1,)),
             'independent actual J common-future receipt')
        need(child < n, 'J pays despite source likelihood above1')
    need(F(729, 512)/F(3, 2) == F(243, 256), 'relative likelihood recovers child slope')
    report['J_interface_hostile'] = {'source': 799, 'child': 759, 'join': 1139,
        'source_likelihood': '729/512', 'child_likelihood': '3/2',
        'relative_likelihood': '243/256', 'actual_join_controls': 32}
    print(json.dumps(report, indent=2, sort_keys=True))
    print('PASS: exact colour lift and guarded comparisons; no new all-source coverage.')


if __name__ == '__main__':
    main()
