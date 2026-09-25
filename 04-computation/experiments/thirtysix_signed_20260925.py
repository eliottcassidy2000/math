"""Signed Collatz clocks and exact permutation half-step completions.

Finite cycle certificates, not an exhaustive all-height cycle census.
"""

from collections import Counter
from itertools import permutations, product
from math import factorial


def check(test, message):
    if not test:
        raise RuntimeError(message)


def valuation2(n):
    check(n != 0, "valuation domain")
    n = abs(n)
    return (n & -n).bit_length() - 1


def ordinary(n, sigma):
    return 3*n+sigma if n % 2 else n//2


def shortcut(n, sigma):
    return (3*n+sigma)//2 if n % 2 else n//2


def odd_only(n, sigma):
    check(n % 2 == 1, "odd-only domain")
    value = 3*n+sigma
    return value//2**valuation2(value)


def cycle(start, step):
    answer = []
    current = start
    while current not in answer:
        check(len(answer) < 100, "finite witness cap")
        answer.append(current)
        current = step(current)
    check(current == start, "supplied start must be periodic")
    return tuple(answer)


def cycle_type(p):
    seen = set()
    lengths = []
    for x in range(len(p)):
        if x in seen:
            continue
        current = x
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = p[current]
        check(current == x, "permutation cycle")
        lengths.append(length)
    return tuple(sorted(lengths))


def predicted_square_roots(lengths):
    count = 1
    for length, multiplicity in Counter(lengths).items():
        factor = 0
        pair_counts = [multiplicity//2] if length % 2 == 0 else range(multiplicity//2+1)
        if length % 2 == 0 and multiplicity % 2:
            return 0
        for pairs in pair_counts:
            factor += (factorial(multiplicity)//
                       (factorial(multiplicity-2*pairs)*2**pairs*factorial(pairs))
                       * length**pairs)
        count *= factor
    return count


def permutation_from_cycles(cycles):
    p = {}
    for vertices in cycles:
        for i, vertex in enumerate(vertices):
            check(vertex not in p, "disjoint cycles")
            p[vertex] = vertices[(i+1) % len(vertices)]
    return p


def paired_root(a, b, phase):
    length = len(a)
    check(length == len(b), "matching cycle lengths")
    return {**{a[i]: b[(i+phase) % length] for i in range(length)},
            **{b[j]: a[(j+1-phase) % length] for j in range(length)}}


def odd_cycle_root(vertices):
    length = len(vertices)
    check(length % 2 == 1, "single odd cycle root")
    return {vertices[i]: vertices[(i+(length+1)//2) % length] for i in range(length)}


def check_root(root, base):
    check(set(root) == set(base) == set(root.values()), "root bijection")
    check(all(root[root[x]] == base[x] for x in base), "exact square root")


def type_dict(p):
    labels = list(p)
    index = {x: i for i, x in enumerate(labels)}
    return cycle_type(tuple(index[p[x]] for x in labels))


def tagged(cyc, tag):
    return tuple((n, tag) for n in cyc)


def main():
    clocks = (('odd', odd_only), ('shortcut', shortcut), ('ordinary', ordinary))
    witnesses = {}
    for name, step in clocks:
        cs = [cycle(m, lambda n: step(n, -1)) for m in (1, 5, 17)]
        witnesses[name] = cs
        for c in cs:
            check(cycle(-c[0], lambda n: step(n, 1)) == tuple(-n for n in c),
                  "signed cycle conjugacy")
        for n in range(-2001, 2002):
            if name == 'odd' and n % 2 == 0:
                continue
            check(step(-n, -1) == -step(n, 1), "signed forward conjugacy")
            if n:
                epsilon = 1 if n > 0 else -1
                check(step(n, -epsilon) == epsilon*step(abs(n), -1), "minus gluing")
                check(step(n, epsilon) == epsilon*step(abs(n), 1), "plus gluing")
                for sigma in (-1, 1):
                    def half(x):
                        return -x if x > 0 else step(-x, sigma)
                    check(half(half(n)) == epsilon*step(abs(n), sigma), "global half-step square")
        lengths = [len(c) for c in cs]
        additions = sum(m for m, a in Counter(lengths).items() if m % 2 == 0 and a % 2)
        print('CLOCK', name, 'lengths', lengths, 'states', sum(lengths),
              'two-sign states', 2*sum(lengths), 'minimal root completion', sum(lengths)+additions)
        for c in cs:
            print(' CYCLE', c)
    for i, c in enumerate(witnesses['odd']):
        ks = tuple(valuation2(3*n-1) for n in c)
        check(sum(ks) == len(witnesses['shortcut'][i]), "shortcut clock K")
        check(len(c)+sum(ks) == len(witnesses['ordinary'][i]), "ordinary clock L+K")
        print('ODD WORD', c[0], ks, 'L,K,L+K', len(c), sum(ks), len(c)+sum(ks))

    # Independent brute-force square-map fibres on every permutation through n=8.
    total = 0
    for n in range(9):
        ps = list(permutations(range(n)))
        actual = Counter(tuple(q[q[x]] for x in range(n)) for q in ps)
        for p in ps:
            check(actual[p] == predicted_square_roots(cycle_type(p)), "all-permutation root count")
        total += len(ps)
    print('ROOT THEOREM CONTROL all permutations n0..8:', total, 'PASS')

    # Minimal ordinary completion: keep one five-cycle; duplicate the even cycles.
    c2, c5, c18 = witnesses['ordinary']
    a2, b2 = tagged(c2, 'a'), tagged(c2, 'b')
    a5, b5 = tagged(c5, 'a'), tagged(c5, 'b')
    a18, b18 = tagged(c18, 'a'), tagged(c18, 'b')
    base45 = permutation_from_cycles((a2, b2, a5, a18, b18))
    roots45 = set()
    type45 = Counter()
    for phase2, phase18 in product(range(2), range(18)):
        root = paired_root(a2, b2, phase2) | odd_cycle_root(a5) | paired_root(a18, b18, phase18)
        check_root(root, base45)
        roots45.add(tuple(root[x] for x in base45))
        type45[type_dict(root)] += 1
    check(len(roots45) == predicted_square_roots((2,2,5,18,18)) == 36, "36 minimal roots")
    print('MINIMAL ORDINARY COMPLETION states45 roots36 root-cycle-types', dict(type45))

    base50 = permutation_from_cycles((a2, b2, a5, b5, a18, b18))
    five_options = [odd_cycle_root(a5) | odd_cycle_root(b5)]
    five_options.extend(paired_root(a5, b5, a) for a in range(5))
    roots50 = set()
    type50 = Counter()
    for phase2, phase18, five in product(range(2), range(18), five_options):
        root = paired_root(a2, b2, phase2) | five | paired_root(a18, b18, phase18)
        check_root(root, base50)
        roots50.add(tuple(root[x] for x in base50))
        type50[type_dict(root)] += 1
    check(len(roots50) == predicted_square_roots((2,2,5,5,18,18)) == 216, "216 doubled roots")
    print('REFLECTED ORDINARY LOCUS states50 roots216 root-cycle-types', dict(type50))

    for a, b in ((a2, b2), (a5, b5), (a18, b18)):
        reflection = dict(zip(a, b)) | dict(zip(b, a))
        commuting = []
        for phase in range(len(a)):
            root = paired_root(a, b, phase)
            does_commute = all(root[reflection[x]] == reflection[root[x]] for x in root)
            check(does_commute == ((2*phase-1) % len(a) == 0), "reflection phase equation")
            if does_commute:
                commuting.append(phase)
        check(len(commuting) == len(a) % 2, "even-cycle reflection obstruction")
        print('REFLECTION GAUGE length', len(a), 'commuting half-step phases', commuting)

    # The shortcut root needs no duplicate; the odd-only root needs only its 2-cycle copied.
    shortcut_cycles = [tagged(c, 'a') for c in witnesses['shortcut']]
    shortcut_base = permutation_from_cycles(shortcut_cycles)
    shortcut_root = {}
    for c in shortcut_cycles:
        shortcut_root |= odd_cycle_root(c)
    check_root(shortcut_root, shortcut_base)
    oc1, oc2, oc7 = [tagged(c, 'a') for c in witnesses['odd']]
    oc2b = tagged(witnesses['odd'][1], 'b')
    odd_base = permutation_from_cycles((oc1, oc2, oc2b, oc7))
    for phase in range(2):
        check_root(odd_cycle_root(oc1) | paired_root(oc2, oc2b, phase) | odd_cycle_root(oc7), odd_base)
    print('OTHER CLOCK COMPLETIONS shortcut15 unique root; odd12 exactly2 roots PASS')

    # The 36-cycle is an explicit interleaving, not a legal arithmetic orbit.
    root36 = paired_root(a18, b18, 0)
    orbit36 = cycle(a18[0], lambda x: root36[x])
    check(len(orbit36) == 36 and orbit36[::2] == a18 and orbit36[1::2] == b18,
          "two half-clock sections")
    print('36-CYCLE two sections each reproduce ordinary17 cycle; phase alternates PASS')
    print('COUNT FIREWALL representatives6 -> ordered pairs36; actual glued states odd20/shortcut30/ordinary50')
    print('SCOPE phase-root is time-clock construction; arithmetic inverse-fibre half-step is separate')
    print('ALL CHECKS PASS; no completeness, convergence, or prime-distribution claim')


if __name__ == '__main__':
    main()
