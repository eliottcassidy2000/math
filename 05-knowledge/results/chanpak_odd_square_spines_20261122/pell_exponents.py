"""Sharp two-sided Pell defects for all odd Euclid exponent phases."""
from collections import defaultdict
from math import gcd

from pell_directions import MODULUS, PRIMES, advance, phase


EXPONENTS = (1, 5, 21, 29, 33, 43, 49, 61, 71, 75, 77, 85, 103)


def optimum(chart_order, primes, target):
    data = phase.chart(chart_order)
    tables = [phase.local_counts(p, data, target % p) for p in primes]
    sizes = phase.all_sizes(tables)
    best = min(sizes.values())
    winners = {sig for sig, size in sizes.items() if size == best}
    # These marginal sets are necessary; final witnesses retain the full joint test.
    margins = {}
    for i, (p, table) in enumerate(zip(primes, tables)):
        winning_signatures = {sig[i] for sig in winners}
        margins[p] = {direction for direction, sig in table.items() if sig in winning_signatures}

    def optimal(a, b):
        sig = tuple(table[phase.normalized(a, b, p)] for table, p in zip(tables, primes))
        return sig in winners

    return margins, optimal


def first_in_period(seed, w4, w8):
    a, b = seed
    for k in range(480):
        if w4(a,b) and w8(a,b):
            return k
        a,b = advance(a,b,MODULUS)
    assert all(phase.normalized(a,b,p) == phase.normalized(*seed,p) for p in PRIMES)
    return None


def main():
    groups = defaultdict(list)
    squares = {(1, x*x % 29) for x in range(1,29)}
    nonsquares = {(1, x) for x in range(1,29)} - squares
    old_bank = []
    for residue in range(1,112,2):
        e = 113 if residue == 1 else residue
        target = (1 << (e-1))*((1 << e)-1)-1
        m4, w4 = optimum(4, (2,5,17,29), target)
        m8, w8 = optimum(8, (2,5,29,257), target)
        if w4(7,10) and w8(7,10):
            old_bank.append(residue)
        common = m4[29] & m8[29]
        r = e % 28
        if r in (7,11,13,23):
            assert not common
            groups["incompatible at 29"].append(residue)
            continue
        positive_character = r in (3,9,17,25,27)
        assert common == (squares if positive_character else nonsquares)
        assert all(m4[p] & m8[p] == {(1,0)} for p in (2,5))
        slopes = range(1,17) if e%8 == 1 else {x*x%17 for x in range(1,17)}
        assert m4[17] == {(1,s) for s in slopes}

        extrema = []
        receipts = []
        for unit, norm31 in [((1,2),(1,8)), ((1,0),(7,6))]:
            seed, k = unit, first_in_period(unit, w4, w8)
            if k is None:
                seed, k = norm31, first_in_period(norm31,w4,w8)
            assert k is not None
            a, b = seed
            for _ in range(k):
                a,b = 3*a+2*b,4*a+3*b
            assert a > 0 and b > 0 and gcd(a,b) == 1 and w4(a,b) and w8(a,b)
            extrema.append(2*a*a-b*b)
            receipts.append((seed,k))
        expected = (-62,62) if positive_character else (-2, 2 if e%8 == 1 else 62)
        assert tuple(extrema) == expected
        groups[expected].append(residue)
        print("exponent", residue, "extremal defects", expected, "seeds and powers", receipts)
    assert tuple(old_bank) == EXPONENTS
    for key, residues in groups.items():
        print(key, len(residues), "classes:", residues)


if __name__ == "__main__":
    main()
