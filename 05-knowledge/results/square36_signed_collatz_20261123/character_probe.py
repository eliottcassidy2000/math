"""Legendre square detectors and Hamilton path counts for signed Collatz itineraries."""
from functools import lru_cache
from itertools import combinations, permutations
from math import isqrt

DETECTOR = (79, 283, 331)
HEEGNER = (3, 7, 11, 19, 43, 67, 163)


def prime(p):
    return p >= 2 and all(p % q for q in range(2, isqrt(p) + 1))


def character(value, p):
    residue = pow(value % p, (p - 1) // 2, p)
    return -1 if residue == p - 1 else residue


def first_detectors(primes, n):
    false_sums = [x for x in range(3, 2*n) if isqrt(x)**2 != x]
    # A bit survives if the prime accepts this nonsquare as a local square.
    bad = {p: sum((character(x, p) != -1) << j for j, x in enumerate(false_sums))
           for p in primes}
    for k in range(1, len(primes) + 1):
        good = []
        for cols in combinations(primes, k):
            survivors = (1 << len(false_sums)) - 1
            for p in cols:
                survivors &= bad[p]
            if not survivors:
                good.append(cols)
        if good:
            return good
    return []


def square_code(n):
    """Square-class code witnessing the three-character threshold at 143."""
    syndrome = 0
    divisor = 2
    while divisor <= n:
        code = {2: 1, 3: 2, 5: 4, 7: 7}.get(divisor, 6)
        while n % divisor == 0:
            n //= divisor
            syndrome ^= code
        divisor += 1
    return syndrome


def majority_rows(labels):
    rows = []
    for x in labels:
        mask = 0
        for j, y in enumerate(labels):
            if x == y:
                continue
            signs = [character(y - x, p) for p in DETECTOR]
            assert 0 not in signs
            if sum(signs) > 0:
                mask |= 1 << j
        rows.append(mask)
    return tuple(rows)


def ham(rows):
    @lru_cache(None)
    def continuations(unvisited, last):
        if not unvisited:
            return 1
        choices = rows[last] & unvisited
        total = 0
        while choices:
            bit = choices & -choices
            total += continuations(unvisited ^ bit, bit.bit_length()-1)
            choices ^= bit
        return total

    full = (1 << len(rows)) - 1
    return sum(continuations(full ^ (1 << i), i) for i in range(len(rows)))


def isomorphism(first, second, perm):
    return all((first[i] >> j) & 1 == (second[perm[i]] >> perm[j]) & 1
               for i in range(len(first)) for j in range(len(first)))


def minimum_classes():
    first = majority_rows((1, 3, 5, 6, 7))
    second = majority_rows((1, 3, 5, 6, 9))
    assert isomorphism(first, second, (2, 0, 1, 4, 3))
    masses = {11: 0, 13: 0, 15: 0}
    for n in (4, 5):
        pairs = tuple(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            rows = [0]*n
            for j, (a, b) in enumerate(pairs):
                if (mask >> j) & 1:
                    rows[a] |= 1 << b
                else:
                    rows[b] |= 1 << a
            h = ham(rows)
            if n == 4:
                assert h <= 5
            elif h in masses:
                masses[h] += 1
    assert masses == {11: 120, 13: 120, 15: 64}
    upper = majority_rows(tuple(6*x+1 for x in (1, 3, 5, 6, 9)))
    regular = tuple(sum(1 << j for j in range(5) if (j-i) % 5 in (1, 2)) for i in range(5))
    for rows, order in [(upper, 3), (regular, 5)]:
        assert ham(rows) == 15
        assert sum(isomorphism(rows, rows, perm) for perm in permutations(range(5))) == order
    print("order-five labeled masses", masses, "minimum class counts", (1, 1, 2))


def main():
    minimum_classes()
    inert = [p for p in range(3, 332, 4) if prime(p)]
    detectors = first_detectors(inert, 36)
    assert len(detectors) == 7 and all(t[1:] == (283, 331) for t in detectors)
    assert detectors[0] == DETECTOR
    heegner = first_detectors(HEEGNER, 36)
    assert len(heegner) == 4 and all(len(t) == 5 for t in heegner)
    print("least maximal inert prime among optimal triples", detectors)
    print("minimum Heegner subfamilies", heegner)
    assert all((square_code(n) == 0) == (isqrt(n)**2 == n) for n in range(1, 143))
    assert square_code(143) == 0
    assert [(character(t, 79), character(t, 283), character(t, 331)) for t in (2, 3, -3)] == [
        (1, -1, -1), (-1, -1, -1), (1, 1, 1)]
    # The two legal steps share exactly the same action on pair differences.
    for labels, expected in [
        (tuple(range(1, 10)), (2831, 3159)),
        ((1, 3, 5, 6, 7), (13, 11)),
        ((1, 3, 5, 6, 9), (13, 15)),
    ]:
        before = majority_rows([4*x + 1 for x in labels])
        assert before == majority_rows(labels) == majority_rows([4*x + 3 for x in labels])
        after = majority_rows([6*x + 1 for x in labels])
        assert after == majority_rows([6*x + 5 for x in labels])
        pair = (ham(before), ham(after))
        assert pair == expected
        changed = [(x, y) for i, x in enumerate(labels) for j, y in enumerate(labels[i+1:], i+1)
                   if (before[i] ^ after[i]) & (1 << j)]
        if len(labels) == 5:
            assert changed == [(1, 6)]
        print("labels", labels, "counts", pair, "flipped arcs", changed)


if __name__ == "__main__":
    main()
