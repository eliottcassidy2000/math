"""Check the row-parity selector against planted permutations on eight labels."""
from collections import Counter
from fractions import Fraction
from itertools import permutations


def permutation_cycles(perm):
    unseen = set(range(len(perm)))
    cycles = []
    while unseen:
        start = min(unseen)
        point = start
        cycle = []
        while point in unseen:
            cycle.append(point)
            unseen.remove(point)
            point = perm[point]
        assert point == start
        cycles.append(cycle)
    return cycles


def check_cover_types():
    # The first four labels are columns and the last four are symbols.
    fixed = {(i, 4 + i) for i in range(4)}
    rows = {shift: fixed | {(i, 4 + (i + shift) % 4) for i in range(4)}
            for shift in (1, 2)}
    counts, sums = Counter(), Counter()
    for perm in permutations(range(8)):
        cycles = permutation_cycles(perm)
        kind = tuple(sorted(map(len, cycles)))
        if min(kind) < 3:
            continue
        counts[kind] += 1
        directed = {(i, perm[i]) for i in range(8)}
        for shift, union in rows.items():
            # A free fair orientation gives conditional expectation zero.
            term = 1
            for a, b in union:
                if (b, a) in directed:
                    term *= -1
                elif (a, b) not in directed:
                    term = 0
                    break
            sums[kind, shift] += term
    expect = {((8,), 1): Fraction(1, 2520), ((4, 4), 2): Fraction(1, 315)}
    assert counts == {(8,): 5040, (3, 5): 2688, (4, 4): 1260}
    for kind in sorted(counts):
        correlations = [Fraction(sums[kind, shift], counts[kind]) for shift in (1, 2)]
        assert correlations == [expect.get((kind, shift), 0) for shift in (1, 2)]
        print(f"type {kind}: {counts[kind]} permutations, shifts 1/2 "
              f"correlations {correlations[0]}, {correlations[1]}")


def check_trade():
    # In Z/12, exchange symbols 1 and 5 on one of two reflection orbits.
    n = 12
    table = [[(i + j) % n for j in range(n)] for i in range(n)]
    for i in range(n):
        if i % 4 in (0, 1):
            table[i] = [{1: 5, 5: 1}.get(v, v) for v in table[i]]
    assert all(len(set(row)) == n for row in table)
    assert all(len(set(col)) == n for col in zip(*table))
    assert all(table[i][j] == table[j][i] for i in range(n) for j in range(n))
    p, q = [[row.index(mark) for row in table] for mark in (0, 1)]
    # Traverse the two matching rows on disjoint shores directly.
    from_column = {i: n + p[i] for i in range(n)}
    from_symbol = {n + q[i]: i for i in range(n)}
    successor = from_column | from_symbol
    on_two_shores = tuple(sorted(map(len, permutation_cycles(successor))))
    assert on_two_shores == (8, 8, 8)
    print("symmetric Latin trade at B=12, A=1, C=5: doubled cover (8, 8, 8)")


if __name__ == "__main__":
    check_cover_types()
    check_trade()
