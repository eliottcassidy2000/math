"""Finite checks of the physical Mersenne faces and perfect-number neighbors."""

from collections import Counter
from functools import cache
from math import isqrt


@cache
def F(n):
    if n == 0:
        return 1
    q, bit = divmod(n, 2)
    return F(q) if bit else F(q) + q + (q & 1)


def h(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for used, counts in enumerate(dp):
        for tail, count in enumerate(counts):
            if not count:
                continue
            avail = rows[tail] & ~used
            while avail:
                bit = avail & -avail
                avail ^= bit
                dp[used | bit][bit.bit_length()-1] += count
    return sum(dp[-1])


def apex(marks):
    n = max(marks)+1
    rows = [sum(1 << j for j in range(i+1, n)) for i in range(n)] + [0]
    for i in range(n):
        if n-i-1 in marks:
            rows[n] |= 1 << i
        else:
            rows[i] |= 1 << n
    return rows


def prime(n):
    return n >= 2 and all(n % i for i in range(2, isqrt(n)+1))


def check_faces():
    tests = 0
    for k in range(2, 6):
        for j in range(1, 12-2*k):
            x, y, z = j+2*k, j+k, k-1
            base = 1 + 2**x
            delta = 2**y - 2**z
            for subset, expected in (
                ({x}, base), ({x, y}, base+delta), ({x, z}, base-delta),
                ({x, y, z}, base-2**j),
            ):
                rows = apex(subset)
                assert h(rows) == expected
                assert len(rows) == j+2*k+2
                tests += 1
    print(f"literal Hamilton-path DP on {tests} face corners, orders 7 through 13")


def check_factor_criterion():
    valid = 0
    for q in range(3, 2002, 2):
        for e in range(1, q+1, 2):
            if q % e:
                continue
            f = q // e
            if F(e) + F(f-1) == e+1:
                assert (e+1)//2 <= f <= 2*e-1
                valid += 1
    for p in range(3, 501, 2):
        if prime(p):
            expected = p == 3 or (p >= 11 and (p-3) & (p-4) == 0)
            assert (F(p) + F(p-1) == p+1) == expected
    print(f"intercept-one factor tests through q=2001: {valid} ordered factorizations")
    print("prime-square classification agrees for all odd primes through 500")


def inflate(rows, m, k):
    labels = [0] * m + [1] * k + list(range(2, len(rows)))
    return [sum(1 << j for j, b in enumerate(labels)
                if (i < j if a == b else rows[a] >> b & 1))
            for i, a in enumerate(labels)]


def check_readers():
    low = [110, 12, 88, 0, 11, 30, 58]
    high = [12, 61, 80, 4, 9, 29, 59]
    assert h([18, 8, 3, 5, 14]) == 13 and h([16, 9, 3, 5, 14]) == 15
    for r in range(4, 8):
        center = (1 << (r-2)) * ((1 << r)-1)
        assert h(inflate(low, r-3, r-2)) == center-1
        assert h(inflate(high, r-3, r-2)) == center+1
    print("two-port half-neighbor DP, order 2r for r=4..7, and both r=3 controls")


def check_full_offset(rows, factors, width_gap, exponent_shift, cycle):
    a, b = 1, 9
    assert all(rows[u] >> v & 1 for u, v in zip(cycle, cycle[1:]))
    product, left, right, constant = factors
    full = (1 << len(rows)) - 1
    counts = Counter()

    def extend(used, tail, p, q):
        if used == full:
            counts[p, q] += 1
        # The ordinary letters are used exactly once. Repeated letters cannot
        # run forever because exactly one of the two arcs between marks exists.
        avail = (full if tail is None else rows[tail]) & (~used | (1 << a) | (1 << b))
        while avail:
            bit = avail & -avail
            avail ^= bit
            v = bit.bit_length()-1
            extend(used | bit, v, p+(v == a), q+(v == b))

    extend(0, None, 0, 0)
    expected = {(1, 1): constant+4*product+2*(left+right),
                (1, 2): right+2*product, (2, 1): left+2*product, (2, 2): product}
    assert counts == expected, counts
    # Relabel the ports to 0,1 for independent ordinary subset path DP.
    perm = [a, b] + [i for i in range(len(rows)) if i not in (a, b)]
    relabelled = [sum(1 << j for j, old in enumerate(perm) if rows[i] >> old & 1)
                  for i in perm]
    for m in (1, 2):
        exponent = m+exponent_shift
        perfect_form = 2**(exponent-1) * (2**exponent-1)
        assert h(inflate(relabelled, m, m+width_gap)) == perfect_form+constant
    print(f"full Euclid offset {constant}: compressed visits {sorted(counts.items())}")


def check_orders():
    for p, d, s in [(3, 2, 1), (5, 4, 1), (31, 5, 1), (257, 16, 1),
                    (1093, 364, 2), (3511, 1755, 2)]:
        assert prime(p) and pow(2, d, p**s) == 1 and pow(2, d, p**(s+1)) != 1
        assert all(pow(2, e, p) != 1 for e in range(1, d))
        for a in range(1, 4):
            period = d * p**max(0, a-s)
            assert pow(2, period, p**a) == 1
            for q in range(2, max(d, p)+1):
                if period % q == 0 and prime(q):
                    assert pow(2, period//q, p**a) != 1
        print(f"prime {p}, first three width periods: {[d * p**max(0,a-s) for a in range(1,4)]}")


if __name__ == "__main__":
    check_faces()
    check_factor_criterion()
    check_readers()
    check_full_offset([998,1020,1016,977,929,456,656,768,576,32],
                      (256,-144,14,-9), 3, 6, (0,1,2,3,7,8,9,5,6,4,0))
    check_full_offset([990,1020,1016,976,928,137,688,768,608,32],
                      (128,-46,38,-5), 0, 4, (0,1,2,3,4,7,8,6,9,5,0))
    check_orders()
