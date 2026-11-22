"""Count unit marked fibers fixed by converse, at a fixed order and mark."""

import argparse
from collections import Counter
from functools import cache
from itertools import combinations, permutations
from math import comb, factorial

from fibers import paths, strong


def partial_orientations(n):
    """All partial orientations reversed by (01)(23)... with {0,1} missing."""
    theta = [v ^ 1 if v < 2 * (n // 2) else v for v in range(n)]
    orbits = []
    for i, j in combinations(range(n), 2):
        mate = (theta[j], theta[i])
        if (i, j) == (0, 1) or (i, j) > tuple(sorted(mate)):
            continue
        orbits.append(((i, j),) if mate == (i, j) else ((i, j), mate))

    rows = [0] * n
    for orbit in orbits:
        for i, j in orbit:
            rows[i] |= 1 << j
    yield tuple(rows)
    for step in range(1, 1 << len(orbits)):
        # Successive Gray codes differ at the lowest set bit of the step.
        switch = (step & -step).bit_length() - 1
        for i, j in orbits[switch]:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
        yield tuple(rows)


@cache
def packing(n):
    # Digit S counts paths using precisely the subset encoded by S.
    width = factorial(n).bit_length()
    mask = (1 << width) - 1
    unused = tuple(sum(mask << (width * s) for s in range(1 << n)
                       if not (s >> v) & 1) for v in range(n))
    shifts = tuple(width << v for v in range(n))
    return unused, shifts, width * ((1 << n) - 1)


def packed_paths(rows):
    n = len(rows)
    unused, shifts, last = packing(n)
    predecessors = [[u for u in range(n) if (rows[u] >> v) & 1] for v in range(n)]
    counts = [1 << shift for shift in shifts]
    for _ in range(1, n):
        counts = [(sum(counts[u] for u in pred) & mask) << shift
                  for pred, mask, shift in zip(predecessors, unused, shifts)]
    return sum(counts) >> last


def class_counts(counter, orbit_size):
    assert all(count % orbit_size == 0 for count in counter.values())
    return Counter({h: count // orbit_size for h, count in counter.items()})


def quotient_path_parities(partial):
    """Parity of all nonempty simple paths by endpoint in the pair tournament."""
    m = len(partial) // 2
    quotient = [0] * m
    for a, b in combinations(range(m), 2):
        contacts = (partial[2 * a] >> (2 * b)) & 3
        if contacts.bit_count() % 2:
            quotient[a] |= 1 << b
        else:
            quotient[b] |= 1 << a
    incoming = [sum(1 << v for v in range(m) if quotient[v] & (1 << u))
                for u in range(m)]
    endpoints = [0] * (1 << m)
    for s in range(1, 1 << m):
        if not s & (s - 1):
            endpoints[s] = s
            continue
        for u in range(m):
            if (s >> u) & 1:
                predecessors = endpoints[s ^ (1 << u)] & incoming[u]
                if predecessors.bit_count() % 2:
                    endpoints[s] |= 1 << u
    total = 0
    for mask in endpoints:
        total ^= mask
    starts = 0
    for order in permutations(range(m)):
        if all(quotient[a] & (1 << b) for a, b in zip(order, order[1:])):
            starts ^= 1 << order[0]
    assert starts == total
    assert total.bit_count() % 2 == 1
    return total


def census(n):
    scalar, marked, central = Counter(), Counter(), Counter()
    m = n // 2
    size = 1 << ((comb(n, 2) + m) // 2 - 1)
    for count, partial in enumerate(partial_orientations(n), 1):
        first, second = list(partial), list(partial)
        first[0] |= 1 << 1
        second[1] |= 1 << 0
        h, k = sorted((packed_paths(first), packed_paths(second)))
        if n <= 5:
            assert [h, k] == sorted((paths(first), paths(second)))
        if n % 2:
            assert (k - h) // 2 % 2 == (quotient_path_parities(partial) & 1)
        else:
            assert h % 4 == k % 4 == 1
        scalar.update((h, k))
        if k == h + 2:
            marked[h] += 1
            # The bidirected completion is strong exactly for central pairs,
            # including the exceptional low heights 1 and 3.
            both = list(first)
            both[1] |= 1 << 0
            if strong(both):
                central[h] += 1
    assert count == size and sum(scalar.values()) == 2 * size
    if n % 2:
        assert sum(count for h, count in scalar.items() if h % 4 == 3) == size
    marked_orbit = (1 << m) * factorial(m - 1)
    return scalar, class_counts(marked, marked_orbit), class_counts(central, marked_orbit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--through", type=int, choices=range(2, 9), default=7)
    parser.add_argument("--height", type=int, default=27)
    args = parser.parse_args()
    print("order  partial_states  all_fixed_units  all_fixed_central  target_fixed  target_fixed_central")
    for n in range(2, args.through + 1):
        scalar, marked, central = census(n)
        states = sum(scalar.values()) // 2
        print(n, states, sum(marked.values()), sum(central.values()),
              marked[args.height], central[args.height], flush=True)
        if n <= 7:
            print("  fixed central units by lower height", sorted(central.items()))
        for order, h in ((4, 5), (5, 11), (5, 13), (6, 17), (8, 257)):
            if n == order:
                m = n // 2
                group = (1 << m) * factorial(m)
                assert scalar[h] % group == 0
                print(f"  rigid scalar height {h}, fixed classes {scalar[h] // group}", flush=True)


if __name__ == "__main__":
    main()
