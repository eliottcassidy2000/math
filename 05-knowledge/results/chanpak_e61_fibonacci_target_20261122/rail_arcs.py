"""Exact all-arc derivatives of the three-port strong perfect-number lower construction."""
import argparse
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations
from math import comb, prod
from pathlib import Path
import sys

RESULTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RESULTS / "chanpak_support_birth_grammar_20261122"))
from check_perfect_port import H, LOWER, LPORTS, cap_two, inflate


def endpoint_skeletons(rows, ports, cap):
    """Count compressed nonempty paths by ordinary mask, port visits and last port or vertex."""
    current = Counter()
    for v in range(len(rows)):
        visits = tuple(int(v == p) for p in ports)
        current[(0 if v in ports else 1 << v, visits, v)] = 1
    totals = defaultdict(Counter)
    while current:
        following = Counter()
        for (mask, visits, last), weight in current.items():
            totals[mask, visits][last] += weight
            choices = rows[last] & ~mask
            while choices:
                bit = choices & -choices
                choices ^= bit
                v = bit.bit_length() - 1
                if v in ports:
                    i = ports.index(v)
                    if cap is not None and visits[i] == cap:
                        continue
                    new_visits = tuple(n + int(k == i) for k, n in enumerate(visits))
                    following[mask, new_visits, v] += weight
                else:
                    following[mask | bit, visits, v] += weight
        current = following
    return totals


def two_slot_kernel(rows, ports):
    scores = sorted(sum(rows[p] >> q & 1 for q in ports) for p in ports)
    if scores == list(range(len(ports))):
        cap = None  # Port labels are acyclic between successive singleton visits.
    else:
        cap_two(rows, ports)
        cap = 2
    full = (1 << len(rows)) - 1
    ordinary = full ^ sum(1 << p for p in ports)
    ends = endpoint_skeletons(rows, ports, cap)
    starts = endpoint_skeletons(tuple(full ^ r ^ (1 << i) for i, r in enumerate(rows)), ports, cap)
    by_mask = defaultdict(list)
    for (mask, visits), start_counts in starts.items():
        by_mask[mask].append((visits, start_counts))
    kernel = defaultdict(Counter)
    for (mask, left), end_counts in ends.items():
        for right, start_counts in by_mask[ordinary ^ mask]:
            runs = tuple(a+b for a, b in zip(left, right))
            if min(runs) == 0:
                continue
            if cap is not None:
                assert max(runs) <= 3
            for i, e in end_counts.items():
                for j, s in start_counts.items():
                    kernel[i, j][runs] += e*s
    return kernel


@cache
def fillings(width, runs, last, first):
    """Assign increasing block labels to runs; pin the indicated two-slot boundary."""
    fixed = int(last is not None) + int(first is not None)
    free = runs - fixed
    if free < 0 or (last is not None and last == first):
        return 0
    total = 0
    for omitted in range(free + 1):
        available = free - omitted
        choices = (available + int(last is not None and x < last) + int(first is not None and x > first)
                   for x in range(width) if x != last and x != first)
        total += (-1)**omitted * comb(free, omitted) * prod(choices)
    return total


def arc_changes(rows, ports, widths, kernel):
    physical = [(v, rank) for v in range(len(rows)) for rank in range(dict(zip(ports, widths)).get(v, 1))]
    def joined(i, a, j, b):
        return sum(ways * prod(fillings(w, r, a if p == i else None, b if p == j else None)
                               for p, w, r in zip(ports, widths, runs))
                   for runs, ways in kernel[i, j].items())
    for (i, a), (j, b) in combinations(physical, 2):
        delta = joined(j, b, i, a) - joined(i, a, j, b)
        yield ((i, a), (j, b)), delta if i == j or rows[i] >> j & 1 else -delta


def inspect(exponent, kernel):
    widths = (exponent, (exponent-3)//2, (exponent-3)//2)
    changes = list(arc_changes(LOWER, LPORTS, widths, kernel))
    assert len(changes) == comb(2*exponent+2, 2) and all(d % 2 == 0 for _, d in changes)
    if exponent == 5:
        graph = inflate(LOWER, LPORTS, widths)
        old = H(graph)
        names = [(v, a) for v in range(len(LOWER)) for a in range(dict(zip(LPORTS, widths)).get(v, 1))]
        for (a, b), delta in changes:
            i, j = names.index(a), names.index(b)
            other = list(graph)
            other[i] ^= 1 << j
            other[j] ^= 1 << i
            assert H(other) - old == delta
    positive = min((d, edge) for edge, d in changes if d > 0)
    negative = max((d, edge) for edge, d in changes if d < 0)
    print('exponent', exponent, 'arcs', len(changes), 'distinct changes', len(set(d for _,d in changes)))
    print('neutral', sum(d == 0 for _, d in changes), 'unit increases', sum(d == 2 for _, d in changes))
    print('smallest positive and largest negative', positive, negative)
    return changes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent', type=int, default=61)
    args = parser.parse_args()
    if args.exponent < 5 or args.exponent % 2 == 0:
        parser.error('the strong lower rail uses odd exponents at least five')
    kernel = two_slot_kernel(LOWER, LPORTS)
    inspect(5, kernel)
    if args.exponent != 5:
        inspect(args.exponent, kernel)
