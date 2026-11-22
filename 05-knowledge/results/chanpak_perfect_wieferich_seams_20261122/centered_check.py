"""Test the centered inverse against direct contact words and then at higher exponents."""

from itertools import combinations
from math import isqrt

from check import apex, h


def cost(v):
    marks = [i+1 for i in range(v.bit_length()) if v >> i & 1]
    return sum((1 << (b-a-1)) + (b-a == 1) for a, b in combinations(marks, 2))


def direct_small():
    tested = 0
    for p in (3, 5, 7):
        center, top = p-1, 2*p-1
        for v in range(1 << (p-2)):
            for e in (0, 1):
                marks = {top}
                for a in range(1, p-1):
                    if v >> (a-1) & 1:
                        marks.update([center-a, center+a])
                if e:
                    marks.add(center)
                # Weight the actual inversion pairs, including all adjacent positions.
                raw = 1
                for t in marks:
                    for u in range(t):
                        if u not in marks:
                            raw += (1 << (t-u-1)) + (t-u == 1)
                expr = 1 + (1 << top) - 2*(v*v+cost(v)+e*(v+(v % 2)))
                assert raw == expr
                if p <= 5:
                    assert h(apex(marks)) == expr
                tested += 1
    print(f"all {tested} centered masks at p=3,5,7, including literal DP through p=5")


def check_digits(limit):
    # alpha=1/sqrt(2). Track floor(2**d*alpha), its squared remainder, and pair cost.
    v, remainder, weight = 1, 1, 0
    hits = []
    for d in range(1, limit+1):
        if d <= 100:
            x = 1 << (2*d-1)
            assert v == isqrt(x) and remainder == x-v*v and weight == cost(v)
        digit = int(remainder > v)
        next_weight = weight + digit*(v+(v % 2))
        solutions = [e for e in (0, 1) if remainder == weight+e*(v+(v % 2))]
        if d >= 3:
            assert solutions == ([digit] if remainder == next_weight else [])
        if solutions:
            hits.append((d, solutions))
        remainder = 4*(remainder-digit*v)-digit
        v = 2*v+digit
        weight = next_weight
    print(f"centered solutions through d={limit} (odd exponents through {2*limit+1}): {hits}")


if __name__ == "__main__":
    direct_small()
    check_digits(100_000)
