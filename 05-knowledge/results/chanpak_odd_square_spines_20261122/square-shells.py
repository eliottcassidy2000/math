"""Exact weighted Farey antidiagonals selected by distance in odd-square brackets."""
import argparse
from itertools import permutations
from math import gcd, isqrt
from pathlib import Path
import sys

RESULTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RESULTS / "chanpak_farey_complementary_exact_20261122"))
from cube_controller_exact import TARGET, address_height, chart, envelope


def phase_pairs(data, target):
    """All residue pairs at modulus 2h satisfying the parity-corrected height."""
    h, q, c, t, e, r, k = data
    for x in range(2*h):
        slope = q*x + c
        g = gcd(slope, h)
        step = h // g
        inv = pow(slope // g, -1, step)
        for parity in (0, 1):
            delta = target - k - t*x - slope*parity - r*(x & 1) - e*parity
            if delta % (2*g):
                continue
            first = inv*(delta // (2*g)) % step
            for half in range(first, h, step):
                yield x, 2*half + parity


def sum_bounds(data, target, a, b, distance):
    h, q, c, t, e, r, k = data
    if target < k:
        return 1, 0
    top = isqrt(4*a*b*(target-k) // q) + 4*(distance+1)**2
    linear = max(b*(t+h), a*(c+h))
    constant = 4*a*b*(target-k-r-e)
    lo, hi = 0, top+1
    while lo < hi:
        s = (lo+hi)//2
        if q*s*s + 4*linear*s >= constant:
            hi = s
        else:
            lo = s+1
    return lo, top


def search(data, target, pairs, a, b, distance):
    h, q, c, t, e, r, k = data
    assert a > 0 and b > 0 and gcd(a, b) == 1
    scaled = (h*max(a,b), q, a*c, b*t, a*b*e, a*b*r, a*b*k)
    period = 2*h
    inv = pow(a, -1, b)
    buckets = [[] for _ in range(period)]
    for x, y in pairs:
        weighted_sum = a*x + b*y
        buckets[weighted_sum % period].append((x, weighted_sum))
    bottom, top = sum_bounds(data, target, a, b, distance)
    literal = 0
    hits = []
    for s in range(bottom, top+1):
        residues = buckets[s % period]
        if not residues:
            continue
        for left, right in envelope(scaled, a*b*target, s):
            # The envelope is for U=a*x. Discard the non-integral addresses.
            lo, hi = (left+a-1)//a, right//a
            for x0, phase_sum in residues:
                j = ((s-phase_sum)//period * inv) % b
                first = lo + (x0 + period*j-lo) % (period*b)
                for x in range(first, hi+1, period*b):
                    y = (s-a*x)//b
                    literal += 1
                    if address_height(data, x, y) == target:
                        hits.append((x, y))
    return max(0, top-bottom+1), literal, hits


def controls():
    # Both weights share a prime with the modulus in the 5:17 control.
    data = (85, 4, 12, 8, 4, 2, 3)
    cases = 0
    for u, v in [(0, 0), (1, 9), (6, 3), (17, 10), (48, 41)]:
        target = address_height(data, u, v)
        pairs = list(phase_pairs(data, target))
        h, q, c, t, e, r, k = data
        for a, b in ((1, 1), (2, 3), (5, 17), (17, 5)):
            lo, hi = sum_bounds(data, target, a, b, 2)
            expected = {(x, y) for x in range(1+(target-k)//t)
                        for y in range(1+(target-k-t*x)//(c+q*x))
                        if lo <= a*x+b*y <= hi and address_height(data, x, y) == target}
            _, _, hits = search(data, target, pairs, a, b, 2)
            assert set(hits) == expected and len(hits) == len(expected)
            cases += 1
    print(cases, "weighted-sum intervals compared to complete Cartesian inversion")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=int, default=512)
    parser.add_argument("--optimal-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.distance < 0:
        parser.error("the bracket distance must be nonnegative")
    if args.check:
        controls()
    else:
        for n in (4, 8):
            data = chart(n)
            pairs = list(phase_pairs(data, TARGET))
            directions = [(7, 10)] if args.optimal_only else [(1, 1), *permutations((1, 5, 17), 2), (7, 10)]
            for a, b in directions:
                sums, reads, hits = search(data, TARGET, pairs, a, b, args.distance)
                print(n, a, b, sums, reads, hits, flush=True)
