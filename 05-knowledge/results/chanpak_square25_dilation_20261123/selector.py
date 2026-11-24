"""Exact seed data, strong fixed-order construction and absolute-target inverse."""
from functools import cache
from heapq import heappop, heappush
from math import gcd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_tiling_viewpoint_20261123"))
from precedence_bank import SEEDS

OUTER = (2, 4, 9, 35, 15, 23)
DATA = ((5,32,4), (17,224,6), (37,352,6), (77,1904,8), (129,4884,9), (751,17636,9))


def readings(rows):
    full = (1 << len(rows)) - 1

    @cache
    def paths(last, unused):
        if not unused:
            return 1
        choices = unused if last < 0 else unused & rows[last]
        count = 0
        while choices:
            bit = choices & -choices
            choices -= bit
            count += paths(bit.bit_length()-1, unused ^ bit)
        return count

    return paths(-1, full), sum(paths(-1, s)*paths(-1, full ^ s) for s in range(full+1))


def skeleton(rows, ports):
    """Compress every maximal run in a selected module to one port occurrence.

    This enumerator allows every possible number of returns for the given
    templates. Their port-induced graph is acyclic, so recursion terminates.
    """
    full = (1 << len(rows)) - 1
    port_mask = sum(1 << i for i in ports)
    result = {}

    def extend(last, unused, visits):
        if not unused and all(visits):
            result[visits] = result.get(visits, 0) + 1
        choices = unused | port_mask
        if last >= 0:
            choices &= rows[last]
        while choices:
            bit = choices & -choices
            choices ^= bit
            v = bit.bit_length()-1
            counts = tuple(k + (v == p) for p,k in zip(ports, visits))
            extend(v, unused & ~bit, counts)

    extend(-1, full ^ port_mask, (0,) * len(ports))
    return result


def precedence(rows):
    full = (1 << len(rows)) - 1
    pred = [0]*len(rows)

    def extend(last, used):
        if used == full:
            return 1
        choices = (full if last < 0 else rows[last]) & ~used
        count = 0
        while choices:
            bit = choices & -choices
            choices ^= bit
            v = bit.bit_length()-1
            tails = extend(v, used | bit)
            if tails % 2:
                pred[v] ^= used
            count += tails
        return count

    extend(-1, 0)
    return pred


def cycle_parity_at(rows, start):
    def extend(last, used):
        parity = int(used.bit_count() >= 3 and (rows[last] >> start) & 1)
        choices = rows[last] & ~used
        while choices:
            bit = choices & -choices
            choices ^= bit
            parity ^= extend(bit.bit_length() - 1, used | bit)
        return parity

    return extend(start, 1 << start)


def check_substitution(outer, blocks):
    """Compare the four parity rules with literal words of the inflated tournament."""
    positions = [(i, v) for i, block in enumerate(blocks) for v in range(len(block))]
    rows = [0] * len(positions)
    expected = [0] * len(positions)
    exterior = precedence(outer)
    interior = [precedence(block) for block in blocks]
    for x, (i, u) in enumerate(positions):
        for y, (j, v) in enumerate(positions):
            if x == y:
                continue
            if i == j:
                arc = (blocks[i][u] >> v) & 1
                before = (interior[i][u] >> v) & 1
                if len(blocks[i]) == 2:
                    before ^= cycle_parity_at(outer, i)
            else:
                arc = (outer[i] >> j) & 1
                before = 1 - arc
                if len(blocks[i]) == len(blocks[j]) == 1:
                    before = (exterior[i] >> j) & 1
            rows[x] |= arc << y
            expected[x] |= before << y
    assert precedence(rows) == expected
    return rows


def inverse(target):
    candidates = 0
    matches = []

    def branch(i, h, q, exponents):
        nonlocal candidates
        if i == len(DATA):
            candidates += 1
            scale, rem = divmod(target-5*h, 6*q)
            if not rem and scale > 0 and scale & (scale-1) == 0:
                t = scale.bit_length()-1
                if t + sum(e*d for e,(_,_,d) in zip(exponents, DATA)) >= 3:
                    matches.append((exponents, t))
            return
        base, slots, _ = DATA[i]
        e = 0
        while 5*h + 6*q <= target:
            branch(i+1, h, q, exponents + (e,))
            h *= base
            q *= slots
            e += 1

    branch(0, 1, 1, ())
    return candidates, matches


def local_product_group(modulus, ternary):
    # theta is q/(2^d h) after cancelling common powers of three.
    start = (1, 1)
    dist = {start: 0}
    queue = [(0, start)]
    while queue:
        cost, state = heappop(queue)
        if cost != dist[state]:
            continue
        for h,q,d in DATA:
            if gcd(h, modulus) != 1:
                continue
            theta = q * pow((1 << d)*h, -1, ternary) % ternary
            next_state = (state[0]*h % modulus, state[1]*theta % ternary)
            new_cost = cost+d
            if new_cost < dist.get(next_state, float("inf")):
                dist[next_state] = new_cost
                heappush(queue, (new_cost, next_state))
    return dist


def small_wheel():
    costs = {1:0}
    queue = [(0,1)]
    while queue:
        n,s = heappop(queue)
        if n != costs[s]:
            continue
        for h,_,d in DATA:
            v = s*h % 24
            if n+d < costs.get(v,float("inf")):
                costs[v] = n+d
                heappush(queue,(n+d,v))
    assert set(costs) == set(range(1,24,2))
    return {r:5+max(3,costs[5*r%24]) for r in range(1,24,2)}


if __name__ == "__main__":
    one = (0,)
    pair = (2, 0)
    triangle = (2, 4, 1)
    for blocks in [(pair, pair, one), (triangle, one, pair), (triangle, triangle, triangle)]:
        check_substitution(triangle, blocks)
    check_substitution(OUTER, (one, pair, one, one, pair, one))
    for filler in (triangle, SEEDS[5]):
        inflated = check_substitution(OUTER, (filler, one, one, one, one, one))
        h, q = readings(filler)
        assert q % 4 == 0 and readings(inflated)[0] == 5*h + 6*q
    for h,q,d in DATA:
        assert len(SEEDS[h]) == d and readings(SEEDS[h]) == (h,q)
    assert skeleton(OUTER, (0,)) == {(1,):17,(2,):6}
    ext = [5,2,4,3,1]
    pred = precedence(OUTER)
    for i,v in enumerate(ext):
        assert pred[v] & ~1 == sum(1<<w for w in ext[:i])
    assert OUTER[0] == 1<<1
    M = 8*9*11
    group = local_product_group(M, 9)
    assert len(group) == sum(gcd(x,M)==1 for x in range(M))*3
    assert max(group.values()) == 101 and max(small_wheel().values()) == 27
    print("unit/ternary states at modulus", M, len(group), "maximum vertex cost", max(group.values()))
    print("least orders in this wrapper at modulus 24", small_wheel())
    E = 2**60*(2**61-1)
    for target in (E-1, E+1):
        report = inverse(target)
        assert report == (22337, [])
        print("target", target, "exponent vectors and matches", report)
