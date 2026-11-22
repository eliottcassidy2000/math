"""Replay finite physical readers and their exact exponent-61 width obstruction."""
import argparse
from pathlib import Path

from check import packing_slots


def literal_subpaths(rows):
    n = len(rows)
    dp = [[0]*n for _ in range(1 << n)]
    paths = [1]*(1 << n)
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1 << v)
                dp[mask][v] = sum(dp[rest][u] for u in range(n) if rows[u] >> v & 1) if rest else 1
        paths[mask] = sum(dp[mask])
    return paths


def literal_slots(rows):
    paths = literal_subpaths(rows)
    layer = paths
    result = [paths[-1]]
    for j in range(2,5):
        following = layer.copy()
        for mask in range(1,len(layer)):
            sub = mask
            while sub:
                following[mask] += paths[sub] * layer[mask ^ sub]
                sub = (sub-1) & mask
        layer = following
        result.append(layer[-1])
    return result


def delayed_witness():
    old = (62,508,472,496,480,452,129,257,65)
    hosts = [tuple(v | (int(not(mask >> i & 1)) << 9) for i,v in enumerate(old)) + (mask,)
             for mask in (17,81)]
    assert [literal_slots(g) for g in hosts] == [[2637,75216,928341,7121472],
                                                  [2637,75216,927621,7115712]]
    for g in hosts:
        for v in range(len(g)):
            seen = 1 << v
            while True:
                more = seen
                for i,row in enumerate(g):
                    if seen >> i & 1:
                        more |= row
                if more == seen:
                    break
                seen = more
            assert seen == (1 << len(g))-1
    for expected in [[77853,77853], [2230125,2228685]]:
        n = len(hosts[0])
        hosts = [tuple(r | (1 << n) for r in g) + (1 << (n+1), (1 << n)-1) for g in hosts]
        assert [literal_subpaths(g)[-1] for g in hosts] == expected
    print('10-vertex physical witness: neutral at 2637 and 77853, then 2230125 and 2228685.')


def run(last):
    delayed_witness()
    source = Path(__file__).resolve().parents[3] / 'session212/cycle-levels'
    expected = {5: (8,8,8,8), 6: (24,24,24,24), 7: (156,156,139,156),
                8: (1728,1727,1442,1726)}
    for n in range(5,last+1):
        profiles = set()
        for line in (source / f'classes{n}.rows').read_text().splitlines():
            p = packing_slots(tuple(map(int, line.split())))
            profiles.add(tuple(p[k] for k in range(n+1)))
        pairs = {(sum(p), sum(a << k for k,a in enumerate(p))) for p in profiles}
        def count(m):
            return len({2*h + (2**m-1)*(q-h) for h,q in pairs})
        assert (len(profiles),len(pairs),count(1),count(2)) == expected[n]
        if n == 8:
            assert count(3) == len(pairs)
        print(n, 'profiles', len(profiles), 'H values', len({h for h,q in pairs}),
              'H/Q2 pairs', len(pairs), 'widths 1/2/3', [count(m) for m in (1,2,3)])
    e = 2**60 * (2**61-1)
    assert 9*(2**59-1)**2 > e+1 and 40320 < 2**59-1
    assert 64*(2**58-1)**2 > e+1


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--order-eight', action='store_true')
    run(8 if cli.parse_args().order_eight else 7)
