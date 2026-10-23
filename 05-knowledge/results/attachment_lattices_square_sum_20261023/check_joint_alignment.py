"""Exact controls for the simultaneous D4 decoder and odd chart character guard."""
from collections import defaultdict
from itertools import product
from math import gcd
from pathlib import Path
import json
import sys

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


GROUP = tuple(product((False, True), (-1, 1), (-1, 1)))


def act(g, z):
    swap, a, b = g
    x, y = z[::-1] if swap else z
    return a*x, b*y


def canonical(z, w):
    return min((act(g, z), act(g, w)) for g in GROUP)


def invariants(z, w):
    x, y = z
    u, v = w
    return (x*x+y*y, x*x*y*y, u*u+v*v, u*u*v*v,
            x*u+y*v, x*x*x*u+y*y*y*v)


def pair_orbit_count(z, w):
    zs = {act(g, z) for g in GROUP}
    ws = {act(g, w) for g in GROUP}
    return len({canonical(a, b) for a in zs for b in ws})


def main():
    points = tuple(product(range(-4, 5), repeat=2))
    by_invariant = {}
    by_orbit = {}
    truncated = defaultdict(set)
    for z, w in product(points, repeat=2):
        inv = invariants(z, w)
        orbit = canonical(z, w)
        need(by_invariant.setdefault(inv, orbit) == orbit, 'separation')
        need(by_orbit.setdefault(orbit, inv) == inv, 'orbit invariance')
        for g in GROUP:
            need(invariants(act(g, z), act(g, w)) == inv, 'literal invariance')
        truncated[inv[:5]].add(orbit)
    expected = (
        ('generic/generic', (2, 1), (3, 1), 8),
        ('axis/axis', (1, 0), (2, 0), 3),
        ('axis/diagonal', (1, 0), (1, 1), 2),
        ('diagonal/diagonal', (1, 1), (2, 2), 3),
        ('origin/generic', (0, 0), (2, 1), 1),
    )
    counts = {}
    for name, z, w, target in expected:
        value = pair_orbit_count(z, w)
        need(value == target, name)
        counts[name] = value
    z, w, wp = (2, 1), (-2, 4), (2, -4)
    need(invariants(z, w)[:5] == (5, 4, 20, 64, 0), 'specified hostile')
    need(invariants(z, w)[:5] == invariants(z, wp)[:5], 'truncation collision')
    need(invariants(z, w)[5] == -12 and invariants(z, wp)[5] == 12, 'mixed repair')
    sums = [tuple(a+b for a, b in zip(z, other)) for other in (w, wp)]
    need(sums == [(0, 5), (4, -3)], 'additive witness')
    need(canonical((0, 0), sums[0]) != canonical((0, 0), sums[1]), 'sum orbit loss')
    for m in range(1, 100, 2):
        need(gcd(m, 8) == 1, 'odd group has no 2-primary image')
        need([a for a in range(8) if m*a % 8 == 0] == [0], 'cyclic character')
    result = dict(status='PASS', grid=[-4, 4], ordered_pairs=len(points)**2,
                  simultaneous_orbits=len(by_orbit),
                  five_invariant_collisions=sum(len(v)>1 for v in truncated.values()),
                  relative_frame_counts=counts,
                  hostile=dict(z=z, w=w, alternative=wp, sums=sums))
    frozen = Path(__file__).with_name('joint_alignment_receipt.json')
    data = json.dumps(result, indent=2) + '\n'
    if '--write' in sys.argv:
        frozen.write_text(data, encoding='utf-8')
    else:
        need(json.loads(frozen.read_text(encoding='utf-8')) == json.loads(data), 'frozen receipt')
    print(json.dumps(dict(**result, checks=CHECKS)))


if __name__ == '__main__':
    main()
