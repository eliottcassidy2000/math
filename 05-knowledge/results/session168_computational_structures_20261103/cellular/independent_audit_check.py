#!/usr/bin/env python3
"""Independent audit controls. Imports neither the candidate checker nor its output."""
import itertools
import json
import random


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


MASKS = (153, 8345, 8337, 8369, 177, 185)


def physical_round(h, width):
    h = list(h)
    increments = []
    layers = [
        (range(width), lambda s, i: s[i] + s[width+i] + s[width+(i+1)%width]),
        (range(width), lambda s, i: s[i] + s[width+i]),
        (range(width, 2*width), lambda s, i: s[i-width] - s[i]),
        (range(width), lambda s, i: s[i] - s[width+i]),
    ]
    for targets, target_value in layers:
        old = h[:]
        for i in targets:
            d = (target_value(old, i) - old[i]) % 3
            if d == 2:
                d = -1
            need((MASKS[old[i]] ^ MASKS[(old[i]+d)%6]).bit_count() == bool(d), 'literal arc')
            increments.append((i, d))
            h[i] = (old[i] + d) % 6
    return tuple(h), increments


def main():
    pairs = list(itertools.combinations(range(6), 2))
    counts = []
    for mask in MASKS:
        arcs = {(a,b) if mask >> j & 1 else (b,a) for j,(a,b) in enumerate(pairs)}
        count = sum(all((u,v) in arcs for u,v in zip(p,p[1:])) for p in itertools.permutations(range(6)))
        need(count == 29, 'independent Hamilton path enumeration')
        counts.append(count)
    flips = []
    for i in range(6):
        difference = MASKS[i] ^ MASKS[(i+1)%6]
        need(difference.bit_count() == 1, 'seed arc difference')
        flips.append(pairs[difference.bit_length()-1])
    need(flips == [(3,5),(0,4),(1,2)]*2, 'literal arc sequence')
    expected = {1:8,2:8,3:24,4:40,5:328,7:2920,9:72,27:216,81:648}
    periods = {}
    for width, target in expected.items():
        initial = (1,) + (0,)*(2*width-1)
        x = initial
        for time in range(1, target+1):
            a,b = x[:width],x[width:]
            x = b + tuple((a[i]+b[i]+b[(i+1)%width])%3 for i in range(width))
            if x == initial:
                break
        need(x == initial and time == target, 'logical first return')
        periods[width] = time
    rng = random.Random(168)
    cancellation = []
    for width in (1,3,9,27,81,243):
        for _ in range(4):
            initial = tuple(rng.randrange(6) for _ in range(2*width))
            h = initial
            first, second = [], []
            for time in range(8*width):
                h, increments = physical_round(h, width)
                (first if time < 4*width else second).extend(increments)
                if time == 4*width-1:
                    need(tuple(v%3 for v in h) == tuple(-v%3 for v in initial), 'negative half-period')
            need(h == initial, 'physical whole-state return')
            need(second == [(i,-d) for i,d in first], 'termwise signed-increment cancellation')
        cancellation.append({'width':width,'random_physical_starts':4})
    print(json.dumps({'hamilton_counts':counts,'literal_arc_sequence':flips,
                      'logical_first_returns':periods,'cancellation_controls':cancellation},sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
