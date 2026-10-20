"""Independent arc-set replay of the complete endpoint-two case partition.

No imported producer, graph-isomorphism library, or assertion-only guard.
Exhaustive at n=9,10; deterministic samples at n=11..100. No strong-partner
filter: the symbolic theorem does not need one.
"""
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
from random import Random


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def middle_orbits(r):
    seen = set()
    forced, free = [], []
    for i, j in combinations(range(1, r + 1), 2):
        if (i, j) in seen:
            continue
        orbit = {(i, j), (r + 1 - j, r + 1 - i)}
        seen.update(orbit)
        (forced if j == i + 1 else free).append(orbit)
    return forced, free


def middle_from_bits(forced, free, bits):
    arcs = set().union(*forced)
    for orbit, bit in zip(free, bits):
        arcs.update(orbit if bit else {(j, i) for i, j in orbit})
    return arcs


def endpoint_source(n, middle, u):
    z = n - 1
    arcs = set(middle) | {(0, i) for i in range(1, z)}
    arcs |= {(i, z) for i in range(1, z)} | {(z, 0)}
    arcs.remove((0, u)); arcs.add((u, 0))
    v = n - 1 - u
    arcs.remove((v, z)); arcs.add((z, v))
    return frozenset(arcs)


def kmove(arcs, word):
    path_pairs = {frozenset((u, v)) for u, v in zip(word, word[1:])}
    return frozenset((u, v) if frozenset((u, v)) in path_pairs else (v, u)
                     for u, v in arcs)


def is_blue(arcs, word):
    rev = dict(zip(word, reversed(word)))
    return all((rev[v], rev[u]) in arcs for u, v in arcs)


def choose_ledger(n, arcs, u):
    a, b, x, y, c, z = 0, 1, 2, n - 3, n - 2, n - 1
    W = tuple(range(2, c))
    D = tuple(range(3, c))
    if u == c:
        if (x, c) in arcs:
            q = (a, x) + D + (c, z, b)
            s = (a, b, x, c) + D[::-1] + (z,)
            v = (b, z) + D + (c, x, a)
            return 'last_port_chord', (q, s[::-1], v)
        return 'last_port_block', (
            (a, x) + D + (c, z, b),
            D + (a, x, c, z, b),
            D + (a, b, x, c, z))
    if (c, b) in arcs:
        if (c, x) in arcs:
            return 'closed_epsilon0', (
                (c,) + W + (b, z, a),
                (z, c) + W + (b, a),
                (a, z, b) + W + (c,))
        if u == x:
            U = tuple(range(3, y))
            return 'closed_second_port', (
                (a, b, x) + U + (z, y, c),
                (a, x) + U + (z, y, b, c),
                (a, x) + U + (y, b, c, z))
        return 'closed_epsilon1', (
            (c, b, z, a) + W,
            (z, c, b, a) + W,
            (a, z, b) + W + (c,))
    if u == y:
        return 'open_penultimate_port', (
            (a, b, c, z, x) + D,
            (b, c, a, z, x) + D,
            (b, x) + D + (c, a, z))
    middle_tail = tuple(range(2, z))
    return 'open_interior_port', (
        (b,) + middle_tail + (z, a),
        (z,) + middle_tail + (b, a),
        (a, z) + middle_tail + (b,))


def replay(n, arcs, u):
    P = tuple(range(n))
    need(len(arcs) == n * (n - 1) // 2, 'complete pair count')
    need(all(((i, j) in arcs) != ((j, i) in arcs)
             for i, j in combinations(P, 2)), 'tournament pairs')
    need(all((i, i + 1) in arcs for i in range(n - 1)), 'source path')
    need(is_blue(arcs, P), 'source blue control')
    need((n - 1, 0) in arcs, 'source closed')
    need(sum(i == 0 for i, _ in arcs) == n - 3, 'first endpoint score')
    need(sum(i == n - 1 for i, _ in arcs) == 2, 'last endpoint score')
    branch, words = choose_ledger(n, arcs, u)
    state = arcs
    for word in words:
        need(len(word) == n and set(word) == set(P), 'Hamiltonian word')
        need(all((i, j) in state for i, j in zip(word, word[1:])),
             f'directed {n=} {u=} {branch=} {word=}')
        need(not is_blue(state, word), f'black {n=} {u=} {branch=} {word=}')
        state = kmove(state, word)
    need(state == kmove(arcs, P), 'literal labelled complement target')
    return branch


def main():
    exhaustive = {}
    for n in (9, 10):
        forced, free = middle_orbits(n - 2)
        counts = Counter()
        for bits in product((0, 1), repeat=len(free)):
            middle = middle_from_bits(forced, free, bits)
            for u in range(2, n - 1):
                counts[replay(n, endpoint_source(n, middle, u), u)] += 1
        exhaustive[n] = dict(counts=counts, total=sum(counts.values()))
    rng = Random(20261020)
    samples = Counter()
    for n in range(11, 101):
        forced, free = middle_orbits(n - 2)
        for _ in range(7):
            bits = [rng.getrandbits(1) for _ in free]
            middle = middle_from_bits(forced, free, bits)
            u = rng.randrange(2, n - 1)
            samples[replay(n, endpoint_source(n, middle, u), u)] += 1
    # This fixed n=8 source is a hostile to using the fixed-tail A word
    # below its proved order. A failure here is not a metagraph bridge.
    rows = (62, 28, 184, 240, 224, 194, 135, 3)
    hostile = frozenset((i, j) for i, row in enumerate(rows)
                       for j in range(8) if (row >> j) & 1)
    branch, words = choose_ledger(8, hostile, 6)
    state, blue_count = hostile, 0
    for word in words:
        need(all((i, j) in state for i, j in zip(word, word[1:])), 'hostile legal')
        blue_count += is_blue(state, word)
        state = kmove(state, word)
    need(branch == 'last_port_block' and blue_count > 0, 'hostile distinguishes scope')
    report = dict(status='PASS', universe='all marked closed-blue endpoint-two presentations at n=9,10',
                  exhaustive=exhaustive,
                  deterministic_orders=[11, 100], sample_counts=dict(samples),
                  sample_total=sum(samples.values()), hostile_n8_blue_moves=blue_count,
                  no_strong_partner_filter=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
