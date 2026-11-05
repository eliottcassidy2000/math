#!/usr/bin/env python3
"""Exact controls for module-hull/ear recursion. Standard library only."""
from itertools import combinations, permutations
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def vertices(mask):
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask -= b


def tournament(n, code):
    a = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if (code >> k) & 1:
            a[i] |= 1 << j
        else:
            a[j] |= 1 << i
    return a


def induced(a, mask):
    labels = list(vertices(mask))
    return [sum(1 << j for j, y in enumerate(labels) if (a[x] >> y) & 1)
            for x in labels]


def converse(a):
    full = (1 << len(a)) - 1
    return [full ^ row ^ (1 << i) for i, row in enumerate(a)]


def endpoint_tables(a):
    """Terminal removal; all induced H values and ending vectors."""
    n = len(a)
    E = [[0] * n for _ in range(1 << n)]
    H = [0] * (1 << n)
    H[0] = 1
    for mask in range(1, 1 << n):
        for j in vertices(mask):
            previous = mask ^ (1 << j)
            if not previous:
                E[mask][j] = 1
            else:
                E[mask][j] = sum(E[previous][i] for i in vertices(previous)
                                 if (a[i] >> j) & 1)
        H[mask] = sum(E[mask])
    return H, E


def f2_table(H):
    F = [0] * len(H)
    for mask in range(1, len(H)):
        part = (mask - 1) & mask
        while part:
            F[mask] += H[part] * H[mask ^ part]
            part = (part - 1) & mask
    return F


def literal_count(a):
    return sum(all((a[p[i]] >> p[i + 1]) & 1 for i in range(len(a) - 1))
               for p in permutations(range(len(a))))


def is_strong(a, mask):
    if not mask:
        return False
    root = mask & -mask
    for rows in (a, converse(a)):
        reached = root
        while True:
            new = reached
            for i in vertices(reached):
                new |= rows[i] & mask
            if new == reached:
                break
            reached = new
        if reached != mask:
            return False
    return True


def is_module(a, mask):
    full = (1 << len(a)) - 1
    return all((a[z] & mask) in (0, mask) for z in vertices(full ^ mask))


def mixed(a, mask):
    full = (1 << len(a)) - 1
    return [z for z in vertices(full ^ mask) if (a[z] & mask) not in (0, mask)]


def close(a, start, reverse=False):
    current = start
    chain = [current]
    while True:
        candidates = mixed(a, current)
        if not candidates:
            return current, chain
        z = candidates[-1] if reverse else candidates[0]
        current |= 1 << z
        chain.append(current)


def grow_ears(a, start):
    full = (1 << len(a)) - 1
    current = start
    steps = []
    while current != full:
        candidates = mixed(a, current)
        if candidates:
            after = current | (1 << candidates[0])
            steps.append((current, after, 1))
        else:
            require(is_module(a, current), 'stuck core not a module')
            A = [x for x in vertices(full ^ current) if not (a[x] & current)]
            B = [x for x in vertices(full ^ current) if (a[x] & current) == current]
            pairs = [(x, y) for x in A for y in B if (a[x] >> y) & 1]
            require(bool(pairs), 'no escaping ear in strong ambient')
            x, y = pairs[0]
            after = current | (1 << x) | (1 << y)
            steps.append((current, after, 2))
        require(is_strong(a, after), 'ear step lost strongness')
        current = after
    return steps


def add_vertex(a, X):
    n = len(a)
    return [row | ((1 << n) if not ((X >> i) & 1) else 0)
            for i, row in enumerate(a)] + [X]


def cyclic_pair(a):
    n = len(a)
    return [row | (1 << n) for row in a] + [1 << (n + 1), (1 << n) - 1]


def response_data(a, H, E):
    n = len(a)
    _, starts = endpoint_tables(converse(a))
    full = (1 << n) - 1
    B = [[0] * n for _ in a]
    for mask in range(1, full):
        complement = full ^ mask
        for y in vertices(mask):
            for x in vertices(complement):
                B[y][x] += E[mask][y] * starts[complement][x]
    return starts[full], E[full], B


def response(s, e, B, X):
    n = len(s)
    return (sum(s[x] for x in range(n) if (X >> x) & 1)
            + sum(e[y] for y in range(n) if not ((X >> y) & 1))
            + sum(B[y][x] for y in range(n) if not ((X >> y) & 1)
                  for x in range(n) if (X >> x) & 1))


def main():
    counts = dict(tournaments=0, nonempty_subsets=0, induced_monotonicity_pairs=0,
                  strong_core_pairs=0, strict_threshold_pairs=0, ternary_pairs=0,
                  singleton_steps=0, pair_steps=0, literal_counts=0,
                  cyclic_controls=0, cyclic_literal_controls=0, response_checks=0)
    census = []
    for n in range(1, 7):
        level = dict(order=n, tournaments=0, strong_tournaments=0)
        full = (1 << n) - 1
        for code in range(1 << (n * (n - 1) // 2)):
            a = tournament(n, code)
            counts['tournaments'] += 1
            level['tournaments'] += 1
            H, E = endpoint_tables(a)
            F = f2_table(H)
            strong = [False] + [is_strong(a, m) for m in range(1, full + 1)]
            modules = [m for m in range(1, full + 1) if is_module(a, m)]
            if strong[full]:
                level['strong_tournaments'] += 1
            for core in range(1, full + 1):
                counts['nonempty_subsets'] += 1
                require(F[core] >= (core.bit_count() - 1) * H[core], 'split bound')
                require(F[core] >= (1 << core.bit_count()) - 2, 'subset bound')
                # Every ordered subset pair is covered independently of strongness.
                sub = core
                while sub:
                    require(H[sub] <= H[core] and F[sub] <= F[core], 'induced monotonicity')
                    counts['induced_monotonicity_pairs'] += 1
                    sub = (sub - 1) & core
                small, chain = close(a, core)
                large, _ = close(a, core, True)
                expected = full
                for mod in modules:
                    if mod & core == core:
                        expected &= mod
                require(small == large == expected, 'canonical module hull mismatch')
                require(is_module(a, small), 'hull not module')
                if strong[core]:
                    require(all(strong[m] for m in chain), 'singleton closure lost strongness')
                if not (strong[core] and strong[full]):
                    continue
                counts['strong_core_pairs'] += 1
                if H[full] < 3 * H[core] + F[core]:
                    counts['strict_threshold_pairs'] += 1
                    require(small == full, 'strict-threshold singleton chain absent')
                if H[core] > 1 and H[full] in (3 * H[core] + 2, 3 * H[core] + 4):
                    counts['ternary_pairs'] += 1
                    require(small == full, 'ternary ancestry failure')
                r = 0
                lower = H[core]
                for before, after, size in grow_ears(a, core):
                    if size == 1:
                        counts['singleton_steps'] += 1
                        require(H[before] <= H[after], 'singleton H decreased')
                    else:
                        counts['pair_steps'] += 1
                        require(H[after] == 3 * H[before] + F[before], 'pair count identity')
                        lower *= core.bit_count() + 2 * r + 2
                        r += 1
                require(H[full] >= lower, 'pair-depth lower bound')
            if n <= 5:
                require(literal_count(a) == H[full], 'literal count mismatch')
                counts['literal_counts'] += 1
                C = cyclic_pair(a)
                HC, _ = endpoint_tables(C)
                require(HC[-1] == 3 * H[full] + F[full], 'cyclic formula')
                require(is_strong(C, (1 << (n + 2)) - 1), 'cyclic strongness')
                require(not is_strong(C, full | (1 << n)), 'cyclic first singleton strong')
                require(not is_strong(C, full | (1 << (n + 1))), 'cyclic second singleton strong')
                counts['cyclic_controls'] += 1
                if n <= 4:
                    require(literal_count(C) == HC[-1], 'cyclic literal mismatch')
                    counts['cyclic_literal_controls'] += 1
                s, e, B = response_data(a, H, E)
                require(sum(map(sum, B)) == F[full], 'B total mismatch')
                for X in range(1 << n):
                    child = add_vertex(a, X)
                    actual, _ = endpoint_tables(child)
                    require(response(s, e, B, X) == actual[-1], 'singleton response mismatch')
                    if strong[full]:
                        require(is_strong(child, (1 << (n + 1)) - 1) == (0 < X < full),
                                'mixed iff strong mismatch')
                    counts['response_checks'] += 1
        census.append(level)
    A = [16, 1, 3, 7, 14]
    B = [20, 1, 2, 7, 14]
    hostile = []
    for a in (A, B):
        H, E = endpoint_tables(a)
        F = f2_table(H)
        require(is_strong(a, 31) and H[-1] == 9 and F[-1] == 66, 'scalar hostile metadata')
        spectrum = sorted({endpoint_tables(add_vertex(a, X))[0][-1] for X in range(1, 31)})
        hostile.append(dict(rows=a, H=H[-1], F2=F[-1], strong_child_spectrum=spectrum))
    require(hostile[0]['strong_child_spectrum'] == [15,17,19,23,25,27,29,33,37,41], 'hostile A')
    require(hostile[1]['strong_child_spectrum'] == [15,17,19,23,29,31,33,37,43], 'hostile B')
    # At equality, a C3 core in the regular T5 can still have full module hull.
    regular = [sum(1 << ((i+d) % 5) for d in (1,2)) for i in range(5)]
    HR, _ = endpoint_tables(regular)
    FR = f2_table(HR)
    core = next(m for m in range(1,32) if m.bit_count() == 3 and is_strong(regular,m))
    require(HR[-1] == 15 == 3*HR[core]+FR[core], 'equality positive count')
    require(close(regular,core)[0] == 31, 'equality positive chain')
    proper_hull_hostile = [22,20,24,19,32,15]
    seed = (1 << 0) | (1 << 2) | (1 << 3)
    require(is_strong(proper_hull_hostile, 63), 'proper-hull ambient strongness')
    require(is_strong(proper_hull_hostile, seed), 'proper-hull core strongness')
    require(not is_module(proper_hull_hostile, seed), 'proper-hull core unexpectedly module')
    require(close(proper_hull_hostile, seed)[0] == 15, 'proper-hull closure')
    require(endpoint_tables(proper_hull_hostile)[0][-1] == 37, 'proper-hull H')
    # Singleton h>1 caveat: H5 is a 3*1+2 target, but no singleton strong step exists.
    H5 = [2,12,9,1]
    require(endpoint_tables(H5)[0][-1] == 5 and is_strong(H5,15), 'H5 boundary')
    require(close(H5,1)[0] == 1, 'singleton boundary')
    receipt = dict(status='PASS', scope='All labelled tournaments of orders 1 through 6; '
                   'all nonempty subsets; every one-vertex response from parents through order 5.',
                   counts=counts, census=census, scalar_compression_hostile=hostile,
                   equality_with_growth=dict(rows=regular,core_mask=core,H=15),
                   proper_hull_hostile=dict(rows=proper_hull_hostile,core_mask=seed,hull_mask=15,H=37),
                   all_order_status='Proved in PROOF.md, not inferred from finite checks.',
                   spectrum_status='OPEN')
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    Path(__file__).with_name('receipt.json').write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
