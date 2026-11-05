#!/usr/bin/env python3
"""Standalone exact controls for the intact-module threshold.
All proof obligations remain enabled under python -O. No inherited bank.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def labelled(n):
    pairs = list(itertools.combinations(range(n), 2))
    for word in range(1 << len(pairs)):
        rows = [0] * n
        for k, (i, j) in enumerate(pairs):
            u, v = (i, j) if word >> k & 1 else (j, i)
            rows[u] |= 1 << v
        yield tuple(rows)


def all_h(rows):
    """Endpoint subset DP returns H for all induced vertex sets."""
    n = len(rows)
    d = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for v in range(n):
        d[1 << v][v] = 1
    for mask in range(1, 1 << n):
        h[mask] = sum(d[mask])
        for u in range(n):
            if not d[mask][u]:
                continue
            todo = rows[u] & ~mask
            while todo:
                bit = todo & -todo
                v = bit.bit_length() - 1
                d[mask | bit][v] += d[mask][u]
                todo -= bit
    return h


def literal_h(rows):
    return sum(all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in itertools.permutations(range(len(rows))))


def substitute(q, s):
    """Expand quotient vertex zero into S; all other vertices stay singleton."""
    n = len(s)
    out = list(s) + [0] * (len(q) - 1)
    for i in range(1, len(q)):
        y = n + i - 1
        if q[0] >> i & 1:
            for u in range(n):
                out[u] |= 1 << y
        else:
            out[y] |= (1 << n) - 1
        for j in range(1, len(q)):
            if q[i] >> j & 1:
                out[y] |= 1 << (n + j - 1)
    return tuple(out)


def reachable(rows, source):
    seen = 1 << source
    while True:
        more = seen
        for u, row in enumerate(rows):
            if seen >> u & 1:
                more |= row
        if more == seen:
            return seen
        seen = more


def cyclic_at_zero(rows):
    reach = [reachable(rows, i) for i in range(len(rows))]
    return any(reach[0] >> i & 1 and reach[i] & 1
               for i in range(1, len(rows)))


def bin_value(h, n, d):
    total = 0
    for assignment in itertools.product(range(d), repeat=n):
        masks = [0] * d
        for v, slot in enumerate(assignment):
            masks[slot] |= 1 << v
        total += math.prod(h[mask] for mask in masks)
    return total


def ordered_join(a, b):
    n = len(a)
    suffix = ((1 << len(b)) - 1) << n
    return tuple(row | suffix for row in a) + tuple(row << n for row in b)


def main():
    small = [s for n in range(1, 5) for s in labelled(n)]
    D = (2, 12, 9, 17, 7)
    TT5 = tuple(sum(1 << j for j in range(i + 1, 5)) for i in range(5))
    inputs = small + [D, TT5]
    quotients = [q for n in range(1, 5) for q in labelled(n)]
    triangle = (2, 4, 1)
    counts = dict(input_tournaments=len(inputs), labelled_inputs_through4=len(small),
                  quotient_tournaments=len(quotients), contexts=0,
                  cyclic_contexts=0, acyclic_contexts=0,
                  nontransitive_context_exclusions=0, literal_cyclic_checks=0)
    profiles = []
    for s in inputs:
        n = len(s)
        hsub = all_h(s)
        h = hsub[-1]
        full = (1 << n) - 1
        f2 = sum(hsub[mask] * hsub[full ^ mask] for mask in range(1, full))
        require(h % 2 == 1, 'Hamilton-path parity control')
        require(f2 >= max((n - 1) * h, (1 << n) - 2), 'two-cover lower bounds')
        cs = substitute(triangle, s)
        exact = 3 * h + f2
        require(literal_h(cs) == all_h(cs)[-1] == exact, 'independent cyclic identity')
        require(all(reachable(cs, i) == (1 << len(cs)) - 1
                    for i in range(len(cs))), 'cyclic extension strongness')
        counts['literal_cyclic_checks'] += 1
        profiles.append([list(s), h, f2, exact])
        for q in quotients:
            out = substitute(q, s)
            actual = all_h(out)[-1]
            cyclic = cyclic_at_zero(q)
            triangle_exists = any(q[0] >> a & 1 and q[a] >> b & 1 and q[b] & 1
                                  for a in range(1, len(q))
                                  for b in range(1, len(q)) if a != b)
            require(cyclic == triangle_exists, 'root SCC / triangle equivalence')
            if cyclic:
                require(actual >= exact, 'sharp cyclic threshold')
                counts['cyclic_contexts'] += 1
            else:
                require(actual % h == 0, 'singleton quotient SCC factor')
                counts['acyclic_contexts'] += 1
            if h > 1:
                require(actual not in (3 * h + 2, 3 * h + 4), 'literal ternary exclusion')
                counts['nontransitive_context_exclusions'] += 1
            counts['contexts'] += 1
    hd = all_h(D)
    bins = [bin_value(hd, len(D), d) for d in (1, 2, 3)]
    require(bins == [13, 104, 471], 'primitive H13 profile')
    require(math.gcd(*bins) == 1, 'primitive H13 content')
    strong5 = (2, 12, 9, 1)
    require(all_h(strong5)[-1] == 5, 'H5 seed')
    exception = ordered_join(TT5, strong5)
    require(all_h(exception)[-1] == 5, 'H1 exception')
    require(all_h(substitute(triangle, (2, 0)))[-1] == 5, 'two-vertex equality')
    receipt = dict(status='PASS', finite_universe=counts,
                   primitive_hostile=dict(rows=D, bins123=bins, F2=78, threshold=117,
                                          blocked_children=[41, 43]),
                   transitive_exception=dict(rows=exception, H=5, retained_module_order=5),
                   exact_input_profiles=profiles,
                   scope='Finite controls only; all-order theorem proved in PROOF.md')
    raw = (json.dumps(receipt, indent=2, sort_keys=True) + '\n').encode()
    Path(__file__).with_name('receipt.json').write_bytes(raw)
    print(json.dumps(dict(status='PASS', **counts,
                         receipt_sha256=hashlib.sha256(raw).hexdigest()), sort_keys=True))


if __name__ == '__main__':
    main()
