#!/usr/bin/env python3
"""Independent Session243 audit: no producer or inherited-verifier imports.

Rows are outgoing-neighbour masks. Counts are exact integers. The primary
method is a memoized first-vertex recurrence. The second uses literal
permutations and every split. Slot coefficients supply an independent source
response, without the rooted U statistic or producer transfer formula.
"""
from functools import lru_cache
from itertools import permutations, combinations
from hashlib import sha256
import json
import sys


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def vertices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def validate(rows):
    n = len(rows)
    require(all(0 <= x < 1 << n and not (x >> i & 1)
                for i, x in enumerate(rows)), 'invalid row')
    require(all(((rows[i] >> j) ^ (rows[j] >> i)) & 1
                for i in range(n) for j in range(i)), 'not a tournament')


def decode(n, code):
    rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if code & 1 else (j, i)
            rows[u] |= 1 << v
            code >>= 1
    return tuple(rows)


def transitive(n):
    return tuple(((1 << n) - 1) ^ ((1 << (i + 1)) - 1) for i in range(n))


def table(rows):
    n = len(rows)
    @lru_cache(None)
    def first(mask, v):
        tail = mask ^ (1 << v)
        return 1 if not tail else sum(first(tail, u) for u in vertices(rows[v] & tail))
    h = [1] + [sum(first(mask, v) for v in vertices(mask))
               for mask in range(1, 1 << n)]
    return h, first


def gadget(rows):
    h, _ = table(rows)
    full = (1 << len(rows)) - 1
    return h[full], sum(h[s] * h[full ^ s] for s in range(full + 1))


def profile(rows, root):
    h, first = table(rows)
    n = len(rows)
    full = (1 << n) - 1
    reverse = tuple(sum(1 << j for j in range(n) if rows[j] >> i & 1) for i in range(n))
    _, last = table(reverse)
    return (h[full], last(full, root),
            sum(first(s, root) * h[full ^ s] for s in range(full + 1) if s >> root & 1),
            h[full ^ (1 << root)])


def literal_profile(rows, root):
    """No dynamic programming: enumerate permutations of every induced subset."""
    n = len(rows)
    full = (1 << n) - 1
    h, starts = [0] * (full + 1), [0] * (full + 1)
    h[0] = 1
    terminal = 0
    for size in range(1, n + 1):
        for p in permutations(range(n), size):
            if all(rows[u] >> v & 1 for u, v in zip(p, p[1:])):
                mask = sum(1 << v for v in p)
                h[mask] += 1
                starts[mask] += p[0] == root
                terminal += size == n and p[-1] == root
    return (h[full], terminal, sum(starts[s] * h[full ^ s] for s in range(full + 1)),
            h[full ^ (1 << root)])


def coefficients(rows, root):
    """Count eligible B-run slots in permutations of A; c[k] is k runs."""
    n = len(rows)
    c = [0] * (n + 2)
    for p in permutations(range(n)):
        eligible = set()
        if p[0] == root:
            eligible.add(0)
        if p[-1] != root:
            eligible.add(n)
        for i in range(1, n):
            if p[i - 1] != root and p[i] == root:
                eligible.add(i)
        mandatory = {i for i in range(1, n) if not rows[p[i - 1]] >> p[i] & 1}
        if not mandatory <= eligible:
            continue
        optional = sorted(eligible - mandatory)
        for k in range(len(optional) + 1):
            for _ in combinations(optional, k):
                c[len(mandatory) + k] += 1
    require(not any(c[3:]), 'single root allowed more than two runs')
    return c[1], c[2]


def ab(state):
    dtotal, end, cover, deleted = state
    return dtotal + 2 * deleted - end - cover, cover - deleted


def graft(rows, root, block):
    n, m = len(rows), len(block)
    ans = list(rows) + [(x << n) | (1 << root) for x in block]
    fresh = ((1 << m) - 1) << n
    for v in range(n):
        if v != root:
            ans[v] |= fresh
    validate(ans)
    return tuple(ans)


def ordinal(left, right):
    n, m = len(left), len(right)
    fresh = ((1 << m) - 1) << n
    return tuple(x | fresh for x in left) + tuple(x << n for x in right)


def reachable(rows, root):
    seen = 1 << root
    while True:
        after = seen
        for v in vertices(seen):
            after |= rows[v]
        if after == seen:
            return seen == (1 << len(rows)) - 1
        seen = after


def strong(rows):
    return bool(rows) and all(reachable(rows, v) for v in range(len(rows)))


def one_check(rows, root, block):
    state = profile(rows, root)
    a, b = ab(state)
    H, Q = gadget(block)
    physical = graft(rows, root, block)
    got = profile(physical, root)
    want = (a * H + b * Q, state[3] * H, b * Q + state[3] * H, state[3] * H)
    require(got == want, ('full state', rows, root, block, got, want))
    require(ab(got) == (H * a, Q * b), 'diagonal transform')
    require(strong(physical) == (len(rows) >= 2 and reachable(rows, root)), 'strong iff')
    return got


def main():
    if '--fail' in sys.argv:
        require(False, 'intentional failure gate')
    seed_data = [((8, 1, 3, 6), 3, -1, (5, 1, 7, 1)),
                 ((8, 1, 3, 6), 2,  1, (5, 1, 9, 3)),
                 ((20, 1, 2, 7, 14), 1, 3, (9, 3, 9, 3)),
                 ((20, 1, 2, 7, 14), 2, 5, (9, 3, 11, 5))]
    report = {'status': 'PASS', 'independent_methods': ['first-vertex subset recurrence',
              'literal permutations on all subsets', 'mandatory/optional B-run slots'], 'seeds': []}
    for rows, root, a, expected in seed_data:
        validate(rows)
        p = profile(rows, root)
        require(p == expected == literal_profile(rows, root), 'source state')
        require(ab(p) == (a, 6), 'source diagonal coordinates')
        c1, c2 = coefficients(rows, root)
        require((c1 - 2 * c2, c2) == (a, 6), 'independent source-slot coefficients')
        require(strong(rows), 'seed strong')
        report['seeds'].append({'rows': rows, 'root': root, 'state': p,
                                'ab': [a, 6], 'run_coefficients': [c1, c2]})

    # Full small universe including singleton, sink and nonstrong source roots.
    small_blocks = [decode(n, k) for n in range(1, 4) for k in range(1 << (n * (n - 1) // 2))]
    tests = 0
    for n in range(1, 5):
        for k in range(1 << (n * (n - 1) // 2)):
            rows = decode(n, k)
            for root in range(n):
                require(profile(rows, root) == literal_profile(rows, root), 'small source literal')
                c1, c2 = coefficients(rows, root)
                a, b = ab(profile(rows, root))
                require((c1 - 2 * c2, c2) == (a, b), 'all small slots')
                for block in small_blocks:
                    one_check(rows, root, block)
                    tests += 1
    report['all_source_order_le4_all_roots_gadget_order_le3'] = tests

    # Four fixed sources against every labelled gadget through order five.
    digest = sha256()
    tests = 0
    gadgets = 0
    for n in range(1, 6):
        for k in range(1 << (n * (n - 1) // 2)):
            block = decode(n, k)
            H, Q = gadget(block)
            require(Q >= (n + 1) * H and Q % 2 == 0, 'cover floor/parity')
            if n >= 2:
                require(Q % 4 == 0, 'four-divisibility')
            values = []
            for rows, root, a, _ in seed_data:
                result = one_check(rows, root, block)
                require(result[0] == 6 * Q + a * H, 'bank scalar')
                values.append(result[0])
                digest.update(repr((rows, root, block, result)).encode())
                tests += 1
            require(all(y - x == 2 * H for x, y in zip(values, values[1:])), 'common gadget gaps')
            gadgets += 1
    report['bank_all_labelled_gadgets_through5'] = {'gadgets': gadgets, 'grafts': tests,
                                                     'receipt_sha256': digest.hexdigest()}

    report['transitive_controls'] = []
    for m in range(1, 9):
        block = transitive(m)
        require(gadget(block) == (1, 1 << m), 'transitive gadget pair')
        values, sizes = [], []
        for rows, root, a, _ in seed_data:
            result = one_check(rows, root, block)
            require(result[0] == 6 * (1 << m) + a, 'window')
            values.append(result[0]); sizes.append(len(rows) + m)
        report['transitive_controls'].append({'m': m, 'H_values': values, 'orders': sizes})

    # Independent permutation output checks; maximal order seven, no transfer row used.
    literal_checks = 0
    for rows, root, _, _ in seed_data:
        for block in [(0,), (2, 0)]:
            output = graft(rows, root, block)
            require(profile(output, root) == literal_profile(output, root), 'literal output')
            literal_checks += 1
    report['literal_output_permutation_checks'] = literal_checks

    # Nested fixed-root graft equals one graft with an ordinal sum, arc for arc.
    seq_blocks = [(0,), (2, 4, 1), (2, 0)]
    report['nested_controls'] = []
    for rows, root, a0, state0 in seed_data:
        out = rows
        joined = ()
        hprod, qprod = 1, 1
        for block in seq_blocks:
            out = graft(out, root, block)
            joined = ordinal(joined, block)
            H, Q = gadget(block)
            hprod *= H; qprod *= Q
            require(out == graft(rows, root, joined), 'literal ordinal fusion')
            p = profile(out, root)
            require(ab(p) == (a0 * hprod, 6 * qprod), 'product law')
            require(p[0] == a0 * hprod + 6 * qprod, 'product scalar')
        report['nested_controls'].append({'root': root, 'seed': rows, 'state': p,
                                         'H_product': hprod, 'Q_product': qprod})

    # Find an accidental Euclid output with intercept !=1; finite witness only.
    hostile = None
    for n in range(2, 6):
        if hostile: break
        for k in range(1 << (n * (n - 1) // 2)):
            if hostile: break
            rows = decode(n, k)
            for root in range(n):
                a, b = ab(profile(rows, root))
                if a == 1: continue
                for m in range(1, 6):
                    if a + b * (1 << m) == 31:
                        out = graft(rows, root, transitive(m))
                        require(profile(out, root)[0] == 31, 'Euclid coincidence literal')
                        hostile = {'source_rows': rows, 'root': root, 'source_ab': [a, b],
                                   'transitive_gadget_order': m, 'actual_H': 31,
                                   'genuine_intercept': a, 'primorial': 30,
                                   'strong': strong(out)}
                        break
                if hostile: break
    require(hostile is not None, 'need actual accidental Euclid witness')
    report['accidental_Euclid_output_not_unit_intercept'] = hostile

    # Deliberate boundaries: empty gadget, root movement, and sink padding.
    rows, root = seed_data[1][:2]
    p = profile(rows, root)
    require(p[0] != sum(ab(p)), 'empty gadget hostile failed')
    out = graft(rows, root, (0,))
    require(ab(profile(out, root)) != ab(profile(out, len(rows))), 'reroot hostile failed')
    padded = ordinal(out, (0,))
    require(gadget(padded)[0] == gadget(out)[0] and not strong(padded), 'padding hostile')
    report['boundary_hostiles'] = {'empty_gadget': {'actual': p[0], 'false_a_plus_b': sum(ab(p))},
        'reroot': {'retained_root_ab': ab(profile(out, root)), 'new_root_ab': ab(profile(out, len(rows)))},
        'sink_padding': {'H_before': gadget(out)[0], 'H_after': gadget(padded)[0], 'strong_after': False}}

    # The unit-intercept argument uses H=1 iff transitive. Exhaustively verify
    # that equivalence through order five, plus the exact v2 gate arithmetically.
    count = 0
    for n in range(1, 6):
        for k in range(1 << (n * (n - 1) // 2)):
            rows = decode(n, k)
            H, Q = gadget(rows)
            no_triangle = all(not (rows[i] >> j & 1 and rows[j] >> l & 1 and rows[l] >> i & 1)
                              for i, j, l in permutations(range(n), 3))
            require((H == 1) == no_triangle, 'unit iff transitive')
            if H == 1:
                require(Q == 1 << n, 'unit response')
            count += 1
    report['unit_gadget_controls'] = count
    print(json.dumps(report, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
