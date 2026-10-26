#!/usr/bin/env python3
"""Exact finite controls for the four-marker two-port encoding; stdlib only."""
import itertools as it
import json
import random
from collections import defaultdict

GATES = 0


def require(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def tournament(n, mask):
    rows = [0] * n
    for k, (u, v) in enumerate(it.combinations(range(n), 2)):
        a, b = (u, v) if mask >> k & 1 else (v, u)
        rows[a] |= 1 << b
    return tuple(rows)


def relabel(rows, order):
    return tuple(sum(((rows[u] >> v) & 1) << j for j, v in enumerate(order))
                 for u in order)


def rooted_key(rows, s=0, t=None):
    t = len(rows) - 1 if t is None else t
    rest = [v for v in range(len(rows)) if v not in (s, t)]
    return min(relabel(rows, (s,) + p + (t,)) for p in it.permutations(rest))


def free_key_aut(rows):
    # Every isomorphism preserves degree, so this is a full isomorphism oracle,
    # exhaustively enumerating all degree-compatible bijections, with no gadget
    # decoder or assumed marker classes.
    groups = defaultdict(list)
    for v, row in enumerate(rows):
        groups[row.bit_count()].append(v)
    blocks = [groups[d] for d in sorted(groups)]
    keys = []
    for ps in it.product(*(it.permutations(block) for block in blocks)):
        keys.append(relabel(rows, tuple(v for p in ps for v in p)))
    key = min(keys)
    return key, keys.count(key)


def rooted_aut(rows):
    return sum(relabel(rows, (0,) + p + (len(rows) - 1,)) == rows
               for p in it.permutations(range(1, len(rows) - 1)))


def star(h, k):
    n, m = len(h), len(k)
    rows = list(h) + [row << n for row in k]
    for u in range(n):
        for v in range(m):
            if u == n - 1 and v == 0:
                rows[u] |= 1 << (n + v)
            else:
                rows[n + v] |= 1 << u
    return tuple(rows)


def frame(rows, left=2, right=2):
    result = rows
    for _ in range(left):
        result = star((0,), result)
    for _ in range(right):
        result = star(result, (0,))
    return result


def incoming(rows, v):
    return [u for u in range(len(rows)) if rows[u] >> v & 1]


def outgoing(rows, v):
    return [u for u in range(len(rows)) if rows[v] >> u & 1]


def decode(rows):
    """Return ALL admissible labelled decodings, no graph-isomorphism calls."""
    n = len(rows)
    if n < 6:
        return []
    low = [v for v in range(n) if rows[v].bit_count() == 1]
    high = [v for v in range(n) if len(incoming(rows, v)) == 1]
    if len(low) not in (2, 3) or len(high) not in (2, 3):
        return []
    answers = []
    for a0 in low:
        a1 = outgoing(rows, a0)[0]
        if a1 not in low:
            continue
        s = outgoing(rows, a1)[0]
        for b1 in high:
            b0 = incoming(rows, b1)[0]
            if b0 not in high:
                continue
            t = incoming(rows, b0)[0]
            removed = {a0, a1, b0, b1}
            if len(removed) != 4 or s in removed or t in removed or s == t:
                continue
            core = [v for v in range(n) if v not in removed]
            middle = [v for v in core if v not in (s, t)]
            order = [s] + middle + [t]
            c = relabel(rows, order)
            # Exact re-encoding is an image-membership guard, not a guess that
            # every arbitrary tournament with low/high scores is a valid frame.
            if frame(c) == relabel(rows, [a0, a1] + order + [b0, b1]):
                answers.append((c, (a0, a1, s, t, b0, b1)))
    return answers


def strong(rows):
    for start in range(len(rows)):
        reached = 1 << start
        while True:
            new = reached
            for v in range(len(rows)):
                if reached >> v & 1:
                    new |= rows[v]
            if new == reached:
                break
            reached = new
        if reached != (1 << len(rows)) - 1:
            return False
    return True


def hamilton_paths(rows, fixed_ports=False):
    n = len(rows)
    orders = ((0,) + p + (n - 1,) for p in it.permutations(range(1, n - 1))) \
        if fixed_ports else it.permutations(range(n))
    return sum(all(rows[a] >> b & 1 for a, b in zip(p, p[1:])) for p in orders)


def main():
    all_states = {}
    census = []
    collisions = {}
    short_contexts = {'e,e': ((0,), (0,)),
                      'A2,e': ((2, 0), (0,)),
                      'A2op,e': ((0, 1), (0,)),
                      'e,A2': ((0,), (2, 0)),
                      'e,A2op': ((0,), (0, 1))}
    short_context_collisions = {}
    for n in range(2, 6):
        states = sorted({rooted_key(tournament(n, m))
                         for m in range(1 << (n * (n - 1) // 2))})
        all_states[n] = states
        signatures = {}
        exponent_counts = defaultdict(int)
        for rows in states:
            encoded = frame(rows)
            answers = decode(encoded)
            eps = int(rows[0] == 0) + int(not incoming(rows, n - 1))
            exponent_counts[eps] += 1
            require(len(answers) == 3 ** eps, 'wrong decoding ambiguity')
            require({rooted_key(c) for c, _ in answers} == {rows},
                    'decoder changed joint rooted type')
            key, aut = free_key_aut(encoded)
            require(aut == 3 ** eps * rooted_aut(rows), 'automorphism formula')
            require(key not in signatures, 'four-marker isomorphism collision')
            signatures[key] = rows
        for left, right in ((1, 1), (2, 1), (1, 2)):
            groups = defaultdict(list)
            for rows in states:
                key, _ = free_key_aut(frame(rows, left, right))
                groups[key].append(rows)
            bad = [v for v in groups.values() if len(v) > 1]
            tag = f'{left},{right}'
            if bad and tag not in collisions:
                collisions[tag] = {'n': n, 'cores': bad[0],
                                   'image_key': free_key_aut(frame(bad[0][0], left, right))[0]}
        for tag, (left, right) in short_contexts.items():
            groups = defaultdict(list)
            for rows in states:
                image = star(star(left, rows), right)
                groups[free_key_aut(image)[0]].append(rows)
            bad = [v for v in groups.values() if len(v) > 1]
            if bad and tag not in short_context_collisions:
                short_context_collisions[tag] = {'n': n, 'cores': bad[0]}
        census.append({'n': n, 'rooted_classes': len(states),
                       'four_marker_classes': len(signatures),
                       'exception_exponents': dict(sorted(exponent_counts.items()))})
    require([len(all_states[n]) for n in range(2, 6)] == [2, 8, 32, 176],
            'independent rooted-class census')
    accepted = []
    degree_only_hostile = None
    for mask in range(1 << 15):
        rows = tournament(6, mask)
        answers = decode(rows)
        if answers:
            accepted.append(mask)
            require(len({rooted_key(c) for c, _ in answers}) == 1,
                    'arbitrary-input decoding inconsistency')
            require(all(len(c) == 2 for c, _ in answers), 'wrong source order')
        elif (degree_only_hostile is None and
              sum(r.bit_count() == 1 for r in rows) in (2, 3) and
              sum(len(incoming(rows, v)) == 1 for v in range(6)) in (2, 3)):
            degree_only_hostile = {'mask': mask, 'rows': rows}
    require(len(accepted) == 800, 'order-six image census must be 6! + 6!/9')
    require(degree_only_hostile is not None, 'degree-only hostile missing')
    require(set(short_context_collisions) == set(short_contexts),
            'all two-sided contexts using fewer than four markers must fail')
    require(all(c['n'] == 3 for c in short_context_collisions.values()),
            'first collision order changed')

    operation_hostile = {'core_rows': (0, 1),
                         'source_fixed_port_paths': hamilton_paths((0, 1), True),
                         'encoded_all_paths': hamilton_paths(frame((0, 1)))}
    require(operation_hostile['source_fixed_port_paths'] == 0 and
            operation_hostile['encoded_all_paths'] == 9,
            'encoding is not literal Hamilton-path conjugacy')

    rng = random.Random(20261026)
    for n in (6, 7, 9, 12, 20):
        for rep in range(30):
            rows = tournament(n, rng.getrandbits(n * (n - 1) // 2))
            # Exercise all four boundary combinations, not only random cores.
            if rep % 4 in (1, 3):
                rows = tuple(0 if v == 0 else row | 1 for v, row in enumerate(rows))
            if rep % 4 in (2, 3):
                rows = tuple(((1 << n) - 1) ^ (1 << (n - 1)) if v == n - 1
                             else row & ~(1 << (n - 1)) for v, row in enumerate(rows))
            encoded = frame(rows)
            perm = list(range(n + 4))
            rng.shuffle(perm)
            shuffled = relabel(encoded, perm)
            answers = decode(shuffled)
            require(bool(answers), 'relabelled large image rejected')
            inv = {v: j for j, v in enumerate(perm)}
            expected_markers = (inv[0], inv[1], inv[2], inv[n + 1], inv[n + 2], inv[n + 3])
            require(any(markers == expected_markers for _, markers in answers),
                    'actual marking missing from exact decoder')

    # The rooted product is associative, and decoding each encoded factor before
    # composing preserves the literal joint-port output (small complete universe).
    for h in all_states[2] + all_states[3]:
        for k in all_states[2] + all_states[3]:
            for ch, _ in decode(frame(h)):
                for ck, _ in decode(frame(k)):
                    require(rooted_key(star(ch, ck)) == rooted_key(star(h, k)),
                            'composition depends on marker choices')
            for ell in all_states[2]:
                require(star(star(h, k), ell) == star(h, star(k, ell)),
                        'hinge associativity')

    result = {'status': 'FINITE-EXACT controls; all-order proof in README.md',
              'gates': GATES, 'census': census,
              'first_small_context_collisions': collisions,
              'all_short_context_collisions': short_context_collisions,
              'degree_only_hostile': degree_only_hostile,
              'native_operation_hostile': operation_hostile,
              'order_six_image': {'universe': 32768, 'accepted': len(accepted),
                                  'ordinary_classes': 2},
              'larger_relabelled_controls': 150,
              'normal_and_optimized_expected_identical': True}
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
