"""Exact local black-detour controls; research construction, scopes in detours.md."""
from collections import Counter
from itertools import combinations


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def arc(t, a, b):
    return (t[a] >> b) & 1


def path(t, p):
    return len(p) == len(t) and len(set(p)) == len(t) and all(
        arc(t, a, b) for a, b in zip(p, p[1:]))


def blue(t, p):
    n = len(p)
    return all(arc(t, p[i], p[j]) == arc(t, p[n - 1 - j], p[n - 1 - i])
               for i, j in combinations(range(n), 2))


def complement(t, p):
    n = len(t)
    need(path(t, p), "complement requires an actual Hamiltonian path")
    result = [((1 << n) - 1) ^ row ^ (1 << a) for a, row in enumerate(t)]
    for a, b in zip(p, p[1:]):
        result[a] ^= 1 << b
        result[b] ^= 1 << a
    return tuple(result)


def strong(t):
    n = len(t)
    for start in range(n):
        seen = 1 << start
        while True:
            expanded = seen
            for v in range(n):
                if seen >> v & 1:
                    expanded |= t[v]
            if expanded == seen:
                break
            seen = expanded
        if seen != (1 << n) - 1:
            return False
    return True


def isomorphic(t, u):
    n = len(t)
    candidates = {v: tuple(w for w in range(n)
                           if t[v].bit_count() == u[w].bit_count()) for v in range(n)}
    order = sorted(range(n), key=lambda v: len(candidates[v]))
    def extend(mapping, used):
        if len(mapping) == n:
            return True
        v = order[len(mapping)]
        for w in candidates[v]:
            if used >> w & 1:
                continue
            if all(arc(t, v, x) == arc(u, w, y) for x, y in mapping.items()):
                mapping[v] = w
                if extend(mapping, used | (1 << w)):
                    return True
                del mapping[v]
        return False
    return extend({}, 0)


def automorphism(t, permutation):
    return all(arc(t, a, b) == arc(t, permutation[a], permutation[b])
               for a, b in combinations(range(len(t)), 2))


def paths(t, start=None, end=None):
    n = len(t)
    def extend(p, used):
        if len(p) == n:
            if end is None or p[-1] == end:
                yield tuple(p)
            return
        for v in range(n):
            if used >> v & 1 or (end == v and len(p) != n - 1):
                continue
            if not p or arc(t, p[-1], v):
                yield from extend(p + [v], used | (1 << v))
    if start is None:
        yield from extend([], 0)
    else:
        yield from extend([start], 1 << start)


def blue_presentations(n):
    free = {(a, b) for a in range(n) for b in range(a + 2, n)}
    orbits = []
    while free:
        edge = min(free)
        a, b = edge
        orbit = {edge, (n - 1 - b, n - 1 - a)}
        orbits.append(orbit)
        free -= orbit
    for mask in range(1 << len(orbits)):
        rows = [0] * n
        for v in range(n - 1):
            rows[v] |= 1 << (v + 1)
        for i, orbit in enumerate(orbits):
            for a, b in orbit:
                if mask >> i & 1:
                    rows[a] |= 1 << b
                else:
                    rows[b] |= 1 << a
        yield tuple(rows)


def rotation(p, edge):
    for i, pair in enumerate(zip(p, p[1:])):
        if pair == edge:
            return p[i + 1:] + p[:i + 1]
    raise RuntimeError("missing shared cut edge")


def local_hub(t, p):
    partner = complement(t, p)
    candidates = [edge for edge in zip(p, p[1:])
                  if not blue(partner, rotation(p, edge))]
    has_black_endpoint_path = False
    for q in paths(t, p[0], p[-1]):
        if blue(t, q):
            continue
        has_black_endpoint_path = True
        qpartner = complement(t, q)
        qedges = set(zip(q, q[1:]))
        for edge in candidates:
            if edge not in qedges:
                continue
            r = rotation(q, edge)
            if blue(qpartner, r):
                continue
            hub = complement(qpartner, r)
            last = rotation(p, edge)
            need(path(hub, last) and not blue(hub, last), "bad hub last path")
            need(complement(hub, last) == partner, "hub target is not literal")
            return (q, r, last), True
    return None, has_black_endpoint_path


def endpoint_release(t, p):
    n = len(t)
    if n < 5 or arc(t, p[0], p[-1]):
        return None
    a, b, z = p[0], p[1], p[-1]
    w = p[2:-1]
    if not arc(t, b, w[-1]) or not arc(t, w[0], z):
        return None
    words = ((b,) + w + (z, a), (z,) + w + (b, a), (a, z) + w + (b,))
    now = t
    for word in words:
        need(path(now, word), "endpoint-release direction")
        if blue(now, word):
            return None
        now = complement(now, word)
    need(now == complement(t, p), "endpoint-release target")
    return words


def cyclic_wrap(middle):
    r = len(middle)
    rows = [sum(1 << v for v in range(1, r + 1))]
    rows += [(row << 1) | (1 << (r + 1)) for row in middle]
    rows += [1]
    return tuple(rows)


def wrap_words(middle):
    r = len(middle)
    a, z = 0, r + 1
    v = tuple(range(1, r + 1))
    b, c = v[0], v[-1]
    if arc(middle, 0, r - 1):
        w = v[1:]
        return 'open', ((b,) + w + (z, a), (z,) + w + (b, a),
                        (a, z) + w + (b,))
    w = v[1:-1]
    if arc(middle, 0, r - 2):
        return 'closed_chord', ((c, b, z, a) + w, (z, c, b, a) + w,
                                (a, z) + v)
    return 'closed_no_chord', ((c,) + w + (b, z, a), (z, c) + w + (b, a),
                               (a, z) + v)


def wrap_replay(middle):
    t = cyclic_wrap(middle)
    p = tuple(range(len(t)))
    need(blue(t, p), 'wrapped middle is blue')
    target = complement(t, p)
    need(strong(t) and strong(target), 'wrapped endpoints must both be strong')
    branch, words = wrap_words(middle)
    now = t
    charges = []
    for q in words:
        need(path(now, q), 'wrapped direction')
        need(not blue(now, q), 'wrapped path must be black')
        charges.append(now[q[0]].bit_count() + now[q[-1]].bit_count() - len(t) + 1)
        now = complement(now, q)
    need(now == target, 'wrapped literal target')
    r = len(middle)
    db, dc = middle[0].bit_count(), middle[-1].bit_count()
    if branch == 'open':
        expected = (db, -1, db - 1)
    elif branch == 'closed_no_chord':
        expected = (dc, -1, dc - 1)
    else:
        dy = middle[-2].bit_count()
        expected = (dc + dy - r + 1, r - 2 - dy, dc - 1)
    need(tuple(charges) == expected, 'wrapped symbolic endpoint charges')
    return branch, tuple(charges)


def main():
    singleton = cyclic_wrap((0,))
    need(all(blue(singleton, p) for p in paths(singleton)), 'C3 source boundary')
    one_target = complement(singleton, tuple(range(3)))
    need(all(blue(one_target, p) for p in paths(one_target)), 'TT3 target boundary')
    four = cyclic_wrap((2, 0))
    need(isomorphic(four, complement(four, tuple(range(4)))), 'r2 merged loop boundary')
    rotational_hostile = (26, 36, 41, 18, 38, 9)
    need(blue(rotational_hostile, tuple(range(6))) and
         blue(rotational_hostile, (1, 2, 3, 4, 5, 0)), 'even cyclic blue hostile')
    need(not automorphism(rotational_hostile, (1, 2, 3, 4, 5, 0)) and
         automorphism(rotational_hostile, (2, 3, 4, 5, 0, 1)),
         'shift-by-two, not shift-by-one')
    print('Sharp controls: C3/TT3 pure, r2 merged loop, even shift-by-two hostile.', flush=True)
    wrapped = Counter()
    for r in range(3, 9):
        row = Counter()
        for middle in blue_presentations(r):
            if r == 3 and not arc(middle, 0, r - 1):
                target = complement(cyclic_wrap(middle), tuple(range(r + 2)))
                need(all(blue(target, p) for p in paths(target)), 'C3 middle pure-target hostile')
                row['pure_target_hostile'] += 1
                continue
            branch, charges = wrap_replay(middle)
            row[branch] += 1
            wrapped[branch] += 1
        print('wrap', r, dict(row), flush=True)
    print('wrapped totals', dict(wrapped), flush=True)
    for n in range(4, 9):
        counts = Counter()
        hostiles = []
        for t in blue_presentations(n):
            p = tuple(range(n))
            need(blue(t, p), "blue generator")
            if arc(t, p[-1], p[0]):
                shift_two = tuple((v + 2) % n for v in range(n))
                need(blue(t, p[1:] + p[:1]) == automorphism(t, shift_two),
                     'whole-cycle symmetry equivalence')
            if not arc(t, p[0], p[-1]) or not strong(t):
                continue
            partner = complement(t, p)
            if not strong(partner):
                continue
            if n <= 5 and (not any(not blue(t, q) for q in paths(t)) or
                          not any(not blue(partner, q) for q in paths(partner))):
                continue
            counts['strong_mixed_open'] += 1
            if isomorphic(t, partner):
                counts['loop'] += 1
                continue
            released = endpoint_release(partner, p)
            if released:
                counts['endpoint_release'] += 1
            result, has_endpoint_path = local_hub(t, p)
            if result:
                counts['hub'] += 1
            else:
                counts['hub_fail'] += 1
                counts['endpoint_flexible' if has_endpoint_path else 'endpoint_frozen'] += 1
                if released:
                    counts['release_repairs_hub_failure'] += 1
                if len(hostiles) < 4:
                    hostiles.append((t, tuple(x.bit_count() for x in t), has_endpoint_path))
        print(n, dict(counts), 'hostiles', hostiles, flush=True)
    print('PASS: directedness, literal targets, symbolic charges, sharp boundaries and guarded hostiles.')


if __name__ == '__main__':
    main()
