"""Check the one-hole lift certificates using exact finite arithmetic."""
from itertools import combinations, product


def compose(first, second):
    return tuple(first[x] for x in second)


def group_closure(generators):
    identity = tuple(range(len(generators[0])))
    found, new = {identity}, [identity]
    while new:
        word = new.pop()
        for step in generators:
            image = compose(step, word)
            if image not in found:
                found.add(image)
                new.append(image)
    return found


def disjoint(first, second):
    return all(a != b for a, b in zip(first, second))


def check_seven_point_groups():
    linear_group = []
    for a, b, c in product(range(1, 8), repeat=3):
        images = tuple((a if x & 1 else 0) ^ (b if x & 2 else 0) ^
                       (c if x & 4 else 0) for x in range(1, 8))
        if len(set(images)) == 7 and 0 not in images:
            linear_group.append(tuple(x - 1 for x in images))
    identity = tuple(range(7))
    derangements = [g for g in linear_group if disjoint(g, identity)]

    def cliques(chosen, candidates):
        if len(chosen) == 6:
            yield chosen
            return
        for idx, g in enumerate(candidates):
            rest = [h for h in candidates[idx+1:] if disjoint(g, h)]
            if len(chosen) + 1 + len(rest) >= 6:
                yield from cliques(chosen + [g], rest)

    count = 0
    for six in cliques([], derangements):
        sharp = {identity, *six}
        assert all(group_closure([g]) == sharp for g in six)
        count += 1
    assert (len(linear_group), len(derangements), count) == (168, 48, 8)
    print("seven-point linear group: 168 elements, 48 derangements, 8 sharp Singer groups")


def trim(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def add(a, b, sign=1):
    size = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) + sign *
                  (b[i] if i < len(b) else 0)) % 5 for i in range(size)])


def square(a):
    result = [0] * (2*len(a) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(a):
            result[i+j] = (result[i+j] + x*y) % 5
    return trim(result)


def remainder(a, b):
    a = a.copy()
    while len(a) >= len(b):
        shift = len(a) - len(b)
        scalar = a[-1] * pow(b[-1], -1, 5) % 5
        for i, x in enumerate(b):
            a[i+shift] = (a[i+shift] - scalar*x) % 5
        trim(a)
    return a


def check_orbits():
    x = [0, 1]
    iterate, total = x, []
    for _ in range(7):
        total = add(total, iterate)
        iterate = add(add([1], iterate, -1), square(iterate), -1)
    a, b = add(iterate, x, -1), total
    assert (len(a)-1, len(b)-1) == (128, 64)
    while b:
        a, b = b, remainder(a, b)
    assert len(a) == 1
    print("period-seven and zero-trace polynomials mod 5: degrees 128,64, gcd 1")


def cover_check(p, maps, expected_order):
    identity = tuple(range(p))
    labels = range(p)
    directed = {(0, i): identity for i in range(1, p)}
    directed.update(maps)
    for (i, j), perm in list(directed.items()):
        directed[j, i] = tuple(perm.index(v) for v in labels)
    for i, j in combinations(labels, 2):
        routes = [directed[i, j]]
        for k in labels:
            if k != i and k != j:
                routes.append(compose(directed[k, j], directed[i, k]))
        for images in zip(*routes):
            assert len(set(images)) == p - 1
    assert len(group_closure(list(maps.values()))) == expected_order
    print(f"{p} branches, {p} sheets, girth at least five, holonomy order {expected_order}")


def check_petersen_puncture():
    # Petersen's vertices are the two-subsets of {0,1,2,3,4}, adjacent if disjoint.
    points = set(combinations(range(5), 2))
    points.remove((0, 1))
    path = [min(points)]

    def extend():
        if len(path) == len(points):
            return set(path[-1]).isdisjoint(path[0])
        for other in sorted(points.difference(path)):
            if set(path[-1]).isdisjoint(other):
                path.append(other)
                if extend():
                    return True
                path.pop()
        return False

    assert extend()
    print("Petersen puncture: a 9-cycle", path)


if __name__ == "__main__":
    check_orbits()
    check_petersen_puncture()
    check_seven_point_groups()
    cover_check(3, {(1, 2): (1, 2, 0)}, 3)
    cover_check(5, {
        (1, 2): (1, 2, 3, 4, 0), (1, 3): (2, 0, 4, 1, 3),
        (1, 4): (4, 3, 0, 2, 1), (2, 3): (1, 4, 3, 0, 2),
        (2, 4): (3, 2, 4, 1, 0), (3, 4): (2, 4, 1, 0, 3),
    }, 60)
