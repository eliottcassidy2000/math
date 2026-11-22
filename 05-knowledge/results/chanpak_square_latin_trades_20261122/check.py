"""Exact checks of the two-matching theorem and the five-square family."""

from itertools import combinations
from math import comb, gcd, isqrt


def residue_pairs(a, g):
    return [(i, (a - i) % g) for i in range(g) if i < (a - i) % g]


def component_sizes(b, a, c, pairs, chosen):
    g = 2 * len(pairs)
    target = [a] * g
    for index in chosen:
        i, j = pairs[index]
        target[i] = target[j] = c

    # Walk actual ordinary edges; label 0 represents vertex b.
    seen = bytearray(b)
    sizes = []
    for start in range(b):
        if seen[start]:
            continue
        seen[start] = 1
        stack = [start]
        size = 0
        while stack:
            x = stack.pop()
            size += 1
            for y in ((-x) % b, (target[x % g] - x) % b):
                if x != y:
                    assert (x or b) + (y or b) in (b, a, a+b, c, c+b)
                if not seen[y]:
                    seen[y] = 1
                    stack.append(y)
        sizes.append(size)
    return sorted(sizes)


def predicted_sizes(b, a, c, g, ell):
    d = gcd(a + 2 * ell * ((c - a) // g), b // g)
    return [b // d] + [2 * b // d] * ((d - 1) // 2)


def direct_path_count(b, a, c):
    # Enumerate paths using ordinary sums, with b fixed as the first endpoint.
    allowed = {b, a, c, b+a, b+c}
    neighbors = {x: [y for y in range(1, b+1) if y != x and x+y in allowed]
                 for x in range(1, b+1)}
    required = {tuple(sorted((x, b-x))) for x in range(1, b) if x != b-x}

    def extend(path, seen):
        if path[-1] == b // 2:
            if len(path) != b:
                return 0
            edges = {tuple(sorted(pair)) for pair in zip(path, path[1:])}
            return int(required <= edges)
        return sum(extend(path + [y], seen | {y})
                   for y in neighbors[path[-1]] if y not in seen)

    return extend([b], {b})


def check_small():
    states = 0
    systems = 0
    for b in range(4, 33, 2):
        for a, c in combinations([x for x in range(b) if gcd(x, b) == 1], 2):
            g = gcd(b, c-a)
            pairs = residue_pairs(a, g)
            assert len(pairs) == g // 2
            good = 0
            for ell in range(len(pairs) + 1):
                expected = predicted_sizes(b, a, c, g, ell)
                good += (expected == [b]) * comb(len(pairs), ell)
                for chosen in combinations(range(len(pairs)), ell):
                    assert component_sizes(b, a, c, pairs, chosen) == expected
                    states += 1
            if b <= 12:
                assert direct_path_count(b, a, c) == good
            systems += 1
    return systems, states


def check_example():
    b, a, c = 840**2, 41**2, 559**2
    assert all(isqrt(x)**2 == x for x in (b, a, c, a+b, c+b))
    g = gcd(b, c-a)
    pairs = residue_pairs(a, g)
    for ell in (0, 1, 2100, 2106):
        actual = component_sizes(b, a, c, pairs, range(ell))
        assert actual == predicted_sizes(b, a, c, g, ell)
        print(f"840, layer {ell}: component sizes {actual}")
    assert comb(4200, 2100) > 2**4187


def check_family():
    for k in range(4, 103, 2):
        m, n = 2*k+1, k
        u, v, s = m*m-n*n, 2*m*n, m*m+n*n
        r = u+v
        b = 2*u*v*r*s
        a, c = sorted((((u*s)**2 - (v*r)**2)**2, ((v*s)**2 - (u*r)**2)**2))
        assert 0 < a < c < b**2
        assert all(isqrt(x)**2 == x for x in (a+b**2, c+b**2))
        g = gcd(b**2, c-a)
        assert g == 4*u*v*r*s*s and b**2 // g == u*v*r
        assert c-a == abs(v**4-u**4)*(r**4-s**4)
        assert gcd((a+c)//2, b**2//g) == 1
    print("parameter identities: all 50 even k from 4 through 102")


def check_isolated_layers():
    for j in range(10, 40):
        n = 42*j+22
        m = 2*n+21
        u, v, s = m*m-n*n, 2*m*n, m*m+n*n
        r = u+v
        a, c = sorted((((u*s)**2 - (v*r)**2)**2, ((v*s)**2 - (u*r)**2)**2))
        b = 2*u*v*r*s
        g = gcd(b*b, c-a)
        modulus = b*b // g
        delta = (c-a) // g
        assert a < c < b*b and modulus % 21 == 0
        assert a % 21 == 1 and delta % 21 == 37 % 21
        # Construct the CRT offset without factoring modulus: remove all 3,7 parts.
        other = modulus
        for prime in (3, 7):
            while other % prime == 0:
                other //= prime
        q = (other * ((6 * pow(other, -1, 21)) % 21)) % modulus
        ell = g//4+q
        assert gcd(a + 2*ell*delta, modulus) == 1
        assert (a + 2*(ell-1)*delta) % 7 == 0
        assert (a + 2*(ell+1)*delta) % 3 == 0
    print("isolated-layer CRT test: j=10 through 39")


if __name__ == "__main__":
    systems, states = check_small()
    print(f"even B through 32: {systems} unit pairs, {states} trade choices")
    check_example()
    check_family()
    check_isolated_layers()
