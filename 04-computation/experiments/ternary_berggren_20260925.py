"""Exact Berggren recursion, consecutive-edge antichains, and ray clocks.

Standard library only. All checks remain active under python -O.
"""

from fractions import Fraction
from math import gcd, isqrt


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def valuation(n, prime):
    check(n != 0, "valuation domain")
    n = abs(n)
    count = 0
    while n % prime == 0:
        n //= prime
        count += 1
    return count


def odd_step(n, sign=1):
    value = 3 * n + sign
    k = valuation(value, 2)
    return value // 2 ** k, k


def root_pair(x, y):
    check(x != y and x > 0 and y > 0 and x % 2 and y % 2, "edge domain")
    check(gcd(x, y) == 1, "primitive edge")
    return max(x, y), min(x, y)


def triple(pair):
    s, t = pair
    return s * t, (s * s - t * t) // 2, (s * s + t * t) // 2


def triple_inverse(value):
    a, b, c = value
    s, t = isqrt(c + b), isqrt(c - b)
    check(s * s == c + b and t * t == c - b and s * t == a, "triple inverse")
    return s, t


def child(pair, branch):
    s, t = pair
    return ((s + 2 * t, t), (2 * s + t, s), (2 * s - t, s))[branch]


def parent(pair):
    s, t = pair
    check(pair != (3, 1), "root has no positive parent")
    if s > 3 * t:
        return s - 2 * t, t
    if s > 2 * t:
        return t, s - 2 * t
    return t, 2 * t - s


def ancestor(first, second):
    """Inclusive ancestry via the root-pair inverse cones."""
    while second[0] >= first[0] and second[1] >= first[1]:
        if second == first:
            return True
        if second == (3, 1):
            return False
        second = parent(second)
    return False


TRIPLE_MATRICES = (
    ((1, -2, 2), (2, -1, 2), (2, -2, 3)),
    ((1, 2, 2), (2, 1, 2), (2, 2, 3)),
    ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3)),
)
LORENTZ = (1, 1, -1)
TRIPLE_INVERSES = tuple(tuple(tuple(LORENTZ[i] * m[j][i] * LORENTZ[j]
                                       for j in range(3)) for i in range(3))
                          for m in TRIPLE_MATRICES)


def matvec(matrix, vector):
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def triple_ancestor(first, second):
    """Independent ancestry via inverse Lorentz matrices on triple vertices."""
    while second[2] >= first[2]:
        if second == first:
            return True
        candidates = [matvec(m, second) for m in TRIPLE_INVERSES]
        candidates = [p for p in candidates if min(p) > 0 and p[2] < second[2]]
        if second == (3, 4, 5):
            check(not candidates, "triple root boundary")
            return False
        check(len(candidates) == 1, "unique positive triple parent")
        second = candidates[0]
    return False


def cw_pair(k):
    a, b = 1, 1
    for bit in bin(k)[3:]:
        a, b = (a, a + b) if bit == "0" else (a + b, b)
    return a, b


def cw_index(a, b):
    bits = []
    while a != b:
        if a < b:
            bits.append(0)
            b -= a
        else:
            bits.append(1)
            a -= b
    check(a == 1, "Calkin-Wilf primitive domain")
    k = 1
    for bit in reversed(bits):
        k = 2 * k + bit
    return k


def reflect_index(k):
    return 3 * 2 ** (k.bit_length() - 1) - 1 - k


def cw_children(k):
    return 2 * k - 1, 4 * reflect_index(k) + 3, 4 * k + 3


def cw_parent(k):
    check(k != 3, "ordinal root boundary")
    if k % 4 == 1:
        return (k + 1) // 2
    if k % 8 == 3:
        return reflect_index((k - 3) // 4)
    check(k % 8 == 7, "ordinal inverse cone")
    return (k - 3) // 4


def matmul(a, b):
    return tuple(tuple(sum(a[i][j] * b[j][k] for j in range(2))
                       for k in range(2)) for i in range(2))


def projective_invariant(matrix):
    a, b = matrix[0]
    c, d = matrix[1]
    return Fraction((a + d) ** 2, a * d - b * c)


def ray_height(k0, j):
    return 2 ** (k0 - 1) * (4 ** j - 1) // 3


def ray_index(k0, height):
    scale = 2 ** (k0 - 1)
    if height < 0 or height % scale:
        return None
    power = 1 + 3 * (height // scale)
    if power & (power - 1):
        return None
    exponent = power.bit_length() - 1
    return exponent // 2 if exponent % 2 == 0 else None


def ternary_phase(k0, height, depth):
    """Inverse in Z/3^depth, by three-way refinement rather than orbit search."""
    phase, modulus = 0, 1
    for _ in range(depth):
        candidates = [phase + digit * modulus for digit in range(3)]
        modulus *= 3
        matches = [j for j in candidates if (ray_height(k0, j) - height) % modulus == 0]
        check(len(matches) == 1, "unique ternary inverse digit")
        phase = matches[0]
    return phase


def main():
    # An intrinsic tree enumeration checked through two independent charts.
    level = [((3, 1), 3)]
    nodes = 0
    for depth in range(9):
        next_level = []
        for pair, ordinal in level:
            s, t = pair
            check(s > t > 0 and s % 2 and t % 2 and gcd(s, t) == 1, "root-pair domain")
            value = triple(pair)
            check(value[0] ** 2 + value[1] ** 2 == value[2] ** 2, "Pythagorean")
            check(gcd(gcd(value[0], value[1]), value[2]) == 1, "primitive triple")
            check(triple_inverse(value) == pair, "lossless root chart")
            m, n = (s + t) // 2, (s - t) // 2
            check(cw_pair(ordinal) == (m, n) and cw_index(m, n) == ordinal, "CW conjugacy")
            if depth < 8:
                for branch, k in enumerate(cw_children(ordinal)):
                    out = child(pair, branch)
                    check(parent(out) == pair and cw_parent(k) == ordinal, "both inverses")
                    check(matvec(TRIPLE_MATRICES[branch], value) == triple(out), "matrix independent child")
                    next_level.append((out, k))
            nodes += 1
        level = next_level
    print("TREE_CHARTS", nodes, "all Berggren words depth0..8 PASS")

    ordinals = 0
    for k in range(1, 65536):
        m, n = cw_pair(k)
        good = m > n and (m - n) % 2 == 1
        check(good == (k % 6 in (3, 5)), "exact ordinal domain")
        if good:
            check(cw_index(m, n) == k, "ordinal inverse")
            ordinals += 1
    print("ORDINAL_DOMAIN", ordinals, "k1..65535, exactly k=3,5 mod6 PASS")

    # New scope: arbitrary ancestor distance, not just immediate adjacency.
    kinds = {"rise-rise": 0, "rise-fall": 0, "fall-rise": 0, "fall-fall": 0}
    independent = 0
    for x in range(3, 100002, 2):
        y, _ = odd_step(x)
        if y == 1:
            continue
        z, _ = odd_step(y)
        first, second = root_pair(x, y), root_pair(y, z)
        check(not ancestor(first, second) and not ancestor(second, first), "consecutive plus antichain")
        kind = ("rise" if y > x else "fall") + "-" + ("rise" if z > y else "fall")
        kinds[kind] += 1
        if x <= 10001:
            check(not triple_ancestor(triple(first), triple(second))
                  and not triple_ancestor(triple(second), triple(first)), "independent triple antichain")
            independent += 1
    print("CONSECUTIVE_PLUS", sum(kinds.values()), "odd x3..100001, nondegenerate", kinds)
    print("INDEPENDENT_TRIPLE_ANCESTRY", independent, "odd x<=10001 PASS")
    y, k = odd_step(27, -1)
    z, j = odd_step(y, -1)
    check((y, z, k, j) == (5, 7, 4, 1), "minus hostile orbit")
    check(child(child((7, 5), 0), 0) == (27, 5), "minus strict ancestry")
    check(triple_ancestor(triple((7, 5)), triple((27, 5))), "minus independent ancestry")
    check(odd_step(7, -1)[0] == 5, "minus unoriented two-cycle")
    print("MINUS_HOSTILES 27->5->7: pair(27,5)=B1^2(7,5); 5<->7 shares one triangle PASS")

    matrices = (((1, 2), (0, 1)), ((2, 1), (1, 0)), ((2, -1), (1, 0)))
    level = [((1, 0), (0, 1))]
    matrix_count = 0
    for depth in range(9):
        for matrix in level:
            check(projective_invariant(matrix).denominator == 1, "unimodular invariant")
            matrix_count += 1
        if depth < 8:
            level = [matmul(m, old) for old in level for m in matrices]
    affine_count = 0
    for r in range(1, 17):
        for total in range(r, 65):
            for sign in (1, -1):
                matrix = ((3 ** r, sign), (0, 2 ** total))
                invariant = projective_invariant(matrix)
                check(invariant.denominator == 3 ** r * 2 ** total, "coprime invariant denominator")
                affine_count += 1
    check(projective_invariant(((4, 1), (0, 1))) == Fraction(25, 4), "sibling invariant")
    print("PROJECTIVE_OBSTRUCTION", matrix_count, "unimodular words;", affine_count,
          "affine slope controls r1..16,K=r..64,both signs PASS")

    # Exact inverse of the inherited lacunary ray, and its 3-adic sidecar.
    residues = comparisons = 0
    for k0 in range(2, 10):
        for depth in range(1, 7):
            modulus = 3 ** depth
            images = {ray_height(k0, j) % modulus for j in range(modulus)}
            check(len(images) == modulus, "ray ternary residue bijection")
            residues += modulus
        for j in range(81):
            h = ray_height(k0, j)
            check(ray_index(k0, h) == j, "exact ray inverse")
            check(ray_height(k0, j + 1) == 4 * h + 2 ** (k0 - 1), "adding-machine conjugacy")
            for i in range(j):
                check(valuation(h - ray_height(k0, i), 3) == valuation(j - i, 3), "ray 3-adic isometry")
                comparisons += 1
        for h in range(2001):
            idx = ray_index(k0, h)
            check((idx is not None) == (h in {ray_height(k0, j) for j in range(8)}), "hostile ray heights")
        for h in range(243):
            phase = ternary_phase(k0, h, 5)
            check((ray_height(k0, phase) - h) % 243 == 0, "constructive ternary inverse")
            next_phase = ternary_phase(k0, 4 * h + 2 ** (k0 - 1), 5)
            check((next_phase - phase) % 243 == 1, "finite adding-machine square")
    print("RAY_TERNARY", residues, "complete residues, k0=2..9, depth1..6;", comparisons, "valuation pairs PASS")
    print("TERNARY_INVERSE 1944 height classes at depth5, k0=2..9; inverse digits and commuting squares PASS")

    returns = 0
    for sign in (1, -1):
        for y in range(1, 102, 2):
            for k in range(2, 9):
                numerator = 2 ** k * y - sign
                if numerator % 3:
                    continue
                x = numerator // 3
                if x <= y:
                    continue
                start = root_pair(x, y)
                gap = 2 ** (k - 1)
                for h in range(1, gap + 1):
                    candidate = x + 2 * h * y
                    quotient, remainder = divmod(3 * candidate + sign, y)
                    legal = not remainder and quotient > 0 and quotient & (quotient - 1) == 0
                    check(legal == (h == gap), "exact first return on ray")
                next_x = 4 * x + sign
                check((next_x, y) == (start[0] + 2 * gap * start[1], start[1]), "ray induced sibling")
                returns += 1
    print("RAY_FIRST_RETURN", returns, "all signed valid fibres y<=101,k2..8; every intermediate height checked PASS")
    print("ALL CHECKS PASS; no Collatz convergence implication")


if __name__ == "__main__":
    main()
