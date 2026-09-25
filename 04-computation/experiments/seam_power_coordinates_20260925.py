"""Exact controls for scale/power transport and the surviving odd layers.

Standard library. Universes and independent controls are printed in main.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign(x):
    return (x > 0) - (x < 0)


def ratio(a, b):
    return Fraction(a**b, b**a)


def primitive_arrow(a, b):
    g = gcd(a, b)
    p, q = a // g, b // g
    return p**q > q**p


def cycle_counts(tournament, length):
    n = len(tournament)
    cycles = []
    for vertices in combinations(range(n), length):
        first, rest = vertices[0], vertices[1:]
        for suffix in permutations(rest):
            cycle = (first,) + suffix
            if all(tournament[cycle[i]][cycle[(i + 1) % length]] for i in range(length)):
                cycles.append(cycle)
    return cycles


def main():
    pairs = 0
    ties = []
    for a, b in combinations(range(1, 65), 2):
        d, s = b - a, a + b
        c, h = Fraction(s, 2), Fraction(d, 2)
        r = ratio(a, b)
        check(a + b == 2 * a + d == 2 * b - d, "sum coordinates")
        check(a * b == a*a + a*d == b*b - b*d, "product coordinates")
        check(a*b == c*c - h*h, "center product")
        check(r == Fraction(a**d, 1) / Fraction(b, a)**a, "gap ratio")
        check(b**4 - a**4 == 8*c**3*h + 8*c*h**3, "odd quartic part")
        check(b**4 + a**4 == 2*c**4 + 12*c*c*h*h + 2*h**4, "even quartic part")
        check((a**4)**b == a**(4*b) == (a**b)**4, "power transport A")
        check((b**4)**a == b**(4*a) == (b**a)**4, "power transport B")
        check(Fraction((4*a)**b, (4*b)**a) == 4**d * r, "scaled bases")
        check(ratio(4*a, 4*b) == (4**d * r)**4, "scaled arguments")
        expected = -1 if a == 1 or (a, b) == (2, 3) else (0 if (a, b) == (2, 4) else 1)
        check(sign(r - 1) == expected, "raw comparison classification")
        if r == 1:
            ties.append((a, b))
        # Squared arguments keep finite exact arithmetic moderate at this bound.
        expected_power = -1 if a == 1 else 1
        check(sign(a**(b*b) - b**(a*a)) == expected_power, "square-level stabilization")
        pairs += 1
    print(f"All {pairs} pairs 1<=A<B<=64: sum/product/gap identities, quartic parity, scale transports PASS")
    print("Raw equality pairs:", ties)
    print("Square-level comparison: smaller wins except A=1; all pairs checked")

    for a, b in combinations(range(1, 13), 2):
        check(sign(a**(b**4) - b**(a**4)) == (-1 if a == 1 else 1), "quartic level")
        # GCD-normalized fourth powers: compare actual integer powers.
        expected = not (b % a == 0)
        check(primitive_arrow(a**4, b**4) == expected, "primitive quartic classifier")
    print("All 66 pairs 1<=A<B<=12: exact fourth-power comparisons and gcd-normalized classifier PASS")

    for a in range(1, 65):
        check(ratio(a, 2*a) == Fraction(a, 2)**a, "dyadic tie formula")
    check(ratio(2, 4) == 1 and ratio(8, 16) == 65536, "dilation destroys tie")
    check(ratio(2, 3) < 1 and ratio(8, 12) > 1, "dilation reverses sign")
    check(all(primitive_arrow(u, v) for u, v in ((16, 81), (81, 256), (256, 16))), "quartic cycle")
    print("Dyadic tie: R(A,2A)=(A/2)^A; unique integer tie A=2")
    print("Scale hostile: R(2,4)=1, R(8,16)=65536; R(2,3)<1 but R(8,12)>1")
    print("GCD-normalized quartic cycle: 16 -> 81 -> 256 -> 16")

    for a in range(1, 65):
        for b in range(1, 65):
            left = a*a - 2*a*b + 2*b*b
            right = a*a + 2*a*b + 2*b*b
            check(a**4 + 4*b**4 == left*right, "Sophie Germain identity")
            check((left == 1) == ((a, b) == (1, 1)), "Sophie exceptional factor")
    print("Sophie Germain factorization: all 4096 positive pairs <=64; only (1,1) has unit factor")

    # Formal series coefficients after removing the linear log(c) term.
    # (1+z)log(1-z) - (1-z)log(1+z), independently multiplied.
    degree = 21
    log_minus = [Fraction(0)] + [-Fraction(1, j) for j in range(1, degree + 1)]
    log_plus = [Fraction(0)] + [Fraction((-1)**(j+1), j) for j in range(1, degree + 1)]
    coeff = [Fraction(0)] * (degree + 1)
    for j in range(1, degree + 1):
        coeff[j] = log_minus[j] - log_plus[j] + log_minus[j-1] + log_plus[j-1]
        expected = Fraction(0) if j % 2 == 0 else (
            Fraction(-2) if j == 1 else -Fraction(2 * j - 1, ((j-1)//2) * j))
        check(coeff[j] == expected, "odd log expansion coefficient")
    check(all(coeff[j] < 0 for j in range(3, degree + 1, 2)), "higher odd terms survive")
    print("Log-gap formal series independently multiplied through degree21; even terms zero, all higher odd terms negative")
    print("First nonlinear coefficients:", [coeff[j] for j in (3, 5, 7, 9)])

    # A positive odd-cycle control: the regular five-vertex tournament.
    t = [[(j-i) % 5 in (1, 2) for j in range(5)] for i in range(5)]
    c3, c5 = cycle_counts(t, 3), cycle_counts(t, 5)
    paths = [p for p in permutations(range(5)) if all(t[p[i]][p[i+1]] for i in range(4))]
    check((len(c3), len(c5), len(paths)) == (5, 2, 15), "regular H5 census")
    check(len(paths) == 1 + 2*len(c3) + 2*len(c5), "OCF exact control")
    check(1 + 2*len(c3) != len(paths), "discarding odd cycles above3 fails")
    print("Regular H5: 5 directed triangles, 2 directed5cycles, 15 Hamiltonian paths")
    print("Triangle-only OCF gives11, not15; quintic layer contributes4")
    print("ALL CHECKS PASSED; universal statements use companion proofs, not finite extrapolation")


if __name__ == "__main__":
    main()
