"""Exact controls for the level-22 oldform bridge; no search for Collatz coverage.

Run with Python or Python -O.  No third-party modules or mutation on import.
The modular-form identifications are cited in the companion note; finite
q-series agreement is a control, not the proof of modularity.
"""

from fractions import Fraction
from math import gcd


U2 = ((-2, 1), (-2, 0))
GOLDEN = ((0, 1), (1, 1))
SWAP = ((0, 1), (1, 0))
IDENTITY = ((1, 0), (0, 1))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def natural(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name)


def matmul(left, right, modulus=None):
    answer = tuple(tuple(sum(left[i][k] * right[k][j] for k in range(2))
                         for j in range(2)) for i in range(2))
    if modulus is not None:
        answer = tuple(tuple(x % modulus for x in row) for row in answer)
    return answer


def matpow(matrix, exponent, modulus=None):
    natural(exponent, "nonnegative integer matrix exponent")
    answer = IDENTITY
    while exponent:
        if exponent & 1:
            answer = matmul(answer, matrix, modulus)
        matrix = matmul(matrix, matrix, modulus)
        exponent >>= 1
    return answer


def eta_coefficients(bound):
    """Coefficients through q^bound of q prod(1-q^d)^(2+2[11|d])."""
    natural(bound, "positive integer series bound", 1)
    coefficients = [1] + [0] * (bound - 1)
    for degree in range(1, bound):
        for _ in range(2 + 2 * (degree % 11 == 0)):
            for j in range(bound - 1, degree - 1, -1):
                coefficients[j] -= coefficients[j - degree]
    return [0] + coefficients


def eta_by_log_derivative(bound):
    """Independent recurrence for the Euler product's logarithmic derivative."""
    natural(bound, "positive integer series bound", 1)
    c = [0] * bound
    for d in range(1, bound):
        for n in range(d, bound, d):
            c[n] += d * (2 + 2 * (d % 11 == 0))
    a = [1]
    for n in range(1, bound):
        numerator = -sum(c[j] * a[n - j] for j in range(1, n + 1))
        require(numerator % n == 0, "integral Euler recurrence")
        a.append(numerator // n)
    return [0] + a


def phi(n):
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def genus_data(level):
    """The standard Gamma0 genus formula, checked only at squarefree 11,22."""
    if type(level) is not int or level not in (11, 22):
        raise ValueError("this control is scoped to levels 11 and 22")
    primitive = sum(gcd(gcd(a, b), level) == 1
                    for a in range(level) for b in range(level))
    index = primitive // phi(level)
    divisors = [d for d in range(1, level + 1) if level % d == 0]
    cusps = sum(phi(gcd(d, level // d)) for d in divisors)
    elliptic2 = sum((x*x + 1) % level == 0 for x in range(level))
    elliptic3 = sum((x*x + x + 1) % level == 0 for x in range(level))
    genus = 1 + Fraction(index, 12) - Fraction(elliptic2, 4)
    genus -= Fraction(elliptic3, 3) + Fraction(cusps, 2)
    return index, cusps, elliptic2, elliptic3, genus


def point_count(p):
    """Literal enumeration of E: y^2+y=x^3-x^2-10x-20 over F_p."""
    return 1 + sum((y*y + y - x*x*x + x*x + 10*x + 20) % p == 0
                   for x in range(p) for y in range(p))


def word_data(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("tuple of positive integer valuations")
    p, q, b = 1, 1, 0
    for a in word:
        p, q, b = 3*p, q * 2**a, 3*b + q
    return p, q, b


def exact_cylinder(word):
    p, q, b = word_data(word)
    return ((q - b) * pow(p, -1, 2*q)) % (2*q), 2*q


def replay(source, word):
    if type(source) is not int or source <= 0 or source % 2 == 0:
        raise ValueError("positive odd integer source")
    word_data(word)
    x = source
    for a in word:
        if x == 1:
            raise ValueError("no padded ROOT edge")
        z = 3*x + 1
        actual = (z & -z).bit_length() - 1
        if actual != a:
            raise ValueError("valuation guard failed")
        x = z >> a
    return x


def main():
    checks = 0

    def check(condition, message):
        nonlocal checks
        require(condition, message)
        checks += 1

    bound = 512
    a = eta_coefficients(bound)
    check(a == eta_by_log_derivative(bound), "two independent eta expansions")
    check(a[1:12] == [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1],
          "inherited corrected coefficients")
    f2 = [a[j // 2] if j % 2 == 0 else 0 for j in range(513)]
    for n in range(1, 257):
        f2n = a[n // 2] if n % 2 == 0 else 0
        check(a[2*n] == -2*a[n] - 2*f2n, "U2 f=-2f-2f2")
        check(f2[2*n] == a[n], "U2 f2=f")
    for p in (2, 3, 5, 7, 13, 17, 19, 23, 29, 31):
        check(a[p] == p + 1 - point_count(p), "good-prime Frobenius trace")
    check(genus_data(11) == (12, 2, 0, 0, 1), "level11 genus")
    check(genus_data(22) == (36, 4, 0, 0, 2), "level22 genus")
    check(matpow(U2, 4) == ((-4, 0), (0, -4)), "four-step negative scalar")
    check(matpow(U2, 8) == ((16, 0), (0, 16)), "eight-step scalar")
    check(matmul(matmul(SWAP, U2, 3), SWAP, 3) == GOLDEN,
          "golden module after explicit basis change mod3")
    for exponent in range(1, 9):
        check((matpow(U2, exponent, 3) == IDENTITY) == (exponent == 8),
              "exact order eight")
    orbits = set()
    for point in ((a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)):
        orbit, x = [], point
        for _ in range(8):
            orbit.append(x)
            x = tuple(sum(U2[i][j] * x[j] for j in range(2)) % 3 for i in range(2))
        check(x == point and len(set(orbit)) == 8, "every nonzero phase has period8")
        orbits.add(frozenset(orbit))
    check(len(orbits) == 1, "one nonzero phase cycle")
    check((-2) % 9 != 1 % 9, "trace blocks mod9 golden conjugacy")
    check(2 % 9 != (-1) % 9, "determinant blocks mod9 golden conjugacy")
    columns = ((0, 1), (1, 0), (-1, -1), (1, 2),
               (0, -1), (-1, 0), (1, 1), (-1, -2))
    for j, column in enumerate(columns):
        image = tuple(sum(U2[i][k] * column[k] for k in range(2)) for i in range(2))
        weight = 1 if j % 2 == 0 else 2
        check(image == tuple(weight*x for x in columns[(j+1) % 8]),
              "positive eight-state lift intertwines AJ=JP")
        check(columns[(j+4) % 8] == tuple(-x for x in column),
              "positive cone has nonpointed quotient")
        state, scale = j, 1
        for _ in range(8):
            scale *= 1 if state % 2 == 0 else 2
            state = (state + 1) % 8
        check(state == j and scale == 16, "positive lift eighth power")
    check(abs(columns[0][0]*columns[1][1] - columns[0][1]*columns[1][0]) == 1,
          "surjective integral quotient")
    for v in ((x, y) for x in range(-5, 6) for y in range(-5, 6)):
        fourth = tuple(sum(matpow(U2, 4)[i][j] * v[j] for j in range(2)) for i in range(2))
        check(fourth == tuple(-4*x for x in v), "cone obstruction identity")
    check(word_data((1, 2)) == (9, 8, 5), "ordered carry12")
    check(word_data((2, 1)) == (9, 8, 7), "ordered carry21")
    check(exact_cylinder((1, 2)) == (11, 16), "legal cylinder12")
    check(exact_cylinder((2, 1)) == (9, 16), "legal cylinder21")
    for t in range(64):
        check(replay(11 + 16*t, (1, 2)) == 13 + 18*t, "literal12 control")
        check(replay(9 + 16*t, (2, 1)) == 11 + 18*t, "literal21 control")
    check(replay(3, (1,)) == replay(13, (3,)) == 5,
          "equal target, opposite original-source payment")
    for bad in (True, 1.0, 0, 2):
        try:
            replay(bad, ())
        except ValueError:
            check(True, "malformed source rejected")
        else:
            raise AssertionError("malformed source accepted")
    print("level22 oldform bridge: FINITE-EXACT controls; modular identification CITED")
    print("eta coefficients q1..q11:", a[1:12])
    print("two eta expansions agree through q512; U2 relations through q256")
    print("Gamma0(11): index12 cusps2 elliptic0/0 genus1")
    print("Gamma0(22): index36 cusps4 elliptic0/0 genus2; old dimension2, new0 (cited)")
    print("U2 matrix:", U2, "; characteristic polynomial X^2+2X+2")
    print("U2^4=-4I; U2^8=16I; no nonzero pointed invariant real cone")
    print("mod3: golden after SWAP; order8; one nonzero eight-point orbit")
    print("mod9: trace and determinant obstruct golden conjugacy")
    print("positive eight-state cycle P: P^8=16I, AJ=JP; J is integrally surjective")
    print("J maps the positive cone onto the whole real plane; positivity is lost at quotient")
    print("ordered Collatz carriers12=(9,8,5),21=(9,8,7); guarded APs distinct")
    print("same target5: source3 grows, source13 shrinks")
    print("No claim of new positive Collatz weights or universal ROOT coverage")
    print("checks:", checks)


if __name__ == "__main__":
    main()
