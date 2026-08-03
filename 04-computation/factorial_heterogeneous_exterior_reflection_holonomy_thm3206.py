#!/usr/bin/env python3
"""Exact finite-field controls for THM-3206.

The proof in the theorem is symbolic.  This companion independently exhausts
the off-discriminant parameter universes over F_p for p=3,5,7,11,13 and
checks longer heterogeneous words over the three smallest fields.  It uses
only integer arithmetic modulo p.
"""

from itertools import product


PRIMES = (3, 5, 7, 11, 13)


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def mmul(a, b, p):
    rows = len(a)
    mid = len(b)
    cols = len(b[0])
    require(len(a[0]) == mid, "matrix dimensions")
    return [
        [sum(a[i][k] * b[k][j] for k in range(mid)) % p
         for j in range(cols)]
        for i in range(rows)
    ]


def mscale(c, a, p):
    return [[c * x % p for x in row] for row in a]


def msub(a, b, p):
    return [
        [(a[i][j] - b[i][j]) % p for j in range(len(a[0]))]
        for i in range(len(a))
    ]


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def det2(a, p):
    return (a[0][0] * a[1][1] - a[0][1] * a[1][0]) % p


def trace2(a, p):
    return (a[0][0] + a[1][1]) % p


PI = [[0, 1, 0], [0, 0, 1]]
IOTA = [[1, -1], [1, 0], [0, 1]]
E0 = [[1], [0], [0]]
LEFT = [[1, -1, 1]]


def delta(s, v, p):
    return (1 - 4 * s * v) % p


def fmat(s, v, p):
    return [
        [-(1 + 2 * v) % p, 2 * v % p],
        [-2 * (s + v + 1) % p, (2 * v + 1) % p],
    ]


def emat(s, v, p):
    return [
        [0, (2 * s + 1) % p, -1 % p],
        [0, -(1 + 2 * v) % p, 2 * v % p],
        [0, -2 * (s + v + 1) % p, (2 * v + 1) % p],
    ]


def legendre(delta_value, p):
    chi = pow(delta_value, (p - 1) // 2, p)
    require(chi in (1, p - 1), "Euler criterion")
    return chi


def dmat(s, v, p):
    scalar = (-legendre(delta(s, v, p), p) * s) % p
    return mscale(scalar, emat(s, v, p), p), scalar


def states(p):
    return [
        (s, v)
        for s in range(1, p)
        for v in range(p)
        if delta(s, v, p) != 0
    ]


def is_scalar_2(a, p):
    return a[0][1] % p == 0 and a[1][0] % p == 0 \
        and (a[0][0] - a[1][1]) % p == 0


def check_one_blocks():
    checks = 0
    require(mmul(PI, IOTA, 101) == eye(2), "split carrier")
    for p in PRIMES:
        for s, v in states(p):
            de = delta(s, v, p)
            f = fmat(s, v, p)
            e = emat(s, v, p)
            require(e == mmul(mmul(IOTA, f, p), PI, p),
                    "iota F pi factorization")
            require(trace2(f, p) == 0, "reflection trace")
            require(det2(f, p) == (-de) % p, "reflection determinant")
            require(mmul(f, f, p) == mscale(de, eye(2), p),
                    "reflection square")
            require(mmul(e, E0, p) == [[0], [0], [0]],
                    "common right kernel")
            require(mmul(LEFT, e, p) == [[0, 0, 0]],
                    "common left conormal")
            checks += 1
    return checks


def check_pairs():
    checks = 0
    scalar_checks = 0
    for p in PRIMES:
        ss = states(p)
        for (s1, v1), (s2, v2) in product(ss, repeat=2):
            f1 = fmat(s1, v1, p)
            f2 = fmat(s2, v2, p)
            q = mmul(f2, f1, p)
            tau = (-2 * (2 * s1 * v2 + 2 * s2 * v1 - 1)) % p
            de1 = delta(s1, v1, p)
            de2 = delta(s2, v2, p)
            determinant = de1 * de2 % p
            require(trace2(q, p) == tau, "pair trace")
            require(det2(q, p) == determinant, "pair determinant")
            q2 = mmul(q, q, p)
            require(
                msub(msub(q2, mscale(tau, q, p), p),
                     mscale(-determinant, eye(2), p), p)
                == [[0, 0], [0, 0]],
                "pair Cayley-Hamilton",
            )

            comm = msub(mmul(f2, f1, p), mmul(f1, f2, p), p)
            discriminant = (tau * tau - 4 * determinant) % p
            require(det2(comm, p) == (-discriminant) % p,
                    "commutator discriminant")

            scalar = is_scalar_2(q, p)
            require(scalar == ((s1, v1) == (s2, v2)),
                    "sharp scalar pair classification")
            checks += 1
            scalar_checks += 1
    return checks, scalar_checks


def check_words():
    checks = 0
    regimes = ((3, 5), (5, 4), (7, 3))
    for p, max_length in regimes:
        ss = states(p)
        for length in range(1, max_length + 1):
            for word in product(ss, repeat=length):
                f_product = eye(2)
                e_product = eye(3)
                d_product = eye(3)
                scalar_product = 1
                determinant_product = 1
                for s, v in word:
                    f = fmat(s, v, p)
                    e = emat(s, v, p)
                    d, scalar = dmat(s, v, p)
                    f_product = mmul(f, f_product, p)
                    e_product = mmul(e, e_product, p)
                    d_product = mmul(d, d_product, p)
                    scalar_product = scalar * scalar_product % p
                    determinant_product = (
                        -determinant_product * delta(s, v, p)
                    ) % p

                carrier = mmul(mmul(IOTA, f_product, p), PI, p)
                require(e_product == carrier, "word carrier factorization")
                require(d_product == mscale(scalar_product, carrier, p),
                        "word physical scalar factor")
                require(det2(f_product, p) == determinant_product,
                        "word determinant")
                require(determinant_product != 0, "word invertibility")
                require(mmul(e_product, E0, p) == [[0], [0], [0]],
                        "word common right kernel")
                require(mmul(LEFT, e_product, p) == [[0, 0, 0]],
                        "word common left conormal")
                # The last two rows/columns of the carrier are f_product,
                # so its image is exactly the common plane and has rank two.
                require([row[1:] for row in e_product[1:]] == f_product,
                        "word exact image plane")
                checks += 1
    return checks


def hostile_control():
    p = 5
    f1 = fmat(1, 1, p)
    f2 = fmat(3, 1, p)
    q = mmul(f2, f1, p)
    q2 = mmul(q, q, p)
    require(delta(1, 1, p) == 2, "hostile first discriminant")
    require(delta(3, 1, p) == 4, "hostile second discriminant")
    require(f1 == [[2, 2], [4, 3]], "hostile F1")
    require(f2 == [[2, 2], [0, 3]], "hostile F2")
    require(q == [[2, 0], [2, 4]], "hostile pair")
    require(trace2(q, p) == 1, "hostile trace")
    require(det2(q, p) == 3, "hostile determinant")
    require(q2 == [[4, 0], [2, 1]], "hostile square")
    require(not is_scalar_2(q, p), "hostile pair nonscalar")
    require(not is_scalar_2(q2, p), "hostile square nonscalar")
    return f1, f2, q, q2


def compact(a):
    return ";".join(",".join(str(x) for x in row) for row in a)


def main():
    one = check_one_blocks()
    pairs, scalar = check_pairs()
    words = check_words()
    f1, f2, q, q2 = hostile_control()
    print("HETEROGENEOUS FACTORIAL EXTERIOR REFLECTION EXACT CONTROL")
    print(f"one_block_parameter_checks={one}")
    print(f"ordered_pair_invariant_checks={pairs}")
    print(f"sharp_scalar_pair_classification_checks={scalar}")
    print(f"heterogeneous_word_checks={words}")
    print(f"p5_F1={compact(f1)}")
    print(f"p5_F2={compact(f2)}")
    print(f"p5_F2F1={compact(q)}")
    print(f"p5_(F2F1)^2={compact(q2)}")
    print("common_kernel=e0")
    print("common_image=y0-y1+y2=0")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
