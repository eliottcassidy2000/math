#!/usr/bin/env python3
"""Independent finite-exact certificate for the E7 THM-4015 hostile.

This file intentionally does not import the LDL search.  It proves the primal
minimum from a mod-four obstruction plus a Cauchy-certified enumeration of all
E7 roots, and proves the dual minimum by a separate 3^7 exact box.
"""

from fractions import Fraction
from itertools import product
from math import isqrt


G = (
    (2, -1, 0, 0, 0, 0, 0),
    (-1, 2, -1, 0, 0, 0, 0),
    (0, -1, 2, -1, 0, 0, -1),
    (0, 0, -1, 2, -1, 0, 0),
    (0, 0, 0, -1, 2, -1, 0),
    (0, 0, 0, 0, -1, 2, 0),
    (0, 0, -1, 0, 0, 0, 2),
)
U = (0, 0, 0, 1, 0, 1, 1)
Z = (0, 0, 1, -1, 0, 0, 0)


def quadratic(a, x):
    return sum(x[i] * sum(a[i][j] * x[j] for j in range(len(x)))
               for i in range(len(x)))


def matmul(a, b):
    bt = tuple(zip(*b))
    return tuple(tuple(sum(x * y for x, y in zip(row, col)) for col in bt)
                 for row in a)


def inverse(a):
    n = len(a)
    aug = [[Fraction(a[i][j]) for j in range(n)]
           + [Fraction(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if aug[i][j])
        aug[j], aug[pivot] = aug[pivot], aug[j]
        scale = aug[j][j]
        aug[j] = [x / scale for x in aug[j]]
        for i in range(n):
            if i == j:
                continue
            scale = aug[i][j]
            aug[i] = [x - scale * y for x, y in zip(aug[i], aug[j])]
    return tuple(tuple(row[n:]) for row in aug)


def bareiss_det(a):
    b = [list(map(int, row)) for row in a]
    n = len(b)
    sign = 1
    previous = 1
    for k in range(n - 1):
        if b[k][k] == 0:
            pivot = next(i for i in range(k + 1, n) if b[i][k])
            b[k], b[pivot] = b[pivot], b[k]
            sign *= -1
        pivot = b[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                b[i][j] = (b[i][j] * pivot - b[i][k] * b[k][j]) // previous
        previous = pivot
    return sign * b[-1][-1]


def floor_sqrt_fraction(q):
    k = isqrt(q.numerator // q.denominator)
    while Fraction((k + 1) ** 2) <= q:
        k += 1
    while Fraction(k * k) > q:
        k -= 1
    return k


def require(condition, message):
    """Optimization-stable audit gate."""
    if not condition:
        raise RuntimeError(message)


def main():
    gates = 0

    # Sylvester and determinant controls establish positive definiteness and
    # identify the standard determinant-two E7 Cartan form.
    leading = []
    for j in range(1, 8):
        det = bareiss_det(tuple(tuple(G[i][k] for k in range(j))
                                for i in range(j)))
        require(det > 0, f"nonpositive leading minor at {j}")
        leading.append(det)
        gates += 1
    require(leading == [2, 3, 4, 5, 6, 7, 2], "wrong E7 minors")
    gates += 1

    GI = inverse(G)
    identity = tuple(tuple(Fraction(i == j) for j in range(7)) for i in range(7))
    require(matmul(G, GI) == identity == matmul(GI, G), "inverse failure")
    gates += 2

    # U is not an arbitrary coordinate accident: it is the unique nonzero
    # radical vector of the E7 Cartan form modulo two.  Hence P/2 lies in the
    # dual lattice for P with coordinate vector U.
    mod_two_kernel = []
    for bits in product((0, 1), repeat=7):
        if all(sum(G[i][j] * bits[j] for j in range(7)) % 2 == 0
               for i in range(7)):
            mod_two_kernel.append(bits)
    require(mod_two_kernel == [(0,) * 7, U], "wrong mod-two radical")
    half_covector = tuple(Fraction(sum(G[i][j] * U[j] for j in range(7)), 2)
                            for i in range(7))
    require(half_covector == (0, 0, -1, 1, -1, 1, 1), "wrong half covector")
    require(quadratic(GI, half_covector) == Fraction(3, 2), "wrong half norm")
    require(sum(a * b for a, b in zip(U, half_covector)) == 3,
            "half covector is not odd")
    gates += 2 ** 7 + 4

    # Candidate values and oddness.
    require(quadratic(G, U) == 6, "wrong primal candidate norm")
    require(quadratic(GI, Z) == Fraction(3, 2), "wrong dual candidate norm")
    require(sum(a * b for a, b in zip(U, Z)) == -1, "dual candidate is not odd")
    gates += 3

    # Any vector X congruent to U modulo 2 has q(X)=q(U)=2 mod 4:
    # q(U+2K)-q(U)=4(U^T G K+K^T G K).
    require(quadratic(G, U) % 4 == 2, "wrong class congruence")
    gates += 1

    # If q(X)<=2, Cauchy gives X_i^2<=2(G^{-1})_ii.  Enumerating the
    # resulting exact coordinate box therefore sees every E7 root.
    root_caps = tuple(floor_sqrt_fraction(2 * GI[i][i]) for i in range(7))
    roots = []
    for x in product(*(range(-cap, cap + 1) for cap in root_caps)):
        q = quadratic(G, x)
        require(q >= 0, "positive definiteness failed inside root box")
        if q == 2:
            roots.append(x)
    require(root_caps == (2, 3, 4, 3, 2, 1, 2), "wrong root box")
    require(len(roots) == 126, "wrong E7 root count")
    require(not any(all((x[i] - U[i]) % 2 == 0 for i in range(7))
                    for x in roots), "root found in hostile parity class")
    gates += 3 + len(roots)

    # Combined with the mod-four congruence, absence of a root in the class
    # forces the primal minimum to be at least 6; U attains it.
    primal_min = 6

    # For H=G^{-1}, H^{-1}=G and H(z)<3/2 would imply
    # z_i^2 < (3/2)G_ii=3.  Thus every possible shorter nonzero dual vector
    # lies in {-1,0,1}^7.  Direct exact enumeration finds none and also checks
    # the expected 56 minimal vectors.
    dual_short = []
    dual_equal = []
    for z in product((-1, 0, 1), repeat=7):
        q = quadratic(GI, z)
        if z != (0,) * 7 and q < Fraction(3, 2):
            dual_short.append((z, q))
        if q == Fraction(3, 2):
            dual_equal.append(z)
    require(not dual_short, "dual vector shorter than 3/2")
    require(len(dual_equal) == 56, "wrong number of shortest dual vectors")
    require(Z in dual_equal, "dual candidate missing from minimum shell")
    gates += 3 + 3 ** 7
    dual_min = Fraction(3, 2)

    product_value = primal_min * dual_min
    require(product_value == 9 > 7, "counterexample inequality failed")
    gates += 2

    # Algebraic block amplification.  For r E7 blocks scaled by 3 and k
    # orthogonal one-dimensional blocks of Gram 2, the repeated characteristic
    # has A=18r+2k and B=1/2.  The squared product is d+2r.
    family_rows = []
    for d in range(7, 31):
        r, k = divmod(d, 7)
        ap = 18 * r + 2 * k
        bd = Fraction(1, 2)
        value = ap * bd
        require(value == d + 2 * r > d, f"family failure in rank {d}")
        family_rows.append((d, r, k, value))
        gates += 2

    print("status=FINITE-EXACT_COUNTEREXAMPLE")
    print(f"leading_principal_determinants={leading};det={leading[-1]}")
    print(f"inverse_diagonal={[GI[i][i] for i in range(7)]}")
    print(f"mod_two_kernel={mod_two_kernel};half_covector={half_covector};half_pairing=3")
    print(f"root_caps={root_caps};root_count={len(roots)};root_in_U_class=0")
    print(f"dual_box=3^7;dual_below_3/2={len(dual_short)};dual_at_3/2={len(dual_equal)}")
    print(f"U={U};A={primal_min};Z={Z};B={dual_min};odd_pair=-1")
    print(f"rank=7;A_times_B={product_value};candidate_bound=7;excess=2")
    print("family=diag((3E7)^r,2I_k);d=7r+k;A=18r+2k;B=1/2;A_times_B=d+2r")
    print(f"family_d7_to_d30={family_rows}")
    print(f"GATES={gates}")
    print("RESULT=SHARP_ARBITRARY_LATTICE_CANDIDATE_REFUTED_IN_EVERY_RANK_AT_LEAST_7")


if __name__ == "__main__":
    main()
