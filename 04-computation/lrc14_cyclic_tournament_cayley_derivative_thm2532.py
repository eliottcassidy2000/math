#!/usr/bin/env python3
"""Exact referee for THM-2532: the cyclic-tournament Cayley derivative.

The calculation is dependency-free.  It constructs the two aligned Radon
marginals on the centred C_13 potential slice, verifies their oriented
difference in all twelve affine charts, proves the Cayley-transform and
sawtooth-inverse identities over ``Fraction``, checks the augmentation
lattice determinant/Smith invoice, and verifies the exact even polynomial
which recovers the THM-2523 chi_7 operator from the oriented derivative.
All checks remain active under ``python -O``.
"""

from fractions import Fraction as F


P = 13
CHI7 = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def zeros(n: int, m: int) -> list[list[F]]:
    return [[F(0) for _ in range(m)] for _ in range(n)]


def identity(n: int) -> list[list[F]]:
    out = zeros(n, n)
    for i in range(n):
        out[i][i] = F(1)
    return out


def transpose(a: list[list[F]]) -> list[list[F]]:
    return [list(row) for row in zip(*a)]


def matmul(a: list[list[F]], b: list[list[F]]) -> list[list[F]]:
    bt = transpose(b)
    return [
        [sum((x * y for x, y in zip(row, col)), F(0)) for col in bt]
        for row in a
    ]


def matadd(*matrices: list[list[F]]) -> list[list[F]]:
    n = len(matrices[0])
    m = len(matrices[0][0])
    return [
        [sum((a[i][j] for a in matrices), F(0)) for j in range(m)]
        for i in range(n)
    ]


def matscale(a: list[list[F]], scalar: F) -> list[list[F]]:
    return [[scalar * x for x in row] for row in a]


def matrix_power(a: list[list[F]], exponent: int) -> list[list[F]]:
    out = identity(len(a))
    base = a
    e = exponent
    while e:
        if e & 1:
            out = matmul(out, base)
        base = matmul(base, base)
        e //= 2
    return out


def matvec(a: list[list[F]], v: list[F]) -> list[F]:
    return [sum((x * y for x, y in zip(row, v)), F(0)) for row in a]


def dot(v: list[F], w: list[F]) -> F:
    return sum((x * y for x, y in zip(v, w)), F(0))


def determinant(a: list[list[F]]) -> F:
    m = [row[:] for row in a]
    n = len(m)
    answer = F(1)
    for c in range(n):
        pivot = next((r for r in range(c, n) if m[r][c]), None)
        if pivot is None:
            return F(0)
        if pivot != c:
            m[c], m[pivot] = m[pivot], m[c]
            answer = -answer
        value = m[c][c]
        answer *= value
        for r in range(c + 1, n):
            if not m[r][c]:
                continue
            ratio = m[r][c] / value
            for j in range(c + 1, n):
                m[r][j] -= ratio * m[c][j]
    return answer


def rank(a: list[list[F]]) -> int:
    m = [row[:] for row in a]
    rows = len(m)
    cols = len(m[0]) if rows else 0
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if m[i][c]), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        value = m[r][c]
        m[r] = [x / value for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c]:
                ratio = m[i][c]
                m[i] = [x - ratio * y for x, y in zip(m[i], m[r])]
        r += 1
        if r == rows:
            break
    return r


def faddeev_charpoly(a: list[list[F]]) -> list[int]:
    """Descending coefficients of det(xI-a), exactly."""

    n = len(a)
    b = identity(n)
    coefficients = [1]
    for k in range(1, n + 1):
        ab = matmul(a, b)
        coefficient = -sum((ab[i][i] for i in range(n)), F(0)) / k
        require(coefficient.denominator == 1, "nonintegral characteristic coefficient")
        coefficients.append(coefficient.numerator)
        b = [row[:] for row in ab]
        for i in range(n):
            b[i][i] += coefficient
    require(all(x == 0 for row in b for x in row), "Cayley-Hamilton terminal")
    return coefficients


def minor(a: list[list[F]], deleted_row: int, deleted_col: int) -> list[list[F]]:
    return [
        [x for j, x in enumerate(row) if j != deleted_col]
        for i, row in enumerate(a)
        if i != deleted_row
    ]


def shift_matrix(tau: int) -> list[list[F]]:
    """(P_tau f)(v)=f(v+tau)."""

    out = zeros(P, P)
    for v in range(P):
        out[v][(v + tau) % P] = F(1)
    return out


def multiplication_chart(tau: int) -> list[list[F]]:
    """(Q_tau f)(v)=f(tau*v), so Q_tau^T C_1 Q_tau=C_tau."""

    out = zeros(P, P)
    for v in range(P):
        out[v][(tau * v) % P] = F(1)
    return out


def radon_matrix(tau: int) -> list[list[F]]:
    out = matscale(identity(P), F(7))
    for v in range(P):
        for s in range(1, 7):
            out[v][(v - 2 * tau * s) % P] += F(1)
    return out


def tournament_derivative(tau: int) -> list[list[F]]:
    return matadd(radon_matrix(tau), matscale(radon_matrix(-tau), F(-1)))


def sawtooth_inverse(tau: int) -> list[list[F]]:
    out = zeros(P, P)
    for v in range(P):
        for d in range(1, P):
            out[v][(v + tau * d) % P] += F(2 * d - P, P)
    return out


def chi7_operator(tau: int) -> list[list[F]]:
    out = zeros(P, P)
    for r in range(P):
        for s, sign in CHI7.items():
            out[r][r] += F(2 * sign)
            out[r][(r - 2 * tau * s) % P] -= F(sign)
            out[r][(r + 2 * tau * s) % P] -= F(sign)
    return out


def augmentation_basis() -> list[list[F]]:
    out = zeros(P, P - 1)
    for i in range(P - 1):
        out[i][i] = F(1)
        out[P - 1][i] = F(-1)
    return out


def flattened(a: list[list[F]]) -> list[F]:
    return [x for row in a for x in row]


def main() -> None:
    zero = zeros(P, P)
    eye = identity(P)
    ones = [[F(1) for _ in range(P)] for _ in range(P)]
    mean_projection = matscale(ones, F(1, P))
    aug_projection = matadd(eye, matscale(mean_projection, F(-1)))

    c1 = tournament_derivative(1)
    expected_first = [0] + [1 if d % 2 else -1 for d in range(1, P)]
    require([int(x) for x in c1[0]] == expected_first, "canonical tournament row")
    require(transpose(c1) == matscale(c1, F(-1)), "skew adjacency")
    require(all(sum(row, F(0)) == 0 for row in c1), "constant kernel")
    require(rank(c1) == 12, "augmentation rank")

    # The canonical tournament is switching-equivalent to the converse of
    # the transitive tournament in the order 0<1<...<12.
    transitive = zeros(P, P)
    diagonal = zeros(P, P)
    for i in range(P):
        diagonal[i][i] = F(-1 if i % 2 else 1)
        for j in range(i + 1, P):
            transitive[i][j] = F(1)
            transitive[j][i] = F(-1)
    require(
        matmul(diagonal, matmul(c1, diagonal)) == matscale(transitive, F(-1)),
        "switching equivalence",
    )

    expected_charpoly = [
        1, 0, 78, 0, 715, 0, 1716, 0, 1287, 0, 286, 0, 13, 0
    ]
    require(faddeev_charpoly(c1) == expected_charpoly, "E_13 characteristic")
    principal_cofactors = tuple(determinant(minor(c1, i, i)) for i in range(P))
    require(principal_cofactors == (F(1),) * P, "unit principal cofactors")

    basis = augmentation_basis()
    lattice_determinants = []
    inverse_rows = []
    cayley_entries = 0
    inverse_entries = 0
    radon_entries = 0
    bridge_entries = 0
    inverse_bridge_entries = 0

    l1 = chi7_operator(1)
    expected_l_first = [0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1]
    require([int(x) for x in l1[0]] == expected_l_first, "THM-2523 first row")

    for tau in range(1, P):
        p_tau = shift_matrix(tau)
        q_tau = multiplication_chart(tau)
        r_plus = radon_matrix(tau)
        r_minus = radon_matrix(-tau)
        c_tau = tournament_derivative(tau)
        b_tau = sawtooth_inverse(tau)
        l_tau = chi7_operator(tau)

        require(c_tau == matadd(r_plus, matscale(r_minus, F(-1))), "Radon difference")
        require(transpose(c_tau) == matscale(c_tau, F(-1)), "chart skewness")
        require(
            matmul(transpose(q_tau), matmul(c1, q_tau)) == c_tau,
            "multiplicative chart conjugacy",
        )
        require(
            matmul(transpose(q_tau), matmul(l1, q_tau)) == l_tau,
            "Fano chart conjugacy",
        )
        radon_entries += P * P

        cayley_left = matmul(matadd(eye, p_tau), c_tau)
        cayley_right = matadd(p_tau, matscale(eye, F(-1)))
        require(cayley_left == cayley_right, "left Cayley derivative")
        require(matmul(c_tau, matadd(eye, p_tau)) == cayley_right, "right Cayley derivative")
        require(
            matmul(p_tau, matadd(eye, matscale(c_tau, F(-1))))
            == matadd(eye, c_tau),
            "reversible Cayley transform",
        )
        require(determinant(matadd(eye, matscale(c_tau, F(-1)))) == 2**12,
                "Cayley denominator determinant")
        cayley_entries += 3 * P * P

        require(matmul(c_tau, b_tau) == aug_projection, "right sawtooth inverse")
        require(matmul(b_tau, c_tau) == aug_projection, "left sawtooth inverse")
        require(transpose(b_tau) == matscale(b_tau, F(-1)), "skew sawtooth")
        inverse_entries += 2 * P * P
        inverse_rows.append(tuple(b_tau[0]))

        # In the lattice basis e_i-e_12, output coordinates are the first
        # twelve ambient coordinates.
        image_columns = matmul(c_tau, basis)
        coordinate_matrix = [row[:] for row in image_columns[: P - 1]]
        det_coordinate = determinant(coordinate_matrix)
        require(abs(det_coordinate) == P, "augmentation lattice determinant")
        lattice_determinants.append(det_coordinate)
        for column in range(P - 1):
            y = [image_columns[r][column] for r in range(P)]
            require(sum(y, F(0)) == 0, "integer augmentation image")
            first_moment = sum((r * y[r] for r in range(P)), F(0))
            require(first_moment.denominator == 1 and first_moment.numerator % P == 0,
                    "mod-13 moment obstruction")

        require(transpose(c_tau) == matscale(c_tau, F(-1)), "skew norm identity")
        ctc = matmul(transpose(c_tau), c_tau)
        require(ctc == matscale(matmul(c_tau, c_tau), F(-1)), "C-star-C=-C^2")
        require(sum((ctc[i][i] for i in range(P)), F(0)) == 156,
                "sum of squared singular values")

        # Exact even functional calculus bridge to THM-2523.
        powers = {e: matrix_power(c_tau, e) for e in (2, 4, 6, 8, 10, 12)}
        bridge = matadd(
            matscale(powers[12], F(3)),
            matscale(powers[10], F(223)),
            matscale(powers[8], F(1294)),
            matscale(powers[6], F(-2178)),
            matscale(powers[4], F(-10545)),
            matscale(powers[2], F(-5437)),
        )
        require(bridge == matscale(l_tau, F(256)), "even Fano functional calculus")
        l_powers = {e: matrix_power(l_tau, e) for e in range(1, 7)}
        inverse_bridge = matadd(
            matscale(l_powers[1], F(765)),
            matscale(l_powers[2], F(-667)),
            matscale(l_powers[3], F(-130)),
            matscale(l_powers[4], F(102)),
            matscale(l_powers[5], F(5)),
            matscale(l_powers[6], F(-3)),
        )
        require(
            inverse_bridge == matscale(powers[2], F(200)),
            "inverse even functional calculus",
        )
        require(matmul(c_tau, l_tau) == matmul(l_tau, c_tau), "commuting circulants")
        require(chi7_operator(-tau) == l_tau, "Fano converse evenness")
        require(tournament_derivative(-tau) == matscale(c_tau, F(-1)),
                "tournament converse oddness")
        bridge_entries += P * P
        inverse_bridge_entries += P * P

    expected_inverse_first = tuple(
        [F(0)] + [F(2 * d - P, P) for d in range(1, P)]
    )
    require(inverse_rows[0] == expected_inverse_first, "canonical sawtooth row")
    require(set(abs(x) for x in lattice_determinants) == {P}, "Smith determinant invoice")

    # adj(C)=J and the regularized inverse are consequences of the unit
    # cofactors and the sawtooth identity; check them directly as well.
    regularized_inverse = matadd(sawtooth_inverse(1), matscale(ones, F(1, P * P)))
    require(matmul(matadd(c1, ones), regularized_inverse) == eye,
            "regularized inverse")
    require(determinant(matadd(c1, ones)) == P * P, "regularized determinant")

    # Functional-algebra dimensions.  C generates the full rational
    # circulant algebra.  L and C^2 generate the same converse-even part;
    # multiplication by C gives the six-dimensional skew part.
    c_powers = [matrix_power(c1, exponent) for exponent in range(P)]
    l_powers = [matrix_power(l1, exponent) for exponent in range(7)]
    c2 = matrix_power(c1, 2)
    c2_powers = [matrix_power(c2, exponent) for exponent in range(7)]
    skew_module = [matmul(c1, matrix_power(l1, exponent)) for exponent in range(6)]
    full_algebra_rank = rank([flattened(a) for a in c_powers])
    even_l_rank = rank([flattened(a) for a in l_powers])
    even_c2_rank = rank([flattened(a) for a in c2_powers])
    skew_module_rank = rank([flattened(a) for a in skew_module])
    require(full_algebra_rank == 13, "C generates full circulant algebra")
    require((even_l_rank, even_c2_rank, skew_module_rank) == (7, 7, 6),
            "converse-even/skew algebra dimensions")

    # Deterministic augmentation probes check the operator/norm formula away
    # from basis vectors too.
    norm_probes = 0
    for seed in range(48):
        vector = [F(((seed + 5) * (r + 3) + r * r) % 17 - 8) for r in range(P)]
        mean = sum(vector, F(0)) / P
        vector = [x - mean for x in vector]
        if all(x == 0 for x in vector):
            continue
        cv = matvec(c1, vector)
        require(sum(cv, F(0)) == 0, "probe remains centred")
        require(dot(cv, cv) == dot(vector, matvec(matscale(matmul(c1, c1), F(-1)), vector)),
                "probe norm identity")
        recovered = matvec(sawtooth_inverse(1), cv)
        require(recovered == vector, "probe sawtooth recovery")
        norm_probes += 1

    print("THM-2532 exact cyclic-tournament Cayley-derivative referee")
    print(
        "canonical_tournament: first_row="
        f"{tuple(expected_first)} rank=12 switching=-transitive"
    )
    print("cayley_identity: (I+P_tau)C_tau=P_tau-I and P_tau=(I+C_tau)/(I-C_tau)")
    print(
        "spectrum: char=x*(x^12+78*x^10+715*x^8+1716*x^6+"
        "1287*x^4+286*x^2+13) trace(C^T*C)=156"
    )
    print(
        "sawtooth_inverse: first_row=(1/13)*"
        "(0,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11)"
    )
    print(
        "augmentation_lattice: determinants="
        f"{tuple(int(x) for x in lattice_determinants)} Smith=(1^11,13) "
        "image=ker(sum_r r*y_r mod 13)"
    )
    print("cofactors: all 13 principal minors=1 adj(C)=J det(C+J)=169")
    print(
        "fano_bridge: 256*L=3*C^12+223*C^10+1294*C^8-2178*C^6-"
        "10545*C^4-5437*C^2"
    )
    print(
        "inverse_even_bridge: 200*C^2=765*L-667*L^2-130*L^3+"
        "102*L^4+5*L^5-3*L^6"
    )
    print(
        f"algebras: Q[C]=Circ_13 rank={full_algebra_rank}; "
        f"Q[L]=Q[C^2] rank={even_l_rank}; C*Q[L] rank={skew_module_rank}"
    )
    print(
        f"exact_checks: radon_entries={radon_entries} cayley_entries={cayley_entries} "
        f"inverse_entries={inverse_entries} bridge_entries={bridge_entries} "
        f"inverse_bridge_entries={inverse_bridge_entries} "
        f"norm_probes={norm_probes}"
    )
    print(
        "VERIFIED: the oriented aligned-Radon difference is the cyclic-tournament "
        "Cayley derivative; chi7 is an even polynomial shadow and cannot recover "
        "the missing chart sign."
    )


if __name__ == "__main__":
    main()
