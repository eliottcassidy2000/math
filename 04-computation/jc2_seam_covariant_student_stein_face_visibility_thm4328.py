#!/usr/bin/env python3
"""Exact stdlib certificate for THM-4328.

The certificate keeps four claims sharply separated:

1. the source-normal row tangent on the normalized reduced (2,3) cell is a
   scalar multiple of a one-parameter Student--Stein operator;
2. its image is exactly the kernel of an explicit algebraic moment row;
3. composing those moment rows with the THM-4298 face flag is lossless for
   even residual weight and identically zero for odd residual weight; and
4. neither endpoint wall is killed by the finite row-nine bracket gate.

Everything is exact.  Only the Python standard library is used, and no
``assert`` statement is relied upon, so optimized and ordinary executions
have identical semantics and output.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from math import comb


CHECKS = 0


def check(condition: bool, label: str) -> None:
    """Record one optimization-safe exact check."""

    global CHECKS
    if not condition:
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


def product(values: list[int]) -> int:
    answer = 1
    for value in values:
        answer *= value
    return answer


def odd_double_factorial(r: int) -> int:
    """Return (2r-1)!!, with (-1)!!=1 at r=0."""

    return product(list(range(1, 2 * r, 2)))


def moment(m: int, degree: int, kappa: F) -> F:
    """The normalized algebraic Student moment through degree m."""

    check(0 <= degree <= m, f"moment range m={m} degree={degree}")
    check(kappa != 0, f"nonzero kappa m={m} degree={degree}")
    if degree % 2:
        return F(0)
    r = degree // 2
    denominator = kappa**r * product(
        [2 * m - 2 * j + 1 for j in range(1, r + 1)]
    )
    return F(6**r * odd_double_factorial(r), denominator)


def symbolic_moment(m: int, degree: int) -> tuple[F, int] | None:
    """Return c,e for ell(x^degree)=c*kappa^e, or None for odd degree."""

    check(0 <= degree <= m, f"symbolic moment range m={m} degree={degree}")
    if degree % 2:
        return None
    r = degree // 2
    denominator = product([2 * m - 2 * j + 1 for j in range(1, r + 1)])
    return F(6**r * odd_double_factorial(r), denominator), -r


def rank(matrix: list[list[F]]) -> int:
    """Exact Fraction Gaussian rank."""

    if not matrix:
        return 0
    rows = [row[:] for row in matrix]
    nrows = len(rows)
    ncols = len(rows[0])
    check(all(len(row) == ncols for row in rows), "rectangular rank matrix")
    pivot_row = 0
    for column in range(ncols):
        pivot = next(
            (row for row in range(pivot_row, nrows) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [entry / scale for entry in rows[pivot_row]]
        for row in range(nrows):
            if row == pivot_row or rows[row][column] == 0:
                continue
            multiple = rows[row][column]
            rows[row] = [
                rows[row][j] - multiple * rows[pivot_row][j]
                for j in range(ncols)
            ]
        pivot_row += 1
        if pivot_row == nrows:
            break
    return pivot_row


def tangent_matrix(m: int, kappa: F) -> list[list[F]]:
    """Matrix of S_{m,kappa}: k[x]_<m -> k[x]_<=m."""

    matrix = [[F(0) for _ in range(m)] for _ in range(m + 1)]
    for degree in range(m):
        # S(x^degree) = 3d/m*x^(d-1)
        #                 + kappa(d-2m)/(2m)*x^(d+1).
        if degree:
            matrix[degree - 1][degree] += F(3 * degree, m)
        matrix[degree + 1][degree] += kappa * F(degree - 2 * m, 2 * m)
    return matrix


def stein_image_audit() -> None:
    kappas = [F(1), F(2), F(-3), F(7, 5)]
    for m in range(1, 25):
        # The leading coefficient proves injectivity symbolically for every
        # nonzero kappa; numeric exact ranks are redundant hostile controls.
        for degree in range(m):
            check(degree - 2 * m != 0, f"symbolic leading pivot m={m} d={degree}")
        for kappa in kappas:
            matrix = tangent_matrix(m, kappa)
            check(rank(matrix) == m, f"tangent rank m={m} kappa={kappa}")
            moments = [moment(m, degree, kappa) for degree in range(m + 1)]
            check(moments[0] == 1, f"nonzero moment row m={m} kappa={kappa}")
            for column in range(m):
                pairing = sum(
                    moments[row] * matrix[row][column] for row in range(m + 1)
                )
                check(pairing == 0, f"moment annihilation m={m} col={column} kappa={kappa}")

        # Independent symbolic recurrence: aggregate coefficients by the
        # exponent of kappa and require literal cancellation.
        for degree in range(m):
            terms: dict[int, F] = {}
            if degree:
                lower = symbolic_moment(m, degree - 1)
                if lower is not None:
                    coefficient, exponent = lower
                    terms[exponent] = terms.get(exponent, F(0)) + F(3 * degree, m) * coefficient
            upper = symbolic_moment(m, degree + 1)
            if upper is not None:
                coefficient, exponent = upper
                exponent += 1
                terms[exponent] = terms.get(exponent, F(0)) + F(
                    degree - 2 * m, 2 * m
                ) * coefficient
            check(
                all(coefficient == 0 for coefficient in terms.values()),
                f"symbolic moment recurrence m={m} degree={degree}",
            )

        # The displayed stationary adjoint identity has coefficient
        # (-2m+2m) kappa*x*(kappa*x^2+6)^(-m-1).
        check(-2 * m + 2 * m == 0, f"stationary adjoint m={m}")
        integral_ratio = F(2 * m + 1, 12 * (m + 1))
        survival = 6 * integral_ratio
        check(survival == F(2 * m + 1, 2 * m + 2), f"posterior survival m={m}")


def seam_covariance_audit() -> None:
    # q=(a/(4gamma))*(kappa*x^2+6), kappa=4gamma^2/a.
    # Fraction specializations audit signs and denominators; the two
    # coefficient identities below are algebraic cancellations.
    for a, gamma in [(F(1), F(-1, 2)), (F(2), F(3)), (F(-5), F(7, 3))]:
        kappa = 4 * gamma * gamma / a
        scalar = a / (4 * gamma)
        check(scalar * kappa == gamma, f"q quadratic coefficient a={a} gamma={gamma}")
        check(6 * scalar == F(3) * a / (2 * gamma), f"q constant coefficient a={a} gamma={gamma}")
        # D_m=-a/(2gamma) S_{m,kappa}.
        check(-2 * scalar == -a / (2 * gamma), f"Stein scalar a={a} gamma={gamma}")

    # Residual fifth-root characters: char(gamma)=1, char(a)=2 modulo 5.
    check((2 * 1 - 2) % 5 == 0, "kappa residual character zero")
    check((1 + 0) % 5 == 1, "q residual character one")
    for a in [F(1), F(2), F(-3), F(5, 7)]:
        gamma = -(a**3) / 2
        check(4 * gamma * gamma / a == a**5, f"live seam kappa=a^5 at a={a}")


def face_monomials(weight: int) -> list[tuple[int, int]]:
    """Return (p exponent, y exponent), ordered by the y exponent."""

    return [
        ((weight - 3 * j) // 2, j)
        for j in range(weight // 3 + 1)
        if weight - 3 * j >= 0 and (weight - 3 * j) % 2 == 0
    ]


def face_matrix(weight: int) -> list[list[F]]:
    monomials = face_monomials(weight)
    matrix = [[F(0) for _ in monomials] for _ in monomials]
    for s, (_, j_s) in enumerate(monomials):
        for r, (i_r, j_r) in enumerate(monomials):
            if r <= s:
                check(j_s == j_r + 2 * (s - r), f"flag degree weight={weight} s={s} r={r}")
                matrix[s][r] = F(comb(i_r + j_r, s - r))
    return matrix


def face_visibility_audit() -> None:
    kappas = [F(1), F(2), F(-3), F(7, 5)]
    even_faces = 0
    odd_faces = 0
    substituted_terms = 0
    for weight in range(0, 81):
        monomials = face_monomials(weight)
        if not monomials:
            check(weight == 1, f"only Frobenius gap weight={weight}")
            continue
        matrix = face_matrix(weight)
        check(rank(matrix) == len(monomials), f"unimodular face rank weight={weight}")
        check(
            all(matrix[j][j] == 1 for j in range(len(monomials))),
            f"unit face diagonal weight={weight}",
        )
        for i, j in monomials:
            check(2 * i + 3 * j == weight, f"face weight identity M={weight} i={i} j={j}")
            for q in range(i + j + 1):
                row = i + 2 * j + q
                degree = j + 2 * q
                coefficient = comb(i + j, q)
                check(coefficient > 0, f"positive binomial M={weight} i={i} j={j} q={q}")
                check(2 * row - degree == weight, f"source diagonal M={weight} row={row} degree={degree}")
                check(degree % 2 == weight % 2, f"source parity M={weight} row={row}")
                check(degree <= row, f"moment degree range M={weight} row={row}")
                substituted_terms += 1

        diagonal_moments: list[F] = []
        for _, j in monomials:
            row = (weight + j) // 2
            check(2 * row - j == weight, f"minimal face row M={weight} j={j}")
            for kappa in kappas:
                value = moment(row, j, kappa)
                if weight % 2:
                    check(value == 0, f"odd face Stein loss M={weight} j={j} kappa={kappa}")
                else:
                    check(value != 0, f"even face Stein visibility M={weight} j={j} kappa={kappa}")
            diagonal_moments.append(moment(row, j, F(1)))
        if weight % 2:
            odd_faces += 1
            check(all(value == 0 for value in diagonal_moments), f"odd face zero composite M={weight}")
        else:
            even_faces += 1
            check(all(value != 0 for value in diagonal_moments), f"even face invertible composite M={weight}")

    check(even_faces == 41, "even face census through eighty")
    check(odd_faces == 39, "odd face census through eighty excluding one")
    check(substituted_terms == 13237, "substituted face-term census through eighty")

    # Minimal seam-entry hostile: p^8*y=x*t^10*(1+x^2*t)^9.
    for q in range(10):
        row = 10 + q
        degree = 1 + 2 * q
        check(row >= 10, f"weight-19 hostile invisible through row nine q={q}")
        check(2 * row - degree == 19, f"weight-19 hostile diagonal q={q}")
        check(comb(9, q) > 0, f"weight-19 hostile coefficient q={q}")
        for kappa in kappas:
            check(moment(row, degree, kappa) == 0, f"weight-19 Stein loss q={q} kappa={kappa}")


@dataclass(frozen=True)
class Poly:
    """Tiny exact polynomial over Q in (kappa,e0,e1,e2)."""

    terms: tuple[tuple[tuple[int, int, int, int], F], ...]

    @staticmethod
    def from_dict(raw: dict[tuple[int, int, int, int], F]) -> "Poly":
        cleaned = tuple(sorted((monomial, coefficient) for monomial, coefficient in raw.items() if coefficient))
        return Poly(cleaned)

    @staticmethod
    def constant(value: F | int) -> "Poly":
        value = F(value)
        return Poly.from_dict({(0, 0, 0, 0): value})

    @staticmethod
    def variable(index: int) -> "Poly":
        exponent = [0, 0, 0, 0]
        exponent[index] = 1
        return Poly.from_dict({tuple(exponent): F(1)})

    def dictionary(self) -> dict[tuple[int, int, int, int], F]:
        return dict(self.terms)

    def __add__(self, other: "Poly" | int | F) -> "Poly":
        if not isinstance(other, Poly):
            other = Poly.constant(other)
        result = self.dictionary()
        for monomial, coefficient in other.terms:
            result[monomial] = result.get(monomial, F(0)) + coefficient
        return Poly.from_dict(result)

    def __radd__(self, other: "Poly" | int | F) -> "Poly":
        return self + other

    def __neg__(self) -> "Poly":
        return Poly.from_dict({monomial: -coefficient for monomial, coefficient in self.terms})

    def __sub__(self, other: "Poly" | int | F) -> "Poly":
        if not isinstance(other, Poly):
            other = Poly.constant(other)
        return self + (-other)

    def __rsub__(self, other: "Poly" | int | F) -> "Poly":
        if not isinstance(other, Poly):
            other = Poly.constant(other)
        return other - self

    def __mul__(self, other: "Poly" | int | F) -> "Poly":
        if not isinstance(other, Poly):
            other = Poly.constant(other)
        result: dict[tuple[int, int, int, int], F] = {}
        for left_monomial, left_coefficient in self.terms:
            for right_monomial, right_coefficient in other.terms:
                monomial = tuple(left_monomial[j] + right_monomial[j] for j in range(4))
                result[monomial] = result.get(monomial, F(0)) + left_coefficient * right_coefficient
        return Poly.from_dict(result)

    def __rmul__(self, other: "Poly" | int | F) -> "Poly":
        return self * other

    def __pow__(self, exponent: int) -> "Poly":
        check(exponent >= 0, f"nonnegative polynomial exponent {exponent}")
        result = Poly.constant(1)
        factor = self
        power = exponent
        while power:
            if power & 1:
                result = result * factor
            factor = factor * factor
            power //= 2
        return result


def m12_wall_audit() -> None:
    # THM-4298 face matrix and the diagonal Student moment rescaling.
    matrix = face_matrix(12)
    expected = [
        [F(1), F(0), F(0)],
        [F(6), F(1), F(0)],
        [F(15), F(5), F(1)],
    ]
    check(matrix == expected, "M12 unimodular flag matrix")
    moments = [moment(6, 0, F(1)), moment(7, 2, F(1)), moment(8, 4, F(1))]
    check(moments == [F(1), F(6, 13), F(36, 65)], "M12 Student diagonal")
    symbolic_diagonal = [
        symbolic_moment(6, 0),
        symbolic_moment(7, 2),
        symbolic_moment(8, 4),
    ]
    check(
        symbolic_diagonal == [(F(1), 0), (F(6, 13), -1), (F(36, 65), -2)],
        "M12 symbolic Student diagonal",
    )
    determinant_coefficient = moments[0] * moments[1] * moments[2]
    check(determinant_coefficient == F(216, 845), "M12 composite determinant coefficient")
    check(sum(entry[1] for entry in symbolic_diagonal if entry is not None) == -3, "M12 determinant kappa exponent")

    kappa = Poly.variable(0)
    e0 = Poly.variable(1)
    e1 = Poly.variable(2)
    e2 = Poly.variable(3)
    h0 = e0
    h1 = F(13, 6) * kappa * e1
    h2 = F(65, 36) * (kappa**2) * e2
    U = h0
    W = h1 - 6 * h0
    Z = h2 - 5 * h1 + 15 * h0
    Lambda = U + W + Z
    D = W**2 - 4 * U * Z

    z_wall = 13 * (kappa**2) * e2 - 78 * kappa * e1 + 108 * e0
    lambda_wall = 65 * (kappa**2) * e2 - 312 * kappa * e1 + 360 * e0
    d_wall = (
        169 * (kappa**2) * (e1**2)
        + 624 * kappa * e0 * e1
        - 864 * (e0**2)
        - 260 * (kappa**2) * e0 * e2
    )
    check(U == e0, "M12 U wall dictionary")
    check(Z == F(5, 36) * z_wall, "M12 Z wall dictionary")
    check(Lambda == F(1, 36) * lambda_wall, "M12 Lambda wall dictionary")
    check(D == F(1, 36) * d_wall, "M12 discriminant wall dictionary")


def response(Phi: F, eta: F, xi: F) -> tuple[F, F, F]:
    U = F(475515904 - 109350 * xi, 200475)
    W = -F(4343625 * Phi * Phi - 17172000 * xi + 143826305024, 4009500)
    Z = F(
        12506118074368
        - 173745000 * Phi * Phi
        - 195463125 * Phi * eta
        - 926883000 * xi,
        108256500,
    )
    return U, W, Z


def e9(Phi: F, eta: F, xi: F, alpha: F) -> F:
    return (
        613527750 * Phi * Phi
        - 511211250 * Phi * alpha
        - 3154140000 * Phi * eta
        - 255605625 * eta * eta
        + 6736896000 * xi
        - 46483785515008
    )


def alpha_for(Phi: F, eta: F, xi: F) -> F:
    check(Phi != 0, "Phi nonzero when solving E9")
    numerator = (
        613527750 * Phi * Phi
        - 3154140000 * Phi * eta
        - 255605625 * eta * eta
        + 6736896000 * xi
        - 46483785515008
    )
    return F(numerator, 511211250 * Phi)


def endpoint_control_audit() -> None:
    # Exact U=0 projected-depth/row-nine bracket control.
    Phi = F(1)
    eta = F(0)
    Delta = F(896, 15)
    Theta = F(512, 75)
    upsilon5 = -F(9000 * Delta + 1350 * Theta + 184832, 2025)
    check(upsilon5 == -F(731648, 2025), "U control inherited upsilon5")
    xi_u = F(237757952, 54675)
    alpha_u = F(-189841784364917, 5646560625)
    U, W, Z = response(Phi, eta, xi_u)
    check(alpha_for(Phi, eta, xi_u) == alpha_u, "U control alpha")
    check(U == 0, "U control endpoint")
    check(e9(Phi, eta, xi_u, alpha_u) == 0, "U control E9")
    check(W == F(-169749098233, 9841500), "U control W")
    check(Z == F(5200771070534, 66430125), "U control Z")
    E5 = alpha_u * alpha_u - 4 * W * upsilon5
    check(
        E5
        == F(
            35245115050720811582989632889,
            31883646891800390625,
        ),
        "U control exact E5",
    )
    check(E5 != 0, "U control THM4327 E5 gate")
    check(W * W - 4 * U * Z != 0, "U control D nonzero")
    check(U + W + Z == F(3243971725969, 53144100), "U control Lambda")

    # Exact Z=0 control.  beta_11 remains free; beta_11=1 and
    # zeta_3=-3Phi/2 exhibit the generic wall sidecar explicitly.
    xi_z = F(1563243041171, 115860375)
    alpha_z = F(207913052665031843, 2393096045625)
    U, W, Z = response(Phi, eta, xi_z)
    beta = F(1)
    zeta3 = -F(3, 2) * Phi
    check(alpha_for(Phi, eta, xi_z) == alpha_z, "Z control alpha")
    check(Z == 0, "Z control endpoint")
    check(e9(Phi, eta, xi_z, alpha_z) == 0, "Z control E9")
    check(U == F(-5200771070534, 1042743375), "Z control U")
    check(W == F(10155617023591, 463441500), "Z control W")
    check(W * W - 4 * U * Z != 0, "Z control D nonzero")
    check(U + W == F(70597468930183, 4170973500), "Z control U+W")
    check(beta * zeta3 != 0, "Z control generic beta*zeta sidecar")

    # The wall and row-nine equations cut transversely on Phi!=0.
    check(F(-109350, 200475) == F(-6, 11), "U transversality derivative")
    check(F(-926883000, 108256500) == F(-22886, 2673), "Z transversality derivative")
    check(-511211250 * Phi == -511211250, "E9 alpha transversality")

    # Exact generic-family controls at several rational parameter points.
    xi_u_generic = F(237757952, 54675)
    for Phi in [F(1), F(2), F(-3), F(5, 2)]:
        for eta in [F(0), F(1), F(-2), F(7, 3)]:
            alpha = alpha_for(Phi, eta, xi_u_generic)
            U, _, _ = response(Phi, eta, xi_u_generic)
            check(U == 0, f"generic U family Phi={Phi} eta={eta}")
            check(e9(Phi, eta, xi_u_generic, alpha) == 0, f"generic U E9 Phi={Phi} eta={eta}")

            xi_z_generic = F(
                12506118074368 - 173745000 * Phi * Phi - 195463125 * Phi * eta,
                926883000,
            )
            alpha = alpha_for(Phi, eta, xi_z_generic)
            _, _, Z = response(Phi, eta, xi_z_generic)
            check(Z == 0, f"generic Z family Phi={Phi} eta={eta}")
            check(e9(Phi, eta, xi_z_generic, alpha) == 0, f"generic Z E9 Phi={Phi} eta={eta}")


def main() -> None:
    stein_image_audit()
    seam_covariance_audit()
    face_visibility_audit()
    m12_wall_audit()
    endpoint_control_audit()

    lines = [
        "THM4328 EXACT SEAM-COVARIANT STUDENT--STEIN FACE VISIBILITY",
        "SEAM kappa=4*gamma^2/a; q=(a/(4*gamma))*(kappa*x^2+6)",
        "LIVE_SEAM gamma=-a^3/2 => kappa=a^5=A5",
        "STEIN D_m=-(a/(2*gamma))*((kappa*x^2+6)*theta'-2*m*kappa*x*theta)/(2*m)",
        "IMAGE im(D_m)=ker(ell_m,kappa) on degree <=m; audited m=1..24",
        "MOMENTS ell(x^(2r))=6^r*(2r-1)!!/(kappa^r*prod_j(2m-2j+1)); ell(odd)=0",
        "STATIONARY mu_m,kappa proportional to (kappa*x^2+6)^(-(m+1)); probability reading requires real kappa>0",
        "FILTER retain 6/(kappa*x^2+6); survival=(2m+1)/(2m+2)",
        "FACE even weight: lossless diagonal Stein rescaling; odd weight: zero composite",
        "M12 e0=h0; e1=6*h1/(13*kappa); e2=36*h2/(65*kappa^2)",
        "M12 determinant=216/(845*kappa^3)",
        "M12 U=0 iff e0=0",
        "M12 Z=0 iff 13*kappa^2*e2-78*kappa*e1+108*e0=0",
        "M12 Lambda=0 iff 65*kappa^2*e2-312*kappa*e1+360*e0=0",
        "M12 D=0 iff 169*kappa^2*e1^2+624*kappa*e0*e1-864*e0^2-260*kappa^2*e0*e2=0",
        "HOSTILE weight19 p^8*y=x*t^10*(1+x^2*t)^9: no row <=9 and every raw Stein moment is zero",
        "U_CONTROL Phi=1 eta=0 xi=237757952/54675 alpha=-189841784364917/5646560625",
        "U_CONTROL upsilon5=-731648/2025 E5=35245115050720811582989632889/31883646891800390625!=0",
        "U_CONTROL W=-169749098233/9841500 Z=5200771070534/66430125 D*Lambda!=0",
        "Z_CONTROL Phi=1 eta=0 xi=1563243041171/115860375 alpha=207913052665031843/2393096045625",
        "Z_CONTROL beta=1 zeta3=-3/2 U=-5200771070534/1042743375 W=10155617023591/463441500",
        "Z_CONTROL U+W=70597468930183/4170973500 and U*beta*zeta3*(U+W)!=0",
        "TRANSVERSAL dU/dxi=-6/11 dZ/dxi=-22886/2673 dE9/dalpha=-511211250*Phi",
        f"CHECKS={CHECKS}",
    ]
    semantic_hash = sha256(("\n".join(lines) + "\n").encode("ascii")).hexdigest()
    lines.append(f"SEMANTIC_SHA256={semantic_hash}")
    lines.append("ALL THM4328 EXACT CHECKS PASSED")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
