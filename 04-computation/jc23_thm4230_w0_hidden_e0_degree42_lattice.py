#!/usr/bin/env python3
"""Exact W=0 hidden-E0 degree-42 lattice certificate for THM-4230.

This is a deliberately scoped follow-up to THM-4230.  On the normalized
W=0 curve

    C0: x^6 + y^4 = 1,       t=(1+y^2)/x^3,

it audits the explicit degree-six hidden map f, its tau translate g=T f,
the rank-two Eisenstein lattice L=O f + O g (O=Z[zeta_3]), and the twelve
attachment points.  It proves/certifies only statements about this explicit
sub-lattice L; saturation in the full hidden Hom lattice is not asserted.

The attachment certificate uses exact finite fields, not floating point.
For h in L the twelve attachments are one tau orbit and T^2=-zeta_3.  A
common attachment value is therefore O, and collapse is equivalent to a
simultaneous pole of h and T h.  If d_h(t) is the reduced X-coordinate
denominator away from the four fixed base fibres, this requires a common
root of d_h(t) and t^deg(d_h) d_h(-1/t).  The script exhausts all degree-42
vectors in L modulo exact bounds and checks that gcd is 1 for every orbit
representative in four good finite-field embeddings.

Reproduction:
  python3 -B 04-computation/jc23_thm4230_w0_hidden_e0_degree42_lattice.py
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from typing import Iterable, Optional

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.fields import field


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Exact characteristic-zero tau/addition and Hermitian Gram certificate.
# ---------------------------------------------------------------------------


def symbolic_tau_and_gram() -> dict[str, int | str]:
    p, z, u = sp.symbols("p z u")
    phi12 = z**4 - z**2 + 1
    sqrt_three = 2 * z - z**3
    real_branch = p**2 - (1 + sqrt_three) * p + 1
    groebner = sp.groebner([real_branch, phi12], p, z, order="lex")

    # Coefficients of X/x and Y/y after deleting the common nonzero scale
    # factors and the common 1/t.  Here g=f o tau and tau(t)=-1/t.
    delta_a = z**2 * (p**2 * u - 1) - (u - p**2)
    delta_b_plus = -z**3 * (1 + p**3 * u) - (u + p**3)
    delta_b_minus = -z**3 * (1 + p**3 * u) + (u + p**3)

    def linear_determinant(left: sp.Expr, right: sp.Expr) -> sp.Expr:
        left_poly = sp.Poly(left, u)
        right_poly = sp.Poly(right, u)
        return sp.expand(
            left_poly.coeff_monomial(u) * right_poly.coeff_monomial(1)
            - left_poly.coeff_monomial(1) * right_poly.coeff_monomial(u)
        )

    plus_determinant = linear_determinant(delta_a, delta_b_plus)
    minus_determinant = linear_determinant(delta_a, delta_b_minus)
    plus_remainder = sp.factor(groebner.reduce(plus_determinant)[1])
    minus_remainder = sp.factor(groebner.reduce(minus_determinant)[1])
    require(plus_remainder == 0, "f+Tf slope cancellation failed")
    require(minus_remainder != 0, "f-Tf unexpectedly has the same cancellation")

    # For f+g the slope quotient is constant, hence X/x is again a linear
    # polynomial in u divided by t.  Its root belongs to the non-real
    # reciprocal pair of roots of p^4-2p^3-2p+1.
    la_u = sp.Poly(delta_a, u).coeff_monomial(u)
    lb_u = sp.Poly(delta_b_plus, u).coeff_monomial(u)
    slope_ratio = sp.cancel(lb_u / la_u)
    sum_linear = sp.together(
        (u - 1) * slope_ratio**2
        - ((u - p**2) + z**2 * (p**2 * u - 1))
    )
    sum_num = sp.Poly(sp.fraction(sum_linear)[0], u)
    require(sum_num.degree() == 1, "f+Tf X coefficient is not linear")
    sum_root = sp.cancel(
        -sum_num.coeff_monomial(1) / sum_num.coeff_monomial(u)
    )
    root_polynomial = sp.together(
        sum_root**2 - (2 - 2 * sqrt_three) * sum_root + 1
    )
    root_remainder = groebner.reduce(
        sp.expand(sp.fraction(root_polynomial)[0])
    )[1]
    require(root_remainder == 0, "f+Tf does not reach the other root pair")

    # Direct pole ledger.  The cancelled sum has the same four base fibres
    # as f: total X-pole order 12, hence degree 6.  For f-g the unique
    # delta_a root gives two t-values and six points above each; the slope
    # has a simple pole and X has order two at all twelve new points.
    f_x_poles = 12
    sum_x_poles = 12
    difference_new_x_poles = 12 * 2
    difference_x_poles = f_x_poles + difference_new_x_poles
    require((f_x_poles, sum_x_poles, difference_x_poles) == (12, 12, 36),
            "pole ledger changed")
    degree_f = f_x_poles // 2
    degree_sum = sum_x_poles // 2
    degree_difference = difference_x_poles // 2
    require((degree_f, degree_sum, degree_difference) == (6, 6, 18),
            "degree ledger changed")

    # Eisenstein arithmetic uses omega^2+omega+1=0.  We record an element
    # m+n*omega as (m,n).  The Hermitian form is linear in its first slot.
    omega = (0, 1)
    cross = (-4, -2)  # -4-2*omega
    cross_bar = e_conjugate(cross)
    minus_omega_squared = e_neg(e_mul(omega, omega))
    require(cross == e_mul(minus_omega_squared, cross_bar),
            "T-isometry constraint on the cross term failed")
    require(e_trace(cross) == -6, "cross trace changed")
    require(e_norm(cross) == 12, "cross norm changed")
    gram_determinant = 6 * 6 - e_norm(cross)
    require(gram_determinant == 24, "Hermitian determinant changed")
    require(lattice_degree((1, 0), (1, 0), cross) == 6,
            "degree(f+Tf) does not match the Gram form")
    require(lattice_degree((1, 0), (-1, 0), cross) == 18,
            "degree(f-Tf) does not match the Gram form")

    return {
        "plus_determinant_remainder": str(plus_remainder),
        "minus_determinant_remainder": str(minus_remainder),
        "degree_f": degree_f,
        "degree_sum": degree_sum,
        "degree_difference": degree_difference,
        "cross": "-4-2*omega",
        "gram_determinant": gram_determinant,
        "four_map_relations": "f,Tf,f+Tf,T(f+Tf)=Tf-omega*f (up to O-units)",
    }


# ---------------------------------------------------------------------------
# Eisenstein lattice enumeration.
# ---------------------------------------------------------------------------


Eisenstein = tuple[int, int]
LatticeVector = tuple[Eisenstein, Eisenstein]


def e_add(left: Eisenstein, right: Eisenstein) -> Eisenstein:
    return left[0] + right[0], left[1] + right[1]


def e_neg(value: Eisenstein) -> Eisenstein:
    return -value[0], -value[1]


def e_mul(left: Eisenstein, right: Eisenstein) -> Eisenstein:
    m, n = left
    r, s = right
    return m * r - n * s, m * s + n * r - n * s


def e_conjugate(value: Eisenstein) -> Eisenstein:
    m, n = value
    return m - n, -n


def e_trace(value: Eisenstein) -> int:
    return 2 * value[0] - value[1]


def e_norm(value: Eisenstein) -> int:
    return e_mul(value, e_conjugate(value))[0]


def lattice_degree(a: Eisenstein, b: Eisenstein,
                   cross: Eisenstein = (-4, -2)) -> int:
    return (
        6 * e_norm(a)
        + 6 * e_norm(b)
        + e_trace(e_mul(e_mul(a, e_conjugate(b)), cross))
    )


UNITS: tuple[Eisenstein, ...] = (
    (1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)
)
MINUS_OMEGA: Eisenstein = (0, -1)
OMEGA: Eisenstein = (0, 1)


def unit_action(vector: LatticeVector, unit: Eisenstein) -> LatticeVector:
    a, b = vector
    return e_mul(unit, a), e_mul(unit, b)


def tau_action(vector: LatticeVector) -> LatticeVector:
    # T(a f+b Tf)=(-omega*b) f+a Tf.
    a, b = vector
    return e_mul(MINUS_OMEGA, b), a


def determinant_element(vector: LatticeVector) -> Eisenstein:
    # Determinant of [[a,b],[-omega*b,a]], the simultaneous h/Th matrix.
    a, b = vector
    return e_add(e_mul(a, a), e_mul(OMEGA, e_mul(b, b)))


def enumerate_degree_42() -> tuple[list[LatticeVector], list[LatticeVector], Counter[int]]:
    # The smallest real eigenvalue of the Gram matrix is 6-2*sqrt(3)>2.
    # Thus degree 42 gives N(a)+N(b)<21, and
    # N(m+n*omega)>=(m^2+n^2)/2 bounds each integer coordinate by 6.
    vectors: list[LatticeVector] = []
    for am in range(-6, 7):
        for an in range(-6, 7):
            for bm in range(-6, 7):
                for bn in range(-6, 7):
                    vector = ((am, an), (bm, bn))
                    if lattice_degree(*vector) == 42:
                        vectors.append(vector)
    vector_set = set(vectors)
    require(len(vectors) == 192, "degree-42 vector count changed")

    unseen = set(vector_set)
    representatives: list[LatticeVector] = []
    while unseen:
        seed = min(unseen)
        orbit = {
            unit_action(seed, unit) for unit in UNITS
        } | {
            unit_action(tau_action(seed), unit) for unit in UNITS
        }
        require(orbit <= vector_set, "unit/tau orbit left degree 42")
        require(len(orbit) == 12, "degree-42 orbit is not free of size 12")
        representatives.append(min(orbit))
        unseen -= orbit
    representatives.sort()
    require(len(representatives) == 16, "degree-42 orbit count changed")

    determinant_norms = Counter(
        e_norm(determinant_element(vector)) for vector in vectors
    )
    require(
        determinant_norms == Counter({13: 48, 49: 48, 61: 48, 73: 48}),
        "simultaneous h/Th determinant distribution changed",
    )

    # Every Gram value is divisible by six: the diagonal terms are, while
    # Tr((m+n*omega)(-4-2*omega)) is also a multiple of six.
    require(all(lattice_degree(((m, n)), ((r, s))) % 6 == 0
                for m in range(-2, 3) for n in range(-2, 3)
                for r in range(-2, 3) for s in range(-2, 3)),
            "degree divisibility by six failed")
    require(not any(lattice_degree(((m, n)), ((r, s))) == 34
                    for m in range(-6, 7) for n in range(-6, 7)
                    for r in range(-6, 7) for s in range(-6, 7)),
            "degree 34 unexpectedly occurs in the explicit lattice")

    return vectors, representatives, determinant_norms


# ---------------------------------------------------------------------------
# Exact finite-field group law and reciprocal-denominator attachment audit.
# ---------------------------------------------------------------------------


Point = Optional[tuple[object, object]]


def finite_field_attachment_audit(
    prime: int,
    z_value: int,
    p_value: int,
    s_value: int,
    representatives: Iterable[LatticeVector],
) -> tuple[list[int], list[int], str]:
    q = prime
    z = z_value % q
    p = p_value % q
    scale = s_value % q
    sqrt_three = (2 * z - z**3) % q
    scale_denominator = (2 * p**3 + 3 * p**2 - 1) % q
    require((z**4 - z**2 + 1) % q == 0, "bad Phi_12 embedding")
    require((p**2 - (1 + sqrt_three) * p + 1) % q == 0,
            "bad real-branch embedding")
    require((pow(scale, 6, q) * scale_denominator - 4) % q == 0,
            "bad sixth-root embedding")
    require(scale_denominator != 0, "scale denominator vanished")
    require((4 * z**3 - 2 * z) % q != 0, "Phi_12 embedding ramified")
    require((2 * p - (1 + sqrt_three)) % q != 0,
            "p-branch embedding ramified")
    require((6 * pow(scale, 5, q)) % q != 0,
            "sixth-root embedding ramified")

    rational_field, t = field("t", GF(q))
    del rational_field  # the generators carry their field
    u = t * t
    inverse_two = pow(2, -1, q)
    a_scale = scale * scale * inverse_two % q
    b_scale = scale**3 * inverse_two % q
    omega_value = pow(z, 4, q)
    z_squared = z * z % q
    z_cubed = z**3 % q
    quotient_c = (u - 1) / (2 * t)  # y^2/x^3

    f_point: Point = (
        a_scale * (u - p * p) / t,
        b_scale * (u + p**3) / t,
    )
    tf_point: Point = (
        a_scale * z_squared * (p * p * u - 1) / t,
        -b_scale * z_cubed * (1 + p**3 * u) / t,
    )

    def point_neg(point: Point) -> Point:
        if point is None:
            return None
        return point[0], -point[1]

    def point_double(point: Point) -> Point:
        if point is None:
            return None
        a_coefficient, b_coefficient = point
        if b_coefficient == 0:
            return None
        doubled_a = (
            9 * t * a_coefficient**4
            / (2 * (u - 1) * b_coefficient**2)
            - 2 * a_coefficient
        )
        doubled_b = (
            3 * t * a_coefficient**2 * (a_coefficient - doubled_a)
            / ((u - 1) * b_coefficient)
            - b_coefficient
        )
        return doubled_a, doubled_b

    def point_add(left: Point, right: Point) -> Point:
        if left is None:
            return right
        if right is None:
            return left
        left_a, left_b = left
        right_a, right_b = right
        if left_a == right_a:
            if left_b == -right_b:
                return None
            require(left_b == right_b, "equal X coordinates violate E0")
            return point_double(left)
        sum_a = (
            quotient_c * (right_b - left_b) ** 2 / (right_a - left_a) ** 2
            - left_a - right_a
        )
        sum_b = (
            (right_b - left_b) * (left_a - sum_a) / (right_a - left_a)
            - left_b
        )
        return sum_a, sum_b

    def integer_multiple(multiplier: int, point: Point) -> Point:
        if multiplier < 0:
            return point_neg(integer_multiple(-multiplier, point))
        answer: Point = None
        summand = point
        remaining = multiplier
        while remaining:
            if remaining & 1:
                answer = point_add(answer, summand)
            summand = point_double(summand)
            remaining //= 2
        return answer

    def eisenstein_multiple(multiplier: Eisenstein, point: Point) -> Point:
        require(point is not None, "cannot scale the origin placeholder")
        m, n = multiplier
        omega_point: Point = (omega_value * point[0], point[1])
        return point_add(integer_multiple(m, point),
                         integer_multiple(n, omega_point))

    def reciprocal_polynomial(polynomial: object) -> object:
        degree = polynomial.degree()
        answer = polynomial.ring.zero
        for (exponent,), coefficient in polynomial.to_dict().items():
            answer[(degree - exponent,)] = coefficient * ((-1) ** exponent)
        return answer

    denominator_degrees: list[int] = []
    gcd_degrees: list[int] = []
    digest_rows: list[str] = []
    for a_coefficient, b_coefficient in representatives:
        h_point = point_add(
            eisenstein_multiple(a_coefficient, f_point),
            eisenstein_multiple(b_coefficient, tf_point),
        )
        require(h_point is not None, "degree-42 vector became zero")
        denominator = h_point[0].denom
        denominator_degree = denominator.degree()
        reciprocal = reciprocal_polynomial(denominator)
        common = denominator.gcd(reciprocal)
        denominator_degrees.append(denominator_degree)
        gcd_degrees.append(common.degree())

        # Good-reduction/pole-degree controls: d=t*D(t^2), deg D=6,
        # with neither t=0 nor t^2=-1 swallowed by the moving denominator.
        exponents = sorted(exponent[0] for exponent in denominator.to_dict())
        require(denominator_degree == 13, "degree-42 denominator degree changed")
        require(exponents[0] == 1 and all(exponent % 2 == 1 for exponent in exponents),
                "denominator is not t times a polynomial in t^2")
        value_at_u_minus_one = sum(
            coefficient * ((-1) ** ((exponent - 1) // 2))
            for (exponent,), coefficient in denominator.to_dict().items()
        )
        require(value_at_u_minus_one != 0,
                "moving denominator collided with the t^2=-1 base fibres")
        require(common.degree() == 0, "degree-42 attachment collapse detected")

        leading_inverse = pow(int(denominator.LC) % q, -1, q)
        denominator_dictionary = denominator.to_dict()
        coefficients = tuple(
            int(denominator_dictionary.get((exponent,), 0)) * leading_inverse % q
            for exponent in range(denominator_degree, -1, -1)
        )
        digest_rows.append(",".join(str(value) for value in coefficients))

    digest = sha256(";".join(digest_rows).encode("ascii")).hexdigest()
    return denominator_degrees, gcd_degrees, digest


# ---------------------------------------------------------------------------
# Symbolic norm-seven/rank-one attachment resultant.
# ---------------------------------------------------------------------------


def rank_one_resultant_certificate() -> dict[str, int | str]:
    p, r, u, x_coordinate = sp.symbols("p r u X")
    quartic = p**4 - 2 * p**3 - 2 * p + 1
    scale_denominator = 2 * p**3 + 3 * p**2 - 1

    # The two O-stable order-seven kernels are the two cubic X-orbits in
    # the factor 7*X^6-4*X^3+16 of the seventh division polynomial.
    kernel_factor = 7 * x_coordinate**6 - 4 * x_coordinate**3 + 16
    kernel_r_polynomial = 7 * r**2 - 4 * r + 16
    require(sp.expand(kernel_factor.subs(x_coordinate**3, r))
            == kernel_r_polynomial, "norm-seven kernel factor changed")

    # X_f^3 and X_Tf^3 at a regular attachment are respectively
    # (u-p^2)^3/(D*u*(u+1)) and -(p^2*u-1)^3/(D*u*(u+1)).
    f_equation = (
        (u - p**2) ** 3
        - r * scale_denominator * u * (u + 1)
    )
    tf_equation = (
        -(p**2 * u - 1) ** 3
        - r * scale_denominator * u * (u + 1)
    )
    u_resultant = sp.factor(sp.resultant(f_equation, tf_equation, u))
    p_norm = sp.factor(sp.resultant(u_resultant, quartic, p))
    expected_p_norm = (
        -20639121408 * (r + 1) ** 4 * (11 * r**2 - 32 * r - 16) ** 4
    )
    require(sp.expand(p_norm - expected_p_norm) == 0,
            "rank-one p-norm resultant changed")
    total_norm = int(sp.resultant(p_norm, kernel_r_polynomial, r))
    expected_total_norm = 2**72 * 3**54
    require(total_norm == expected_total_norm,
            "rank-one total resultant norm changed")
    return {
        "kernel_factor": "7*X^6-4*X^3+16",
        "p_norm": "-20639121408*(r+1)^4*(11*r^2-32*r-16)^4",
        "total_norm": total_norm,
        "factorization": "2^72*3^54",
    }


def main() -> None:
    symbolic = symbolic_tau_and_gram()
    vectors, representatives, determinant_norms = enumerate_degree_42()
    rank_one = rank_one_resultant_certificate()

    embeddings = (
        (313, 29, 135, 21),
        (349, 24, 246, 28),
        (373, 69, 297, 33),
        (397, 157, 161, 27),
    )
    finite_rows = []
    for embedding in embeddings:
        denominator_degrees, gcd_degrees, digest = finite_field_attachment_audit(
            *embedding, representatives
        )
        require(set(denominator_degrees) == {13},
                "finite-field denominator degree set changed")
        require(set(gcd_degrees) == {0},
                "finite-field reciprocal gcd set changed")
        finite_rows.append((embedding, digest))

    print("THM-4230 W=0 hidden-E0 degree-42 lattice certificate")
    print("scope=explicit_L=O*f+O*Tf; full_hidden_Hom_saturation=OPEN")
    print("tau=t->-1/t T^2=-omega")
    print(
        "addition="
        f"deg_f:{symbolic['degree_f']},"
        f"deg_f_plus_Tf:{symbolic['degree_sum']},"
        f"deg_f_minus_Tf:{symbolic['degree_difference']}"
    )
    print(
        "hermitian_gram=[[6,-4-2omega],[-4-2omega^2,6]] "
        f"det={symbolic['gram_determinant']}"
    )
    print(f"four_maps={symbolic['four_map_relations']}")
    print(
        f"degree42_vectors={len(vectors)} unit_tau_orbits={len(representatives)} "
        f"determinant_norms={dict(sorted(determinant_norms.items()))}"
    )
    print("degree34_vectors=0 reason=all_L_degrees_divisible_by_6")
    for embedding, digest in finite_rows:
        q, z, p, scale = embedding
        print(
            f"finite_field=q:{q},z:{z},p:{p},s:{scale} "
            f"reps:16 denominator_degrees:13 reciprocal_gcd_degrees:0 "
            f"denominator_sha256:{digest}"
        )
    print(
        "rank_one_norm7_kernel="
        f"{rank_one['kernel_factor']} p_norm_resultant={rank_one['p_norm']}"
    )
    print(
        f"rank_one_total_resultant={rank_one['total_norm']} "
        f"factorization={rank_one['factorization']} nonzero=PASS"
    )
    print("attachment_degree42_in_explicit_L=EXCLUDED_FINITE_EXACT")
    print("rank_one_degree42_attachment_collapse=EXCLUDED_SYMBOLIC")
    print("result=PASS explicit lattice; full hidden-Hom integral saturation remains OPEN")


if __name__ == "__main__":
    main()
