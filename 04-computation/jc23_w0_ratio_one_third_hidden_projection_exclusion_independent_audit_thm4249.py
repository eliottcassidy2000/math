#!/usr/bin/env python3
"""Clean-room exact audit of THM-4249's degree-42, R=1/3 lane.

This standard-library-only program independently reconstructs all 3,168
degree-42 residual vectors, their 132 mu_6-by-C_12 orbits, and one complete
twelve-node attachment orbit.  It then checks a uniform hidden-projector
obstruction at the common ratio R=1/3 and a split good specialization to
F_397.

For a residual vector m=b*f+c*g+d*h, put

    H=P_L(m)=A*f+B*g,
    A=2*b+omega^2*d,  B=2*c+d,  q(H)=12*K,

where K is one of 13,11,10,7,5.  The E[3] ratio R=1/3 lies in every one of
the corresponding d-kernels, so it contributes one incidence for every one
of the 132 residual map orbits, not only for N(d)=3.

If m collapses the twelve attachments, H does too.  The exact source action
T*f=g, T*g=-omega*f has T^6=-1 and T^8=omega on the hidden plane, so the
common value of H is O.  The first two nodes then give

    A*F+B*G=O,       -omega*B*F+A*G=O,

where F=f(Q_0), G=g(Q_0).  Thus Delta*F=O for
Delta=A^2+omega*B^2.

The obstruction is uniform in all 3,168 vectors.  For
pi=omega^2-1, O/(pi)=F_3 and omega=1.  If pi divided Delta, then
A^2+B^2=0 in F_3, forcing pi|A and pi|B.  Dividing H by pi would give
hidden degree 4*K, impossible because the hidden Gram takes values in 6*Z
and none of the five K-values is divisible by three.  Hence 3 does not
divide N(Delta).  At the exact split specialization below F has additive
order 18.  Delta*F=O would imply
N(Delta)*F=O, hence 18|N(Delta), a contradiction.

The per-orbit direct evaluations are retained as hostile controls, not as the
load-bearing proof.  This is a FINITE-EXACT good-reduction certificate
relative to THM-4230/4241/4249.  It excludes all 132 degree-42 incidences at
the common ratio R=1/3.  It does not close the other marked ratios, W=0,
M=12, seam entry, or JC(2).
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import isqrt


E = tuple[int, int]                 # m+n*omega, omega^2+omega+1=0
Vector = tuple[E, E, E, E]          # coefficients of (u,f,g,h)
Point = tuple[int, int] | None

ZERO: E = (0, 0)
ONE: E = (1, 0)
OMEGA: E = (0, 1)
OMEGA2: E = (-1, -1)
MINUS_OMEGA: E = (0, -1)
PI: E = (-2, -1)                   # omega^2-1
UNITS: tuple[E, ...] = (
    ONE, (-1, 0), OMEGA, MINUS_OMEGA, OMEGA2, (1, 1),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def e_add(left: E, right: E) -> E:
    return left[0] + right[0], left[1] + right[1]


def e_neg(value: E) -> E:
    return -value[0], -value[1]


def e_sub(left: E, right: E) -> E:
    return left[0] - right[0], left[1] - right[1]


def e_scale(multiplier: int, value: E) -> E:
    return multiplier * value[0], multiplier * value[1]


def e_mul(left: E, right: E) -> E:
    m, n = left
    r, s = right
    return m * r - n * s, m * s + n * r - n * s


def e_conjugate(value: E) -> E:
    return value[0] - value[1], -value[1]


def e_norm(value: E) -> int:
    answer = e_mul(value, e_conjugate(value))
    require(answer[1] == 0, "Eisenstein norm left Z")
    return answer[0]


def e_trace(value: E) -> int:
    return 2 * value[0] - value[1]


def e_mod_pi(value: E) -> int:
    """O/(omega^2-1)=F_3, with omega mapped to 1."""

    return (value[0] + value[1]) % 3


def hidden_degree(left: E, right: E) -> int:
    cross = e_trace(e_mul(e_mul(left, e_conjugate(right)), (-4, -2)))
    return 6 * e_norm(left) + 6 * e_norm(right) + cross


def hidden_coefficients(vector: Vector) -> tuple[E, E]:
    _, b, c, d = vector
    return e_add(e_scale(2, b), e_mul(OMEGA2, d)), e_add(e_scale(2, c), d)


def full_degree(vector: Vector) -> int:
    a, _, _, d = vector
    left, right = hidden_coefficients(vector)
    numerator = hidden_degree(left, right)
    require(numerator % 4 == 0, "glued hidden numerator lost divisibility by four")
    return 4 * e_norm(a) + e_norm(d) + numerator // 4


def vector_scale(scalar: E, vector: Vector) -> Vector:
    return tuple(e_mul(scalar, value) for value in vector)  # type: ignore[return-value]


def tau(vector: Vector) -> Vector:
    """T*u=-omega*u, T*f=g, T*g=-omega*f, T*h=omega^2*h-omega*f."""

    a, b, c, d = vector
    return (
        e_mul(MINUS_OMEGA, a),
        e_add(e_mul(MINUS_OMEGA, c), e_mul(MINUS_OMEGA, d)),
        b,
        e_mul(OMEGA2, d),
    )


def rho(vector: Vector) -> Vector:
    """Precomposition by (x,y)->(-x,-y).

    Here rho fixes u,f,g, negates v, and 2h=v+omega^2*f+g, hence
    rho(h)=omega^2*f+g-h.
    """

    a, b, c, d = vector
    return (
        a,
        e_add(b, e_mul(OMEGA2, d)),
        e_add(c, d),
        e_neg(d),
    )


def orbit(vector: Vector) -> set[Vector]:
    answer: set[Vector] = set()
    current = vector
    for _ in range(12):
        answer.update(vector_scale(unit, current) for unit in UNITS)
        current = tau(current)
    require(current == vector, "T^12 failed")
    return answer


def enumerate_degree42_residual() -> tuple[
    set[Vector], list[Vector], Counter[int], Counter[tuple[int, int]], int
]:
    """Reconstruct the complete characteristic-zero degree-42 residual.

    The hidden Hermitian matrix has eigenvalues 6+-sqrt(12).  Since
    sqrt(12)<4, q(A,B)>2(N(A)+N(B)); q=156 gives N(A)+N(B)<78.
    The identity 4N(m+n*omega)=3m^2+(m-2n)^2 (and its symmetric version)
    then makes the coordinate box [-10,10]^4 exhaustive: an omitted
    coordinate of absolute value at least 11 would alone give N>78.
    """

    require(12 < 16 and 3 * 11 * 11 > 4 * 78,
            "characteristic-zero hidden coordinate bound changed")
    coordinate_values = range(-10, 11)
    hidden_pairs: list[tuple[E, E]] = []
    for am in coordinate_values:
        for an in coordinate_values:
            left = am, an
            for bm in coordinate_values:
                for bn in coordinate_values:
                    right = bm, bn
                    if hidden_degree(left, right) == 156:
                        hidden_pairs.append((left, right))
    require(hidden_pairs, "hidden degree-156 shell vanished")
    max_hidden_coordinate = max(
        abs(coordinate)
        for left, right in hidden_pairs
        for coordinate in left + right
    )
    require(max_hidden_coordinate < 10,
            "hidden shell touched the hostile enumeration boundary")

    # Enumerate every possible d with N(d)<=42.  The projector filters are
    # d!=0, N(d)>=3, and K not in {1,2}; THM-4230 removes the d=0 pure-hidden
    # shell.  We do not seed the five surviving norm profiles.
    d_elements = {
        (m, n)
        for m in range(-8, 9)
        for n in range(-8, 9)
        if e_norm((m, n)) <= 42
    }
    require(not any(max(abs(m), abs(n)) == 8 for m, n in d_elements),
            "d enumeration touched its hostile boundary")

    hidden_by_degree: dict[int, list[tuple[E, E]]] = {}
    for degree in range(0, 157, 12):
        hidden_by_degree[degree] = [
            pair for pair in hidden_pairs if hidden_degree(*pair) == degree
        ]
    # hidden_pairs was initially the top shell; fill the lower shells from
    # the same exhaustive box.  Keeping this second path explicit avoids
    # deriving lower rows by scaling or by importing the parent census.
    for degree in range(0, 157, 12):
        if degree == 156:
            continue
        rows: list[tuple[E, E]] = []
        for am in coordinate_values:
            for an in coordinate_values:
                left = am, an
                for bm in coordinate_values:
                    for bn in coordinate_values:
                        right = bm, bn
                        if hidden_degree(left, right) == degree:
                            rows.append((left, right))
        hidden_by_degree[degree] = rows

    vectors: set[Vector] = set()
    for d in d_elements:
        norm_d = e_norm(d)
        if d == ZERO or norm_d < 3 or (42 - norm_d) % 3:
            continue
        k_hidden = (42 - norm_d) // 3
        if k_hidden < 0 or k_hidden in (1, 2):
            continue
        hidden_q = 12 * k_hidden
        if hidden_q not in hidden_by_degree:
            continue
        omega2_d = e_mul(OMEGA2, d)
        for left, right in hidden_by_degree[hidden_q]:
            differences = (
                left[0] - omega2_d[0], left[1] - omega2_d[1],
                right[0] - d[0], right[1] - d[1],
            )
            if any(value % 2 for value in differences):
                continue
            b = differences[0] // 2, differences[1] // 2
            c = differences[2] // 2, differences[3] // 2
            vector: Vector = (ZERO, b, c, d)
            require(full_degree(vector) == 42, "residual vector lost degree 42")
            vectors.add(vector)

    profile = Counter(
        (e_norm(vector[3]), hidden_degree(*hidden_coefficients(vector)) // 12)
        for vector in vectors
    )
    require(len(vectors) == 3168, "complete degree-42 residual count changed")
    require(profile == Counter({
        (3, 13): 672,
        (9, 11): 576,
        (12, 10): 864,
        (21, 7): 768,
        (27, 5): 288,
    }), "complete degree-42 (N(d),K) profile changed")

    unseen = set(vectors)
    representatives: list[Vector] = []
    orbit_sizes: Counter[int] = Counter()
    orbit_index: dict[Vector, int] = {}
    while unseen:
        representative = min(unseen)
        current_orbit = orbit(representative)
        require(current_orbit <= vectors, "mu6-by-C12 orbit left the lane")
        index = len(representatives)
        for row in current_orbit:
            require(row not in orbit_index, "symmetry orbits overlapped")
            orbit_index[row] = index
        unseen -= current_orbit
        representatives.append(representative)
        orbit_sizes[len(current_orbit)] += 1

    require(len(representatives) == 132
            and orbit_sizes == Counter({24: 132}),
            "complete degree-42 residual orbit census changed")

    # rho exchanges the two radical/node orbits.  It commutes with T,
    # preserves the exact residual set, and in fact fixes each of the 132
    # quotient classes individually.
    for vector in vectors:
        require(rho(rho(vector)) == vector, "rho^2 failed")
        require(rho(tau(vector)) == tau(rho(vector)), "rho and T stopped commuting")
        require(rho(vector) in vectors and full_degree(rho(vector)) == 42,
                "rho left the characteristic-zero lane")
    rho_fixed_classes = sum(
        orbit_index[rho(representative)] == index
        for index, representative in enumerate(representatives)
    )
    require(rho_fixed_classes == 132,
            "rho no longer fixes all 132 quotient classes")
    return vectors, representatives, orbit_sizes, profile, max_hidden_coordinate


# ---------------------------------------------------------------------------
# A split good specialization and independent finite-group arithmetic.
# ---------------------------------------------------------------------------

Q = 397
ZETA12 = 157
HIDDEN_A = 161
HIDDEN_SCALE = 27
SOURCE_X = 15
SOURCE_Y = 28


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            return False
    return True


def inv(value: int) -> int:
    value %= Q
    require(value != 0, "attempted finite-field inversion of zero")
    return pow(value, Q - 2, Q)


def multiplicative_order(value: int) -> int:
    current = 1
    for order in range(1, Q):
        current = current * value % Q
        if current == 1:
            return order
    raise RuntimeError("nonzero finite-field element had no order")


OMEGA_F = pow(ZETA12, 4, Q)
ZETA6 = pow(ZETA12, 2, Q)
ZETA4 = pow(ZETA12, 3, Q)


def point_on_curve(point: Point) -> bool:
    return point is None or (point[1] * point[1] - point[0] ** 3 - 1) % Q == 0


def point_neg(point: Point) -> Point:
    return None if point is None else (point[0], -point[1] % Q)


def point_add(left: Point, right: Point) -> Point:
    require(point_on_curve(left) and point_on_curve(right),
            "elliptic addition received an off-curve point")
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right:
        if (y_left + y_right) % Q == 0:
            return None
        require(y_left == y_right and y_left != 0,
                "same-x elliptic branch was neither doubling nor inverse")
        slope = 3 * x_left * x_left * inv(2 * y_left) % Q
    else:
        slope = (y_right - y_left) * inv(x_right - x_left) % Q
    x_new = (slope * slope - x_left - x_right) % Q
    answer = x_new, (slope * (x_left - x_new) - y_left) % Q
    require(point_on_curve(answer), "elliptic addition left E0")
    return answer


def integer_multiple(multiplier: int, point: Point) -> Point:
    if multiplier < 0:
        return integer_multiple(-multiplier, point_neg(point))
    answer: Point = None
    summand = point
    remaining = multiplier
    while remaining:
        if remaining & 1:
            answer = point_add(answer, summand)
        summand = point_add(summand, summand)
        remaining //= 2
    return answer


def integer_multiple_naive(multiplier: int, point: Point) -> Point:
    if multiplier < 0:
        return integer_multiple_naive(-multiplier, point_neg(point))
    answer: Point = None
    for _ in range(multiplier):
        answer = point_add(answer, point)
    return answer


def omega_action(point: Point) -> Point:
    if point is None:
        return None
    answer = OMEGA_F * point[0] % Q, point[1]
    require(point_on_curve(answer), "CM omega action left E0")
    return answer


def eisenstein_multiple(coefficient: E, point: Point) -> Point:
    return point_add(
        integer_multiple(coefficient[0], point),
        integer_multiple(coefficient[1], omega_action(point)),
    )


def eisenstein_multiple_naive(coefficient: E, point: Point) -> Point:
    return point_add(
        integer_multiple_naive(coefficient[0], point),
        integer_multiple_naive(coefficient[1], omega_action(point)),
    )


def point_order(point: Point) -> int:
    # Hasse gives #E(F_Q)<=Q+1+2*sqrt(Q); the +2 makes the integer bound safe.
    bound = Q + 1 + 2 * isqrt(Q) + 2
    current: Point = None
    for order in range(1, bound + 1):
        current = point_add(current, point)
        if current is None:
            return order
    raise RuntimeError("point order exceeded the Hasse bound")


def hidden_maps(point: tuple[int, int]) -> tuple[Point, Point]:
    x_value, y_value = point
    require(x_value % Q != 0, "hidden-map x denominator vanished")
    t_value = (1 + y_value * y_value) * inv(x_value ** 3) % Q
    require(t_value != 0, "hidden-map t denominator vanished")
    half = inv(2)
    t_squared = t_value * t_value % Q
    f_point: Point = (
        HIDDEN_SCALE ** 2 * half * x_value
        * (t_squared - HIDDEN_A ** 2) * inv(t_value) % Q,
        HIDDEN_SCALE ** 3 * half * y_value
        * (t_squared + HIDDEN_A ** 3) * inv(t_value) % Q,
    )
    g_point: Point = (
        HIDDEN_SCALE ** 2 * half * x_value * ZETA6
        * (HIDDEN_A ** 2 * t_squared - 1) * inv(t_value) % Q,
        -HIDDEN_SCALE ** 3 * half * y_value * ZETA4
        * (1 + HIDDEN_A ** 3 * t_squared) * inv(t_value) % Q,
    )
    require(point_on_curve(f_point) and point_on_curve(g_point),
            "explicit hidden map missed E0")
    return f_point, g_point


def tau_point(point: tuple[int, int]) -> tuple[int, int]:
    return ZETA6 * point[0] % Q, ZETA4 * point[1] % Q


def rho_point(point: tuple[int, int]) -> tuple[int, int]:
    return -point[0] % Q, -point[1] % Q


def verify_good_specialization() -> tuple[Point, Point, list[tuple[int, int]]]:
    require(is_prime(Q) and Q not in (2, 3), "bad specialization characteristic")
    require((ZETA12 ** 4 - ZETA12 ** 2 + 1) % Q == 0
            and multiplicative_order(ZETA12) == 12,
            "zeta_12 specialization changed")
    require((OMEGA_F * OMEGA_F + OMEGA_F + 1) % Q == 0
            and OMEGA_F != 1,
            "target CM root changed")

    sqrt_three = (2 * ZETA12 - ZETA12 ** 3) % Q
    quartic = HIDDEN_A ** 4 - 2 * HIDDEN_A ** 3 - 2 * HIDDEN_A + 1
    compatible_quadratic = (
        HIDDEN_A ** 2 - (1 + sqrt_three) * HIDDEN_A + 1
    )
    scale_denominator = (2 * HIDDEN_A ** 3 + 3 * HIDDEN_A ** 2 - 1) % Q
    require(sqrt_three == 377 and sqrt_three * sqrt_three % Q == 3,
            "sqrt(3) embedding changed")
    require(quartic % Q == 0 and compatible_quadratic % Q == 0,
            "hidden quartic root is incompatible with zeta_12")
    require(scale_denominator != 0
            and (HIDDEN_SCALE ** 6 * scale_denominator - 4) % Q == 0,
            "hidden scale relation changed")

    # The triangular defining system has nonzero diagonal Jacobian.  This is
    # the explicit Hensel/good-prime audit for the simultaneous algebraic
    # constants and node radicals used by the characteristic-zero formulas.
    hensel_derivatives = (
        4 * ZETA12 ** 3 - 2 * ZETA12,
        4 * HIDDEN_A ** 3 - 6 * HIDDEN_A ** 2 - 2,
        2 * HIDDEN_A - (1 + sqrt_three),
        6 * HIDDEN_SCALE ** 5 * scale_denominator,
        24 * SOURCE_X ** 5,
        16 * SOURCE_Y ** 3,
    )
    require(all(value % Q != 0 for value in hensel_derivatives),
            "specialization tuple lost a simple root")

    require(pow(SOURCE_X, 3, Q) == inv(2)
            and pow(SOURCE_Y, 2, Q) == sqrt_three * inv(2) % Q,
            "chosen node radical branches changed")
    require((SOURCE_X ** 6 + SOURCE_Y ** 4 - 1) % Q == 0,
            "chosen node missed C0")
    require(SOURCE_X ** 6 * inv(SOURCE_Y ** 4) % Q == inv(3),
            "chosen marked ratio is not 1/3")

    nodes: list[tuple[int, int]] = []
    current = (SOURCE_X, SOURCE_Y)
    for _ in range(12):
        nodes.append(current)
        current = tau_point(current)
    require(current == nodes[0] and len(set(nodes)) == 12,
            "canonical C12 node orbit changed")

    values = [hidden_maps(node) for node in nodes]
    for index, node in enumerate(nodes):
        x_value, y_value = node
        require((x_value ** 6 + y_value ** 4 - 1) % Q == 0,
                "a canonical node left C0")
        require(x_value ** 6 * inv(y_value ** 4) % Q == inv(3),
                "marked ratio changed along the node orbit")
        f_point, g_point = values[index]
        next_f, next_g = values[(index + 1) % 12]
        require(g_point == next_f, "g=f composed with tau failed on a node")
        require(next_g == eisenstein_multiple(MINUS_OMEGA, f_point),
                "Tg=-omega*f failed on a node")
        require(hidden_maps(rho_point(node)) == (f_point, g_point),
                "rho stopped fixing f and g")

    # Exhaust every radical solution over F_397, not only exponent labels.
    x_roots = {value for value in range(Q) if 4 * pow(value, 6, Q) % Q == 1}
    y_roots = {value for value in range(Q) if 4 * pow(value, 4, Q) % Q == 3}
    all_roots = {(x_value, y_value) for x_value in x_roots for y_value in y_roots}
    canonical = set(nodes)
    rho_image = {rho_point(node) for node in nodes}
    require(len(x_roots) == 6 and len(y_roots) == 4 and len(all_roots) == 24,
            "attachment radical root count changed")
    require(canonical.isdisjoint(rho_image)
            and canonical | rho_image == all_roots,
            "rho no longer exchanges the two complete C12 node orbits")

    f_zero, g_zero = values[0]
    require(f_zero == (340, 181) and g_zero == (327, 3),
            "specialized hidden base values changed")
    require(point_order(f_zero) == 18 and point_order(g_zero) == 6,
            "hidden base-point orders changed")
    require(g_zero == integer_multiple(15, f_zero),
            "independent integer relation g(Q0)=[15]f(Q0) changed")
    require(integer_multiple(6, f_zero) == (0, 396)
            and integer_multiple(9, f_zero) == (35, 0)
            and integer_multiple(18, f_zero) is None,
            "order-18 witnesses changed")
    return f_zero, g_zero, nodes


def verify_uniform_obstruction(
    vectors: set[Vector], representatives: list[Vector], f_point: Point,
    g_point: Point,
) -> tuple[list[str], Counter[int], str]:
    require(e_norm(PI) == 3 and e_mod_pi(PI) == 0,
            "pi=omega^2-1 lost norm/residue three")
    residue_zeros = {
        (left, right)
        for left in range(3)
        for right in range(3)
        if (left * left + right * right) % 3 == 0
    }
    require(residue_zeros == {(0, 0)},
            "A^2+B^2 acquired a nontrivial zero over F3")

    # Every hidden degree is divisible by six: the cross trace is
    # Tr((-4-2omega)(r+s*omega))=6(s-r).
    for r in range(-3, 4):
        for s in range(-3, 4):
            require(e_trace(e_mul((-4, -2), (r, s))) == 6 * (s - r),
                    "hidden cross-trace divisibility identity changed")
    require(all((4 * k_hidden) % 6 != 0
                for k_hidden in (13, 11, 10, 7, 5)),
            "uniform pi-divisibility contradiction changed")

    # Audit the uniform statement on all 3,168 vectors, before quotienting.
    for vector in vectors:
        left, right = hidden_coefficients(vector)
        hidden_q = hidden_degree(left, right)
        require(hidden_q % 12 == 0, "hidden projector degree lost 12-divisibility")
        k_hidden = hidden_q // 12
        require(k_hidden in (13, 11, 10, 7, 5) and k_hidden % 3 != 0,
                "a degree-42 residual acquired a K divisible by three")
        require(hidden_q // e_norm(PI) == 4 * k_hidden
                and (4 * k_hidden) % 6 != 0,
                "uniform divided hidden-degree contradiction changed")
        delta = e_add(e_mul(left, left), e_mul(OMEGA, e_mul(right, right)))
        delta_residue = e_mod_pi(delta)
        require(delta_residue
                == (e_mod_pi(left) ** 2 + e_mod_pi(right) ** 2) % 3,
                "Delta residue formula changed")
        require(delta_residue != 0 and e_norm(delta) % 3 != 0,
                "uniform cyclic determinant acquired pi divisibility")

        # Exact coefficient elimination of the Q0/Q1 system.
        second_left = e_mul(MINUS_OMEGA, right)
        eliminated_f = e_sub(e_mul(left, left), e_mul(right, second_left))
        eliminated_g = e_sub(e_mul(left, right), e_mul(right, left))
        require(eliminated_f == delta and eliminated_g == ZERO,
                "first cyclic determinant elimination changed")
        other_f = e_add(
            e_mul(e_mul(OMEGA, right), left), e_mul(left, second_left)
        )
        other_g = e_add(
            e_mul(e_mul(OMEGA, right), right), e_mul(left, left)
        )
        require(other_f == ZERO and other_g == delta,
                "second cyclic determinant elimination changed")

    determinant_norms: Counter[int] = Counter()
    rows: list[str] = []
    for index, vector in enumerate(representatives):
        left, right = hidden_coefficients(vector)
        delta = e_add(e_mul(left, left), e_mul(OMEGA, e_mul(right, right)))
        norm_delta = e_norm(delta)
        determinant_norms[norm_delta] += 1

        # Load-bearing order/norm implication, checked in the specialized
        # group: bar(Delta)Delta=N(Delta), so Delta*F=O would imply
        # [N(Delta)]F=O.  The direct calculations are hostile controls.
        delta_f = eisenstein_multiple(delta, f_point)
        delta_g = eisenstein_multiple(delta, g_point)
        require(delta_f == eisenstein_multiple_naive(delta, f_point),
                "double-and-add and naive CM multiplication disagreed")
        require(delta_f is not None and delta_g is not None,
                "a direct residual determinant unexpectedly annihilated F or G")
        require(eisenstein_multiple(e_conjugate(delta), delta_f)
                == integer_multiple(norm_delta, f_point),
                "bar(Delta)Delta=N(Delta) failed in the specialized CM action")
        require(integer_multiple(norm_delta, f_point) is not None,
                "a determinant norm unexpectedly annihilated the order-18 point")
        rows.append(
            f"orbit{index:02d} A={left} B={right} Delta={delta} "
            f"NDelta={norm_delta} DeltaF={delta_f} DeltaG={delta_g}"
        )

    require(determinant_norms == Counter({
        4: 2, 52: 12, 100: 4, 148: 8, 196: 8, 208: 12,
        244: 8, 292: 12, 388: 4, 400: 12, 436: 4, 592: 12,
        628: 4, 676: 14, 724: 4, 772: 4, 916: 4, 964: 4,
    }), "cyclic determinant norm distribution changed")
    row_digest = sha256("\n".join(rows).encode("ascii")).hexdigest()
    return rows, determinant_norms, row_digest


def verify_hidden_common_value_gate() -> None:
    # On the hidden plane T^2=-omega, hence T^6=-1 and T^8=omega.
    for left, right in ((ONE, ZERO), (ZERO, ONE), ((2, -1), (-3, 4))):
        current = left, right
        powers = [current]
        for _ in range(8):
            current = e_mul(MINUS_OMEGA, current[1]), current[0]
            powers.append(current)
        require(powers[6] == (e_neg(left), e_neg(right)), "T^6!=-1 on hidden plane")
        require(powers[8] == (e_mul(OMEGA, left), e_mul(OMEGA, right)),
                "T^8!=omega on hidden plane")

    # Explicit Bezout: 2(1+omega)-omega(1-omega)=1.  A common collapsed
    # hidden value killed by 2 and 1-omega is therefore O.
    bezout = e_add(
        e_scale(2, (1, 1)),
        e_mul(MINUS_OMEGA, (1, -1)),
    )
    require(bezout == ONE, "hidden common-value Bezout identity changed")


def main() -> None:
    vectors, representatives, orbit_sizes, profile, max_hidden_coordinate = (
        enumerate_degree42_residual()
    )
    verify_hidden_common_value_gate()
    f_point, g_point, nodes = verify_good_specialization()
    rows, determinant_norms, row_digest = verify_uniform_obstruction(
        vectors, representatives, f_point, g_point
    )

    print("THM-4249 degree-42 R=1/3 hidden-projection independent audit")
    print("status=FINITE-EXACT_GOOD-REDUCTION_CERTIFICATE_RELATIVE_TO_THM4230_4241_4249")
    print("engine=python_standard_library exact_integer_and_F397_arithmetic no_sampling")
    print(f"char0_profiles={dict(sorted(profile.items()))}")
    print(
        f"char0_universe=vectors:{len(vectors)} orbits:{len(representatives)} "
        f"sizes:{dict(sorted(orbit_sizes.items()))} max_hidden_coordinate:{max_hidden_coordinate}"
    )
    print("char0_box=qH>2*(N(A)+N(B));Nsum<78;coordinates_abs_at_most10;hostile_boundary_clear")
    print("projector=H=A*f+B*g A=2b+omega^2*d B=2c+d")
    print("hidden_common_value=T6:-1 T8:omega Bezout(2,1-omega):1 hence_common_value:O")
    print("cyclic_system=[A,B;-omega*B,A] determinant=Delta:A^2+omega*B^2")
    print("uniform_pi=omega^2-1 norm3 O/pi:F3 omega:1 A2+B2_zero_iff_A=B=0")
    print("uniform_degree_contradiction=pi_divides_Delta_implies_pi_divides_A_B qHprime=4K_not_in_6Z K:13,11,10,7,5")
    print("good_field=F397 affine_source_gate_smooth target_good simultaneous_Hensel_simple_roots=PASS")
    print("embedding=zeta12:157 sqrt3:377 hidden_a:161 hidden_scale:27")
    print("canonical_node=Q0:(15,28) t:379 R:1/3 nodes:12 all_regular_and_distinct")
    print("root_choices=24 two_C12_orbits rho_exchanges rho_fixes_each_of_132_map_classes")
    print(f"fQ={f_point} order:18 gQ={g_point} order:6 relation:gQ=[15]fQ")
    print("order_witnesses=[6]fQ:(0,396) [9]fQ:(35,0) [18]fQ:O")
    print("load_bearing=Delta*fQ=O_implies_[NDelta]fQ=O_implies_18_divides_NDelta_contradicts_3_not_divides_NDelta")
    print(f"cyclic_determinant_norms={dict(sorted(determinant_norms.items()))}")
    for row in rows:
        print(row)
    print(f"hostile_direct_evaluations=DeltaF_nonzero:132/132 DeltaG_nonzero:132/132 rows_sha256:{row_digest}")
    print("consequence=all_132_degree42_R_one_third_incidences_EXCLUDED")
    print("envelope_update=S42_at_most34_ratios degree42_incidences:780_to_648 total_incidences:1644_to_1512")
    print("scope=other_ratio_tests_W0_M12_seam_entry_JC2_DC2_OPEN")
    print("result=PASS")


if __name__ == "__main__":
    main()
