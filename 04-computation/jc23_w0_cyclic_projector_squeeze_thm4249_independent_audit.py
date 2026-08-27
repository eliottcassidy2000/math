#!/usr/bin/env python3
"""Clean-room exact lattice/orbit audit for THM-4249.

This standard-library path independently reconstructs the full THM-4241
Hermitian shells in the normalized O=Z[omega] basis (u,f,g,h), verifies the
C12 action and its three projector polynomials, and checks the complete
integer sieve induced by THM-4249's geometric projector bounds.

The two geometric inputs

    collapse and a != 0  => N(a) >= 7,
    collapse and d != 0  => N(d) >= 3,

and the attachment noncollapse of hidden q=12,24 shells are theorem inputs,
not re-proved by this arithmetic audit.  Their exact consequence universe,
boundary witnesses, hidden orbit representatives, and final orbit counts are
all reconstructed here.  Nothing in this file proves M=12, chart entry, or
JC(2).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction


E = tuple[int, int]                    # m+n*omega, omega^2+omega+1=0
Vector = tuple[E, E, E, E]             # coefficients of (u,f,g,h)
HiddenVector = tuple[E, E]              # coefficients of (f,g)

ZERO: E = (0, 0)
ONE: E = (1, 0)
OMEGA: E = (0, 1)
OMEGA2: E = (-1, -1)
MINUS_OMEGA: E = (0, -1)
MINUS_OMEGA2: E = (1, 1)
UNITS: tuple[E, ...] = (
    ONE, (-1, 0), OMEGA, MINUS_OMEGA, OMEGA2, MINUS_OMEGA2,
)
ZERO_VECTOR: Vector = (ZERO, ZERO, ZERO, ZERO)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def e_add(left: E, right: E) -> E:
    return left[0] + right[0], left[1] + right[1]


def e_neg(value: E) -> E:
    return -value[0], -value[1]


def e_sub(left: E, right: E) -> E:
    return left[0] - right[0], left[1] - right[1]


def e_mul(left: E, right: E) -> E:
    m, n = left
    r, s = right
    return m * r - n * s, m * s + n * r - n * s


def e_conjugate(value: E) -> E:
    return value[0] - value[1], -value[1]


def e_norm(value: E) -> int:
    answer = e_mul(value, e_conjugate(value))
    require(answer[1] == 0, "Eisenstein norm left Q")
    return answer[0]


def e_trace(value: E) -> int:
    return 2 * value[0] - value[1]


def e_scale_integer(multiplier: int, value: E) -> E:
    return multiplier * value[0], multiplier * value[1]


def vector_add(left: Vector, right: Vector) -> Vector:
    return tuple(e_add(left[j], right[j]) for j in range(4))  # type: ignore[return-value]


def vector_scale(scalar: E, vector: Vector) -> Vector:
    return tuple(e_mul(scalar, value) for value in vector)  # type: ignore[return-value]


def hidden_degree(left: E, right: E) -> int:
    """Diagonal of [[6,-4-2w],[-4-2w^2,6]]."""

    cross = e_trace(e_mul(e_mul(left, e_conjugate(right)), (-4, -2)))
    return 6 * e_norm(left) + 6 * e_norm(right) + cross


def split_coordinates(vector: Vector) -> tuple[int, int, int, E, E]:
    """Return N(a), N(d), K, and coefficients of 2*ell.

    Since 2h=v+omega^2*f+g,

      2*ell=(2b+omega^2*d)f+(2c+d)g,
      q(2*ell)=12K,
      q(m)=4N(a)+N(d)+3K.
    """

    a, b, c, d = vector
    hidden_left = e_add(e_scale_integer(2, b), e_mul(OMEGA2, d))
    hidden_right = e_add(e_scale_integer(2, c), d)
    hidden_q = hidden_degree(hidden_left, hidden_right)
    require(hidden_q % 12 == 0, "q(2ell) lost divisibility by 12")
    return e_norm(a), e_norm(d), hidden_q // 12, hidden_left, hidden_right


def full_degree(vector: Vector) -> int:
    norm_a, norm_d, k_hidden, _, _ = split_coordinates(vector)
    return 4 * norm_a + norm_d + 3 * k_hidden


def tau(vector: Vector) -> Vector:
    """Precomposition by the order-twelve source automorphism.

    T u=-omega*u, T f=g, T g=-omega*f,
    T h=omega^2*h-omega*f.
    """

    a, b, c, d = vector
    return (
        e_mul(MINUS_OMEGA, a),
        e_add(e_mul(MINUS_OMEGA, c), e_mul(MINUS_OMEGA, d)),
        b,
        e_mul(OMEGA2, d),
    )


def hidden_tau(vector: HiddenVector) -> HiddenVector:
    left, right = vector
    return e_mul(MINUS_OMEGA, right), left


def polynomial_multiply(left: list[E], right: list[E]) -> list[E]:
    answer = [ZERO] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] = e_add(answer[i + j], e_mul(a, b))
    return answer


def polynomial_scale(scalar: E, polynomial: list[E]) -> list[E]:
    return [e_mul(scalar, coefficient) for coefficient in polynomial]


def apply_tau_polynomial(polynomial: list[E], vector: Vector) -> Vector:
    answer = ZERO_VECTOR
    current = vector
    for coefficient in polynomial:
        answer = vector_add(answer, vector_scale(coefficient, current))
        current = tau(current)
    return answer


# Coefficients are in increasing powers of T.
P_U = polynomial_scale(
    (-1, 0),
    polynomial_multiply([e_neg(OMEGA2), ONE], [OMEGA, ZERO, ONE]),
)
P_V = polynomial_scale(
    e_neg(OMEGA2),
    polynomial_multiply([OMEGA, ONE], [OMEGA, ZERO, ONE]),
)
P_L = polynomial_scale(
    e_neg(OMEGA2),
    polynomial_multiply([e_neg(OMEGA2), ZERO, ONE],
                        [e_neg(OMEGA), ZERO, ONE]),
)


def expected_projectors(vector: Vector) -> tuple[Vector, Vector, Vector]:
    a, b, c, d = vector
    _, _, _, hidden_left, hidden_right = split_coordinates(vector)
    # v=2h-omega^2 f-g.
    p_u = (a, ZERO, ZERO, ZERO)
    p_v = (ZERO, e_mul(e_neg(OMEGA2), d), e_neg(d), e_scale_integer(2, d))
    p_l = (ZERO, hidden_left, hidden_right, ZERO)
    return p_u, p_v, p_l


def enumerate_shells() -> dict[int, set[Vector]]:
    # THM-4241's positive-form bounds permit these boxes.  Boundary exclusion
    # below is part of the independent exact exhaustion certificate.
    elements = [
        (m, n)
        for m in range(-8, 9)
        for n in range(-8, 9)
        if e_norm((m, n)) <= 42
    ]
    require(not any(max(abs(a), abs(b)) == 8 for a, b in elements),
            "Eisenstein element box touched boundary")

    hidden_by_parity: dict[
        tuple[int, int, int, int], list[tuple[E, E, int]]
    ] = defaultdict(list)
    hidden_rows: list[tuple[E, E, int]] = []
    for am in range(-12, 13):
        for an in range(-12, 13):
            for bm in range(-12, 13):
                for bn in range(-12, 13):
                    left, right = (am, an), (bm, bn)
                    degree = hidden_degree(left, right)
                    if degree <= 168:
                        row = (left, right, degree)
                        hidden_rows.append(row)
                        hidden_by_parity[
                            (am % 2, an % 2, bm % 2, bn % 2)
                        ].append(row)
    require(not any(
        max(abs(a), abs(b), abs(c), abs(d)) == 12
        for (a, b), (c, d), _ in hidden_rows
    ), "hidden box touched boundary")

    shells: dict[int, set[Vector]] = {34: set(), 42: set()}
    for d in elements:
        omega2_d = e_mul(OMEGA2, d)
        parity = (
            omega2_d[0] % 2, omega2_d[1] % 2,
            d[0] % 2, d[1] % 2,
        )
        for hidden_left, hidden_right, hidden_q in hidden_by_parity[parity]:
            b = (
                (hidden_left[0] - omega2_d[0]) // 2,
                (hidden_left[1] - omega2_d[1]) // 2,
            )
            c = (
                (hidden_right[0] - d[0]) // 2,
                (hidden_right[1] - d[1]) // 2,
            )
            base_degree = e_norm(d) + hidden_q // 4
            for a in elements:
                degree = base_degree + 4 * e_norm(a)
                if degree in shells:
                    vector = (a, b, c, d)
                    require(full_degree(vector) == degree,
                            "split degree reconstruction failed")
                    shells[degree].add(vector)
    require({degree: len(shell) for degree, shell in shells.items()}
            == {34: 36288, 42: 16992}, "full theta shells changed")
    return shells


def orbit_partition(vectors: set[Vector]) -> tuple[list[Vector], Counter[int]]:
    unseen = set(vectors)
    representatives: list[Vector] = []
    sizes: Counter[int] = Counter()
    while unseen:
        seed = min(unseen)
        orbit: set[Vector] = set()
        current = seed
        for _ in range(12):
            orbit.update(vector_scale(unit, current) for unit in UNITS)
            current = tau(current)
        require(current == seed, "T^12 failed in orbit builder")
        require(orbit <= vectors, "orbit left claimed invariant set")
        unseen -= orbit
        representatives.append(seed)
        sizes[len(orbit)] += 1
    return representatives, sizes


def hidden_shell_orbits(degree: int) -> tuple[set[HiddenVector], list[HiddenVector], Counter[int]]:
    vectors = {
        ((a, b), (c, d))
        for a in range(-5, 6)
        for b in range(-5, 6)
        for c in range(-5, 6)
        for d in range(-5, 6)
        if hidden_degree((a, b), (c, d)) == degree
    }
    require(not any(
        max(abs(a), abs(b), abs(c), abs(d)) == 5
        for (a, b), (c, d) in vectors
    ), f"hidden degree-{degree} box touched boundary")

    unseen = set(vectors)
    representatives: list[HiddenVector] = []
    sizes: Counter[int] = Counter()
    while unseen:
        seed = min(unseen)
        orbit: set[HiddenVector] = set()
        current = seed
        for _ in range(12):
            orbit.update((e_mul(unit, current[0]), e_mul(unit, current[1]))
                         for unit in UNITS)
            current = hidden_tau(current)
        require(current == seed, "hidden T^12 failed")
        require(orbit <= vectors, "hidden orbit left shell")
        unseen -= orbit
        representatives.append(seed)
        sizes[len(orbit)] += 1
    return vectors, representatives, sizes


def canonical_associate(value: E) -> E:
    require(value != ZERO, "zero has no maximal divisibility class")
    return min(e_mul(unit, value) for unit in UNITS)


def eisenstein_divides(divisor: E, dividend: E) -> bool:
    """Exact principal-ideal divisibility in Z[omega]."""

    denominator = e_norm(divisor)
    numerator = e_mul(dividend, e_conjugate(divisor))
    return numerator[0] % denominator == 0 and numerator[1] % denominator == 0


TorsionPoint = tuple[Fraction, Fraction]  # coefficients in Q^2/Z^2


def fraction_mod_one(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


def canonical_torsion(value: tuple[Fraction, Fraction]) -> TorsionPoint:
    return fraction_mod_one(value[0]), fraction_mod_one(value[1])


def torsion_kernel(element: E) -> set[TorsionPoint]:
    """Enumerate element^(-1)O/O in the rational O-plane.

    Multiplication by c=(m,n) has matrix [[m,-n],[n,m-n]] and determinant
    N(c).  Enumerating an N-by-N numerator box and deduplicating gives the
    full kernel of cardinality N(c), independently of elliptic coordinates.
    """

    m, n = element
    norm = e_norm(element)
    answer: set[TorsionPoint] = set()
    for left in range(norm):
        for right in range(norm):
            answer.add(canonical_torsion((
                Fraction((m - n) * left + n * right, norm),
                Fraction(-n * left + m * right, norm),
            )))
    require(len(answer) == norm, "torsion-kernel cardinality changed")
    return answer


def torsion_unit_action(unit: E, point: TorsionPoint) -> TorsionPoint:
    return canonical_torsion(e_mul(unit, point))


def torsion_orbits(points: set[TorsionPoint]) -> list[set[TorsionPoint]]:
    unseen = set(points)
    answer: list[set[TorsionPoint]] = []
    while unseen:
        seed = min(unseen)
        orbit = {torsion_unit_action(unit, seed) for unit in UNITS}
        require(orbit <= points, "torsion union is not unit invariant")
        unseen -= orbit
        answer.append(orbit)
    return answer


def format_counter(counter: Counter[object]) -> str:
    return "{" + ",".join(
        f"{key}:{counter[key]}" for key in sorted(counter)
    ) + "}"


def main() -> None:
    shells = enumerate_shells()

    # Linear action and projector checks are performed on every shell vector,
    # not just on the four coordinate generators.
    for shell in shells.values():
        for vector in shell:
            current = vector
            for _ in range(12):
                current = tau(current)
            require(current == vector, "T^12 changed a shell vector")
            require(full_degree(tau(vector)) == full_degree(vector),
                    "T changed the Hermitian degree")

            computed = (
                apply_tau_polynomial(P_U, vector),
                apply_tau_polynomial(P_V, vector),
                apply_tau_polynomial(P_L, vector),
            )
            expected = expected_projectors(vector)
            require(computed == expected, "projector formula changed")
            reconstructed = vector_add(
                vector_scale((2, 0), computed[0]),
                vector_add(computed[1], computed[2]),
            )
            require(reconstructed == vector_scale((2, 0), vector),
                    "2I=2P_u+P_v+P_L failed")

    raw_orbits: dict[int, tuple[list[Vector], Counter[int]]] = {
        degree: orbit_partition(shell) for degree, shell in shells.items()
    }
    require(len(raw_orbits[34][0]) == 648
            and raw_orbits[34][1] == Counter({72: 432, 24: 216}),
            "degree-34 raw orbit census changed")
    require(len(raw_orbits[42][0]) == 344
            and raw_orbits[42][1] == Counter({72: 186, 24: 142, 12: 16}),
            "degree-42 raw orbit census changed")

    hidden_data = {
        degree: hidden_shell_orbits(degree) for degree in (12, 24)
    }
    for degree in (12, 24):
        vectors, representatives, sizes = hidden_data[degree]
        require(len(vectors) == 24 and len(representatives) == 2
                and sizes == Counter({12: 2}),
                f"hidden degree-{degree} orbit census changed")

    # The geometric projector bounds plus the hidden q=12,24 no-go force a=0.
    # We retain the nonempty arithmetic boundary to show the hidden no-go is
    # load-bearing rather than silently assuming an empty degree equation.
    expected_a_boundary = {
        34: Counter({(7, 3, 1): 576, (7, 0, 2): 288}),
        42: Counter({(9, 3, 1): 288, (9, 0, 2): 144}),
    }
    boundary_representatives: dict[int, dict[tuple[int, int, int], Vector]] = {}
    for degree, shell in shells.items():
        boundary: Counter[tuple[int, int, int]] = Counter()
        witnesses: dict[tuple[int, int, int], Vector] = {}
        for vector in sorted(shell):
            norm_a, norm_d, k_hidden, _, _ = split_coordinates(vector)
            if norm_a >= 7 and (norm_d == 0 or norm_d >= 3):
                key = (norm_a, norm_d, k_hidden)
                boundary[key] += 1
                witnesses.setdefault(key, vector)
        require(boundary == expected_a_boundary[degree],
                f"degree-{degree} nonzero-a boundary changed")
        require(all(key[2] in (1, 2) for key in boundary),
                "nonzero-a branch escaped hidden q12/q24 no-go")
        boundary_representatives[degree] = witnesses

    stage_counts: dict[int, tuple[int, int, int, int]] = {}
    final_sets: dict[int, set[Vector]] = {}
    final_profiles: dict[int, Counter[tuple[int, int]]] = {}
    for degree, shell in shells.items():
        a_zero = {vector for vector in shell if vector[0] == ZERO}
        d_bound = {
            vector for vector in a_zero
            if vector[3] == ZERO or e_norm(vector[3]) >= 3
        }
        hidden_no_go = {
            vector for vector in d_bound
            if split_coordinates(vector)[2] not in (1, 2)
        }
        pure_hidden = {
            vector for vector in hidden_no_go
            if not (degree == 42 and vector[3] == ZERO)
        }
        stage_counts[degree] = (
            len(a_zero), len(d_bound), len(hidden_no_go), len(pure_hidden)
        )
        final_sets[degree] = pure_hidden
        final_profiles[degree] = Counter(
            (split_coordinates(vector)[1], split_coordinates(vector)[2])
            for vector in pure_hidden
        )

    require(stage_counts == {
        34: (5184, 4608, 4224, 4224),
        42: (3600, 3600, 3360, 3168),
    }, "projector sieve stage counts changed")
    require(final_profiles == {
        34: Counter({
            (4, 10): 864, (7, 9): 1248, (13, 7): 768,
            (16, 6): 576, (19, 5): 576, (25, 3): 192,
        }),
        42: Counter({
            (3, 13): 672, (9, 11): 576, (12, 10): 864,
            (21, 7): 768, (27, 5): 288,
        }),
    }, "final (N(d),K) profile changed")

    final_orbits = {
        degree: orbit_partition(vectors)
        for degree, vectors in final_sets.items()
    }
    require(len(final_orbits[34][0]) == 176
            and final_orbits[34][1] == Counter({24: 176}),
            "degree-34 final orbit census changed")
    require(len(final_orbits[42][0]) == 132
            and final_orbits[42][1] == Counter({24: 132}),
            "degree-42 final orbit census changed")

    # Torsion envelopes.  Evaluating P_v on a collapsed C12 orbit gives
    # d*v(Q)=-P, while comparison at tau Q forces (omega^2-1)P=0.  Hence
    # v(Q) lies in E[c], c=d*pi, pi=omega^2-1.  On the explicit visible map,
    # R=U/Z=-1/(X(v(Q))^3+1), and X^3 is the exact quotient by target mu_6.
    pi = e_sub(OMEGA2, ONE)
    require(e_norm(pi) == 3, "pi=omega^2-1 lost norm three")

    maximal_d: dict[int, set[E]] = {}
    maximal_c: dict[int, list[E]] = {}
    expected_c_norms = {
        34: Counter({21: 2, 39: 2, 48: 1, 57: 2, 75: 1}),
        42: Counter({36: 1, 63: 2, 81: 1}),
    }
    for degree, vectors in final_sets.items():
        d_classes = {canonical_associate(vector[3]) for vector in vectors}
        maximal = {
            d for d in d_classes
            if not any(d != other and eisenstein_divides(d, other)
                       for other in d_classes)
        }
        require(all(any(eisenstein_divides(d, top) for top in maximal)
                    for d in d_classes),
                "a residual d-class is not nested in a maximal class")
        maximal_d[degree] = maximal
        c_values = sorted(canonical_associate(e_mul(pi, d)) for d in maximal)
        maximal_c[degree] = c_values
        require(Counter(e_norm(c) for c in c_values) == expected_c_norms[degree],
                f"degree-{degree} maximal c norms changed")

    torsion_unions: dict[int, set[TorsionPoint]] = {}
    torsion_orbit_data: dict[int, tuple[list[set[TorsionPoint]], Counter[int]]] = {}
    intersection_sizes: dict[int, Counter[int]] = {}
    for degree, c_values in maximal_c.items():
        kernels = [torsion_kernel(c) for c in c_values]
        union = set().union(*kernels)
        intersections = Counter(
            len(kernels[i] & kernels[j])
            for i in range(len(kernels))
            for j in range(i + 1, len(kernels))
        )
        intersection_sizes[degree] = intersections
        torsion_unions[degree] = union
        orbits = torsion_orbits(union)
        torsion_orbit_data[degree] = (orbits, Counter(map(len, orbits)))

    require(len(torsion_unions[34]) == 336
            and intersection_sizes[34] == Counter({3: 28}),
            "degree-34 torsion union/intersections changed")
    require(torsion_orbit_data[34][1]
            == Counter({6: 55, 1: 1, 2: 1, 3: 1}),
            "degree-34 torsion orbit profile changed")
    require(len(torsion_unions[42]) == 216
            and intersection_sizes[42] == Counter({9: 6}),
            "degree-42 torsion union/intersections changed")
    require(torsion_orbit_data[42][1]
            == Counter({6: 35, 1: 1, 2: 1, 3: 1}),
            "degree-42 torsion orbit profile changed")

    kernel_pi = torsion_kernel(pi)
    kernel_two = torsion_kernel((2, 0))
    kernel_three = torsion_kernel((3, 0))
    require(len(kernel_pi) == 3 and len(kernel_two) == 4
            and len(kernel_three) == 9,
            "small torsion kernels changed")
    require(set.intersection(*(torsion_kernel(c) for c in maximal_c[34]))
            == kernel_pi,
            "degree-34 common intersection is not E[pi]")
    require(set.intersection(*(torsion_kernel(c) for c in maximal_c[42]))
            == kernel_three,
            "degree-42 common intersection is not E[3]")
    require(kernel_pi <= kernel_three, "E[pi] not contained in E[3]")

    # O is nonaffine; E[pi]\{O} has X^3=0 and R=-1; E[2]\{O}
    # has X^3=-1 and corresponds to Z=0.  Removing these three exceptional
    # unit orbits leaves the claimed finite ratio envelopes.
    excluded_small = kernel_pi | kernel_two
    ratio_orbits: dict[int, list[set[TorsionPoint]]] = {
        degree: [orbit for orbit in torsion_orbit_data[degree][0]
                 if orbit.isdisjoint(excluded_small)]
        for degree in (34, 42)
    }
    require(len(ratio_orbits[34]) == 55
            and Counter(map(len, ratio_orbits[34])) == Counter({6: 55}),
            "degree-34 55-ratio envelope changed")
    require(len(ratio_orbits[42]) == 35
            and Counter(map(len, ratio_orbits[42])) == Counter({6: 35}),
            "degree-42 35-ratio envelope changed")

    three_extra = kernel_three - kernel_pi
    require(len(three_extra) == 6
            and any(three_extra == orbit for orbit in ratio_orbits[42]),
            "unique E[3] extra orbit missing from degree-42 envelope")
    # psi_3=3X(X^3+4): on E[3]\E[pi], X^3=-4, hence
    # R=-1/(X^3+1)=1/3 exactly.
    e3_x_cubed = -4
    e3_ratio = Fraction(-1, e3_x_cubed + 1)
    require(e3_ratio == Fraction(1, 3), "E[3] overlap ratio changed")

    # Retaining each map's actual d, rather than replacing it by a maximal
    # containing ideal, makes the decisive attachment workload much smaller.
    # Source C12 and target mu_6 have already been quotiented in final_orbits.
    incidence_profiles: dict[int, Counter[tuple[int, int]]] = {}
    incidence_totals: dict[int, int] = {}
    for degree in (34, 42):
        profile: Counter[tuple[int, int]] = Counter()
        total = 0
        for vector in final_orbits[degree][0]:
            norm_d = e_norm(vector[3])
            c = e_mul(pi, vector[3])
            candidate_count = sum(
                orbit.isdisjoint(excluded_small)
                for orbit in torsion_orbits(torsion_kernel(c))
            )
            profile[(norm_d, candidate_count)] += 1
            total += candidate_count
        incidence_profiles[degree] = profile
        incidence_totals[degree] = total
    require(incidence_profiles == {
        34: Counter({
            (4, 1): 36, (7, 3): 52, (13, 6): 32,
            (16, 7): 24, (19, 9): 24, (25, 12): 8,
        }),
        42: Counter({
            (3, 1): 28, (9, 4): 24, (12, 5): 36,
            (21, 10): 32, (27, 13): 12,
        }),
    }, "map-specific torsion incidence profile changed")
    require(incidence_totals == {34: 864, 42: 780},
            "map-specific attachment workload changed")

    # For fixed nonzero A^6,B^4 there are 24 root choices.  The diagonal
    # C12 orbit consists exactly of exponent pairs of equal parity.  The
    # source involution rho:(x,y)->(-x,-y), exponent shift (3,2), exchanges
    # it with the other orbit.  Since precomposition by rho permutes each
    # full Hom shell, testing every map on one canonical node orbit is enough.
    all_root_choices = {(left, right) for left in range(6) for right in range(4)}
    canonical_node_orbit = {(j % 6, j % 4) for j in range(12)}
    rho_shifted_orbit = {((left + 3) % 6, (right + 2) % 4)
                         for left, right in canonical_node_orbit}
    require(len(canonical_node_orbit) == 12
            and canonical_node_orbit.isdisjoint(rho_shifted_orbit)
            and canonical_node_orbit | rho_shifted_orbit == all_root_choices,
            "rho no longer exchanges the two attachment C12 orbits")

    print("THM-4249 W=0 cyclic-projector squeeze independent audit")
    print("arithmetic=exact_standard_library no_floating_point no_sampling")
    print("basis=[u,f,g,h] relation=2h=v+omega^2*f+g")
    print("T_action=u:-omega*u;f:g;g:-omega*f;h:omega^2*h-omega*f")
    print(f"projector_coefficients Pu={P_U} Pv={P_V} PL={P_L}")
    print("projectors=Pu:m->a*u Pv:m->d*v PL:m->2*ell PASS_on_all_53280_vectors")
    print("reconstruction=2I=2Pu+Pv+PL T_order=12 degree_invariant=PASS")
    print("degree_identity=q(m)=4N(a)+N(d)+3K q(2ell)=12K")
    print("raw_shells=d34:36288 d42:16992")
    print("raw_orbits=d34:648[sizes24:216,72:432] d42:344[sizes12:16,24:142,72:186]")
    for degree in (12, 24):
        vectors, representatives, sizes = hidden_data[degree]
        print(
            f"hidden_q{degree}=vectors:{len(vectors)} orbits:{len(representatives)} "
            f"sizes:{format_counter(sizes)} reps:{representatives}"
        )
    print(
        "geometric_inputs=collapse_nonzero_a_implies_Na>=7;"
        "collapse_nonzero_d_implies_Nd>=3;hidden_q12_q24_nocollapse"
    )
    for degree in (34, 42):
        print(
            f"hostile_nonzero_a_d{degree}="
            f"profiles:{format_counter(expected_a_boundary[degree])} "
            f"reps:{boundary_representatives[degree]}"
        )
    print("a_zero_raw=d34:5184 d42:3600")
    print("sieve_stages=d34:5184->4608->4224->4224 d42:3600->3600->3360->3168")
    print(
        f"final_profiles_d34={format_counter(final_profiles[34])} "
        f"final_profiles_d42={format_counter(final_profiles[42])}"
    )
    print("final_orbits=d34:176_all_size24 d42:132_all_size24 total:308")
    print("pure_hidden_degree42=192_vectors_16_raw_orbits removed_by_THM4230")
    print(f"torsion_pi=omega^2-1:{pi} norm:{e_norm(pi)} ratio_formula=R:-1/(X^3+1)")
    print(
        f"maximal_c_d34={[(e_norm(c), c) for c in maximal_c[34]]} "
        f"maximal_c_d42={[(e_norm(c), c) for c in maximal_c[42]]}"
    )
    print("ideal_nesting=every_residual_d_class_divides_a_listed_maximal_d PASS")
    print(
        f"torsion_union_d34=points:{len(torsion_unions[34])} "
        f"pair_intersections:{format_counter(intersection_sizes[34])} "
        f"unit_orbits:{format_counter(torsion_orbit_data[34][1])}"
    )
    print(
        f"torsion_union_d42=points:{len(torsion_unions[42])} "
        f"pair_intersections:{format_counter(intersection_sizes[42])} "
        f"unit_orbits:{format_counter(torsion_orbit_data[42][1])}"
    )
    print("common_intersections=d34:E[pi]_size3 d42:E[3]_size9")
    print("excluded_ratios=E[pi]_nonzero:R=-1;E[2]_nonzero:Z=0;origin:nonaffine")
    print("E3_extra=one_mu6_orbit_size6 psi3:3X(X^3+4) R=1/3")
    print("torsion_ratio_envelopes=S34_at_most55 S42_at_most35")
    print(
        f"map_specific_ratio_counts_d34={format_counter(incidence_profiles[34])} "
        f"map_specific_ratio_counts_d42={format_counter(incidence_profiles[42])}"
    )
    print("node_root_orbits=24_choices_two_C12_orbits rho_shift_(3,2)_exchanges PASS")
    print("decisive_attachment_incidences=d34:864 d42:780 total:1644 one_canonical_node_orbit")
    print("scope=conditional_arithmetic_audit_of_named_geometric_filters;M12_entry_JC2_OPEN")
    print("result=PASS")


if __name__ == "__main__":
    main()
