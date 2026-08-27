#!/usr/bin/env python3
"""Exact certificate for THM-4249's W=0 cyclic-projector squeeze.

This is an independent follow-up to THM-4230 and THM-4241.  It uses the
normalized O-basis (u,f,g,h) of M=Hom(J(C0),E0), where

    O = Z[omega],  omega^2+omega+1=0,
    h = (v + omega^2*f + g)/2,

u and v are the two visible degree-four eigenmaps, and O*f+O*g is the
hidden lattice.  The script verifies the integral tau projectors, enumerates
the exact degree-34/42 residual after their geometric fibre constraints, and
classifies the low-degree hidden vectors used at the last exclusion step.

It does not evaluate the remaining mixed attachment maps.  In particular it
does not prove S_34 or S_42 empty, close M=12, or prove JC(2).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product

import sympy as sp


E = tuple[int, int]  # m+n*omega
Vector = tuple[E, E, E, E]  # coefficients in (u,f,g,h)
Matrix = tuple[tuple[E, ...], ...]

ZERO: E = (0, 0)
ONE: E = (1, 0)
OMEGA: E = (0, 1)
OMEGA2: E = (-1, -1)
MINUS_OMEGA: E = (0, -1)
UNITS: tuple[E, ...] = (
    (1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def e_add(left: E, right: E) -> E:
    return left[0] + right[0], left[1] + right[1]


def e_neg(value: E) -> E:
    return -value[0], -value[1]


def e_sub(left: E, right: E) -> E:
    return e_add(left, e_neg(right))


def e_mul(left: E, right: E) -> E:
    m, n = left
    r, s = right
    return m * r - n * s, m * s + n * r - n * s


def e_scale(multiplier: int, value: E) -> E:
    return multiplier * value[0], multiplier * value[1]


def e_conjugate(value: E) -> E:
    m, n = value
    return m - n, -n


def e_trace(value: E) -> int:
    return 2 * value[0] - value[1]


def e_norm(value: E) -> int:
    return value[0] ** 2 - value[0] * value[1] + value[1] ** 2


def e_half(value: E) -> E:
    require(value[0] % 2 == 0 and value[1] % 2 == 0,
            "nonintegral half in the glued basis")
    return value[0] // 2, value[1] // 2


def matrix_identity(size: int) -> Matrix:
    return tuple(tuple(ONE if i == j else ZERO for j in range(size))
                 for i in range(size))


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(e_add(left[i][j], right[i][j])
                       for j in range(len(left)))
                 for i in range(len(left)))


def matrix_neg(matrix: Matrix) -> Matrix:
    return tuple(tuple(e_neg(value) for value in row) for row in matrix)


def matrix_sub(left: Matrix, right: Matrix) -> Matrix:
    return matrix_add(left, matrix_neg(right))


def matrix_scale(scalar: E, matrix: Matrix) -> Matrix:
    return tuple(tuple(e_mul(scalar, value) for value in row)
                 for row in matrix)


def matrix_mul(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return tuple(tuple(
        sum_e(e_mul(left[i][k], right[k][j]) for k in range(size))
        for j in range(size)) for i in range(size))


def matrix_power(matrix: Matrix, exponent: int) -> Matrix:
    answer = matrix_identity(len(matrix))
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = matrix_mul(answer, base)
        base = matrix_mul(base, base)
        remaining //= 2
    return answer


def matrix_transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[j][i] for j in range(len(matrix)))
                 for i in range(len(matrix)))


def matrix_conjugate(matrix: Matrix) -> Matrix:
    return tuple(tuple(e_conjugate(value) for value in row)
                 for row in matrix)


def sum_e(values) -> E:
    answer = ZERO
    for value in values:
        answer = e_add(answer, value)
    return answer


# Columns are the images of u,f,g,h under precomposition by tau.
T: Matrix = (
    (MINUS_OMEGA, ZERO, ZERO, ZERO),
    (ZERO, ZERO, MINUS_OMEGA, MINUS_OMEGA),
    (ZERO, ONE, ZERO, ZERO),
    (ZERO, ZERO, ZERO, OMEGA2),
)
I4 = matrix_identity(4)
T2 = matrix_mul(T, T)

# Hermitian Gram from THM-4241, linear in the first slot.
H: Matrix = (
    ((4, 0), ZERO, ZERO, ZERO),
    (ZERO, (6, 0), (-4, -2), (-2, 2)),
    (ZERO, (-2, 2), (6, 0), (2, -2)),
    (ZERO, (-4, -2), (4, 2), (4, 0)),
)


def verify_tau_and_projectors() -> None:
    require(matrix_power(T, 12) == I4, "T^12 changed")
    trace_matrix = tuple(tuple(ZERO for _ in range(4)) for _ in range(4))
    for exponent in range(12):
        trace_matrix = matrix_add(trace_matrix, matrix_power(T, exponent))
    require(trace_matrix == tuple(tuple(ZERO for _ in range(4))
                                  for _ in range(4)),
            "1+T+...+T^11 is not zero")

    # T is Rosati-unitary: for a form linear in the first slot the identity is
    # T^t H conjugate(T)=H.
    pulled_gram = matrix_mul(matrix_mul(matrix_transpose(T), H),
                            matrix_conjugate(T))
    require(pulled_gram == H, "T lost the Rosati form")

    # Integral spectral projectors.  Applied to
    # m=a*u+b*f+c*g+d*h, these return a*u, d*v, and 2*ell respectively,
    # where ell=b*f+c*g+d(omega^2*f+g)/2.
    p_u = matrix_neg(matrix_mul(
        matrix_sub(T, matrix_scale(OMEGA2, I4)),
        matrix_add(T2, matrix_scale(OMEGA, I4)),
    ))
    p_v = matrix_scale(e_neg(OMEGA2), matrix_mul(
        matrix_add(T, matrix_scale(OMEGA, I4)),
        matrix_add(T2, matrix_scale(OMEGA, I4)),
    ))
    p_hidden = matrix_scale(e_neg(OMEGA2), matrix_mul(
        matrix_sub(T2, matrix_scale(OMEGA2, I4)),
        matrix_sub(T2, matrix_scale(OMEGA, I4)),
    ))

    expected_u: Matrix = (
        (ONE, ZERO, ZERO, ZERO),
        (ZERO, ZERO, ZERO, ZERO),
        (ZERO, ZERO, ZERO, ZERO),
        (ZERO, ZERO, ZERO, ZERO),
    )
    expected_v: Matrix = (
        (ZERO, ZERO, ZERO, ZERO),
        (ZERO, ZERO, ZERO, e_neg(OMEGA2)),
        (ZERO, ZERO, ZERO, (-1, 0)),
        (ZERO, ZERO, ZERO, (2, 0)),
    )
    expected_hidden: Matrix = (
        (ZERO, ZERO, ZERO, ZERO),
        (ZERO, (2, 0), ZERO, OMEGA2),
        (ZERO, ZERO, (2, 0), ONE),
        (ZERO, ZERO, ZERO, ZERO),
    )
    require(p_u == expected_u, "u projector changed")
    require(p_v == expected_v, "v projector changed")
    require(p_hidden == expected_hidden, "hidden projector changed")


CROSS: E = (-4, -2)


def hidden_degree(a: E, b: E) -> int:
    return (
        6 * e_norm(a) + 6 * e_norm(b)
        + e_trace(e_mul(e_mul(a, e_conjugate(b)), CROSS))
    )


def hidden_tau(vector: tuple[E, E]) -> tuple[E, E]:
    a, b = vector
    return e_mul(MINUS_OMEGA, b), a


def hidden_determinant(vector: tuple[E, E]) -> E:
    # det [[a,-omega*b],[b,a]]=a^2+omega*b^2.
    a, b = vector
    return e_add(e_mul(a, a), e_mul(OMEGA, e_mul(b, b)))


def unit_scale_vector(unit: E, vector: Vector) -> Vector:
    return tuple(e_mul(unit, coefficient) for coefficient in vector)  # type: ignore[return-value]


def tau_vector(vector: Vector) -> Vector:
    a, b, c, d = vector
    return (
        e_mul(MINUS_OMEGA, a),
        e_mul(MINUS_OMEGA, e_add(c, d)),
        b,
        e_mul(OMEGA2, d),
    )


def symmetry_orbits(vectors: set[Vector]) -> tuple[list[Vector], Counter[int]]:
    unseen = set(vectors)
    representatives: list[Vector] = []
    sizes: Counter[int] = Counter()
    while unseen:
        seed = min(unseen)
        orbit: set[Vector] = set()
        current = seed
        for _ in range(12):
            orbit.update(unit_scale_vector(unit, current) for unit in UNITS)
            current = tau_vector(current)
        require(orbit <= vectors, "unit/tau orbit left the residual")
        representatives.append(min(orbit))
        sizes[len(orbit)] += 1
        unseen -= orbit
    representatives.sort()
    return representatives, sizes


def classify_hidden_low_degrees() -> dict[int, tuple[int, int, Counter[int]]]:
    result: dict[int, tuple[int, int, Counter[int]]] = {}
    vectors_by_degree: dict[int, set[tuple[E, E]]] = {}
    for degree in (6, 12, 24):
        vectors = {
            ((a, b), (c, d))
            for a, b, c, d in product(range(-5, 6), repeat=4)
            if hidden_degree((a, b), (c, d)) == degree
        }
        require(not any(max(abs(x) for coefficient in vector for x in coefficient) == 5
                        for vector in vectors),
                "hidden low-degree enumeration touched its boundary")
        unseen = set(vectors)
        orbit_count = 0
        while unseen:
            seed = min(unseen)
            tau_seed = hidden_tau(seed)
            orbit = {
                (e_mul(unit, item[0]), e_mul(unit, item[1]))
                for unit in UNITS for item in (seed, tau_seed)
            }
            require(orbit <= vectors, "hidden low-degree orbit escaped")
            unseen -= orbit
            orbit_count += 1
        determinants = Counter(e_norm(hidden_determinant(vector))
                               for vector in vectors)
        vectors_by_degree[degree] = vectors
        result[degree] = len(vectors), orbit_count, determinants

    require(result[6] == (24, 2, Counter({1: 24})),
            "degree-six hidden classification changed")
    require(result[12] == (24, 2, Counter({4: 24})),
            "degree-twelve hidden classification changed")
    require(result[24] == (24, 2, Counter({16: 24})),
            "degree-twenty-four hidden classification changed")
    doubled_six = {
        (e_scale(2, vector[0]), e_scale(2, vector[1]))
        for vector in vectors_by_degree[6]
    }
    require(vectors_by_degree[24] == doubled_six,
            "degree 24 is no longer exactly twice degree 6")
    twice_units = {e_scale(2, unit) for unit in UNITS}
    require({hidden_determinant(vector) for vector in vectors_by_degree[12]}
            <= twice_units,
            "degree-twelve determinant is not associate to 2")
    require({hidden_determinant(vector) for vector in vectors_by_degree[6]}
            <= set(UNITS),
            "degree-six cyclic determinant is not a unit")
    return result


def enumerate_full_residual():
    # q_hidden >= (6-2sqrt(3))(N(A)+N(B)).  A coordinate of absolute value
    # 12 has Eisenstein norm at least 108, so the box below is rigorous for
    # q_hidden<=168.  Retained boundary emptiness is an exact hostile check.
    hidden_by_parity: dict[tuple[int, int, int, int], list[tuple[E, E, int]]] = defaultdict(list)
    for am, an, bm, bn in product(range(-12, 13), repeat=4):
        a_hidden = (am, an)
        b_hidden = (bm, bn)
        degree = hidden_degree(a_hidden, b_hidden)
        if degree <= 168:
            require(degree >= 0, "hidden form lost positivity")
            parity = am % 2, an % 2, bm % 2, bn % 2
            hidden_by_parity[parity].append((a_hidden, b_hidden, degree))
    retained_hidden = [row for rows in hidden_by_parity.values() for row in rows]
    require(not any(max(abs(x) for coefficient in row[:2] for x in coefficient) == 12
                    for row in retained_hidden),
            "hidden enumeration touched coordinate 12")

    # N(m+n*omega)>=3 max(|m|,|n|)^2/4 after minimizing the other
    # coordinate, so coordinate 8 is already outside N<=42.
    elements = {
        (m, n) for m, n in product(range(-8, 9), repeat=2)
        if e_norm((m, n)) <= 42
    }
    require(not any(max(abs(x) for x in value) == 8 for value in elements),
            "Eisenstein norm enumeration touched coordinate 8")
    norm_multiplicity = Counter(e_norm(value) for value in elements)

    total_counts: dict[int, int] = {}
    all_triples: dict[int, Counter[tuple[int, int, int]]] = {}
    a_zero_vectors: dict[int, set[Vector]] = {}
    for target in (34, 42):
        total = 0
        triples: Counter[tuple[int, int, int]] = Counter()
        zero_vectors: set[Vector] = set()
        for d in elements:
            norm_d = e_norm(d)
            omega2_d = e_mul(OMEGA2, d)
            parity = (omega2_d[0] % 2, omega2_d[1] % 2,
                      d[0] % 2, d[1] % 2)
            for a_hidden, b_hidden, degree_hidden in hidden_by_parity[parity]:
                require(degree_hidden % 4 == 0,
                        "glue parity lost degree divisibility by four")
                require(degree_hidden % 6 == 0,
                        "hidden degree lost divisibility by six")
                hidden_quarter = degree_hidden // 4
                base = norm_d + hidden_quarter
                remainder = target - base
                if remainder < 0 or remainder % 4:
                    continue
                norm_a = remainder // 4
                multiplicity = norm_multiplicity[norm_a]
                if not multiplicity:
                    continue
                require(degree_hidden % 12 == 0,
                        "actual glued row lost hidden divisibility by twelve")
                k_hidden = degree_hidden // 12
                total += multiplicity
                triples[(norm_a, norm_d, k_hidden)] += multiplicity
                if norm_a == 0:
                    b = e_half(e_sub(a_hidden, omega2_d))
                    c = e_half(e_sub(b_hidden, d))
                    zero_vectors.add((ZERO, b, c, d))
        total_counts[target] = total
        all_triples[target] = triples
        a_zero_vectors[target] = zero_vectors

    require(total_counts == {34: 36288, 42: 16992},
            "full THM-4241 theta counts changed")
    require(len(a_zero_vectors[34]) == 5184, "raw degree-34 a=0 count changed")
    require(len(a_zero_vectors[42]) == 3600, "raw degree-42 a=0 count changed")

    # If a projected map collapses the twelve-point attachment orbit, then:
    #   a!=0 => N(a)>=7 (six distinct nonzero u-values plus O in ker[a]);
    #   d!=0 => N(d)>=3 (degree 4N(d) has a twelve-point fibre);
    #   ell!=0 => K>=1, where q(2ell)=12K.
    # The only actual nonzero-a arithmetic profiles are therefore the four
    # displayed below.  Low-hidden-degree classification removes all four.
    possible_nonzero_a: dict[int, set[tuple[int, int, int]]] = {}
    for target in (34, 42):
        possible_nonzero_a[target] = {
            triple for triple in all_triples[target]
            if triple[0] >= 7 and (triple[1] == 0 or triple[1] >= 3)
            and triple[2] >= 1
        }
    require(possible_nonzero_a[34] == {(7, 3, 1), (7, 0, 2)},
            "degree-34 nonzero-u profiles changed")
    require(possible_nonzero_a[42] == {(9, 3, 1), (9, 0, 2)},
            "degree-42 nonzero-u profiles changed")

    # Apply the d-projection fibre bound, and inherit THM-4230's exact
    # exclusion of the 192 pure-hidden degree-42 vectors.
    residual: dict[int, set[Vector]] = {}
    def survives_projected_fibres(vector: Vector) -> bool:
        _, b, c, d = vector
        if e_norm(d) < 3:
            return False
        doubled_hidden_a = e_add(e_scale(2, b), e_mul(OMEGA2, d))
        doubled_hidden_b = e_add(e_scale(2, c), d)
        k_hidden = hidden_degree(doubled_hidden_a, doubled_hidden_b) // 12
        # The exact low-degree hidden audit excludes K=1 and K=2.
        return k_hidden >= 3

    residual[34] = {
        vector for vector in a_zero_vectors[34]
        if survives_projected_fibres(vector)
    }
    residual[42] = {
        vector for vector in a_zero_vectors[42]
        if survives_projected_fibres(vector)
    }
    require(len(residual[34]) == 4224, "degree-34 spectral residual changed")
    require(len(residual[42]) == 3168, "degree-42 spectral residual changed")

    distributions: dict[int, Counter[tuple[int, int]]] = {}
    expected_distributions = {
        34: Counter({
            (4, 10): 864, (7, 9): 1248, (13, 7): 768,
            (16, 6): 576, (19, 5): 576, (25, 3): 192,
        }),
        42: Counter({
            (3, 13): 672, (9, 11): 576, (12, 10): 864,
            (21, 7): 768, (27, 5): 288,
        }),
    }
    for target in (34, 42):
        distribution: Counter[tuple[int, int]] = Counter()
        for _, b, c, d in residual[target]:
            doubled_hidden_a = e_add(e_scale(2, b), e_mul(OMEGA2, d))
            doubled_hidden_b = e_add(e_scale(2, c), d)
            degree_hidden = hidden_degree(doubled_hidden_a, doubled_hidden_b)
            require(degree_hidden % 12 == 0, "residual hidden degree changed")
            distribution[(e_norm(d), degree_hidden // 12)] += 1
        require(distribution == expected_distributions[target],
                f"degree-{target} residual distribution changed")
        distributions[target] = distribution

    orbit_data = {}
    for target in (34, 42):
        representatives, sizes = symmetry_orbits(residual[target])
        orbit_data[target] = representatives, sizes
    require((len(orbit_data[34][0]), orbit_data[34][1])
            == (176, Counter({24: 176})),
            "degree-34 orbit quotient changed")
    require((len(orbit_data[42][0]), orbit_data[42][1])
            == (132, Counter({24: 132})),
            "degree-42 orbit quotient changed")

    # Count the remaining coefficient ideals d modulo target units.  Source T
    # changes d only by the unit omega^2, so these are also orbit invariants.
    associate_counts: dict[int, Counter[int]] = {}
    for target in (34, 42):
        d_values = {vector[3] for vector in residual[target]}
        counts: Counter[int] = Counter()
        while d_values:
            seed = min(d_values)
            associates = {e_mul(unit, seed) for unit in UNITS}
            require(associates <= d_values, "d associate class incomplete")
            counts[e_norm(seed)] += 1
            d_values -= associates
        associate_counts[target] = counts
    require(associate_counts[34] == Counter({
        4: 1, 7: 2, 13: 2, 16: 1, 19: 2,
        25: 1,
    }), "degree-34 d-ideal classes changed")
    require(associate_counts[42] == Counter({
        3: 1, 9: 1, 12: 1, 21: 2, 27: 1,
    }), "degree-42 d-ideal classes changed")

    return (total_counts, a_zero_vectors, residual, distributions,
            orbit_data, associate_counts, possible_nonzero_a)


def verify_hidden_boundary_and_ratio() -> int:
    parameter = sp.symbols("a")
    quartic = parameter**4 - 2 * parameter**3 - 2 * parameter + 1
    resultant = int(sp.resultant(quartic, parameter**6 - 1, parameter))
    require(resultant == -108, "hidden two-torsion boundary resultant changed")

    # For f and Tf at a regular attachment (y!=0), simultaneous two-torsion
    # forces t^2=-a^3 and 1+a^3*t^2=0, hence a^6=1.  The resultant excludes
    # this.  The discarded alternative y=0 is exactly the gate endpoint Z=0.
    # This proves the degree-12/24 hidden projector maps cannot collapse on
    # U*Z*(U+Z)!=0.

    # In the remaining degree-42 N(d)=3 lane, d(omega^2-1) is associate to 3.
    # The full nonzero 3-torsion factor is X*(X^3+4).  The X=0 orbit gives
    # R=-1 and is excluded by U+Z!=0.  On the other orbit, for
    # v(Q)=(-A^-2, i B^2/A^3), one has R=A^6/B^4=-1/(X^3+1)=1/3.
    x_coordinate = sp.symbols("X")
    division_three = 3 * x_coordinate * (x_coordinate**3 + 4)
    require(sp.Poly(division_three, x_coordinate).degree() == 4,
            "three-division polynomial changed")
    r = sp.symbols("r")
    marked_ratio = sp.factor(-1 / (r + 1))
    require(marked_ratio.subs(r, 0) == -1,
            "X=0 three-torsion gate ratio changed")
    require(marked_ratio.subs(r, -4) == sp.Rational(1, 3),
            "N(d)=3 marked ratio changed")
    return resultant


def verify_torsion_envelopes() -> tuple[int, int]:
    """Count the exact marked-ratio envelopes after ideal nesting.

    Put pi=omega^2-1.  Every surviving point lies in E[d*pi].  At degree 34,
    E[2*pi] is contained in E[4*pi], leaving maximal ideal norms
    21,21,39,39,48,57,57,75.  Distinct maximal kernels meet in E[pi].

    At degree 42, the pi^2 and pi^3 kernels are contained in E[pi^4],
    leaving norms 36,63,63,81.  Their pairwise intersection is E[pi^2].
    These containments use the PID O=Z[omega] and reverse divisibility of
    annihilator ideals into inclusions of torsion kernels.

    The map P=(X,Y) |-> X^3 identifies exactly the target-mu_6 orbits.
    E[pi] consists of O and the two X=0 points (R=-1, excluded).  If 2
    divides the annihilator, the three nonzero E[2] points have X^3=-1 and
    give the excluded Z=0 endpoint.  All remaining orbits have size six.
    """

    degree34_maximal_norms = (21, 21, 39, 39, 48, 57, 57, 75)
    degree34_contributions = tuple(
        (norm - (6 if norm == 48 else 3)) // 6
        for norm in degree34_maximal_norms
    )
    require(degree34_contributions == (3, 3, 6, 6, 7, 9, 9, 12),
            "degree-34 torsion-envelope contributions changed")
    # All intersections are E[pi], whose points are excluded, so the
    # admissible ratio sets are disjoint.
    degree34_envelope = sum(degree34_contributions)
    require(degree34_envelope == 55, "degree-34 torsion envelope changed")

    degree42_maximal_norms = (36, 63, 63, 81)
    degree42_contributions = (
        (36 - 6) // 6,
        (63 - 3) // 6,
        (63 - 3) // 6,
        (81 - 3) // 6,
    )
    require(degree42_contributions == (5, 10, 10, 13),
            "degree-42 torsion-envelope contributions changed")
    # The four sets share exactly E[pi^2]=E[3].  After the excluded E[pi]
    # orbit is removed, their common part is the single ratio R=1/3.  It was
    # counted four times and must be counted once.
    degree42_envelope = sum(degree42_contributions) - 3
    require(degree42_envelope == 35, "degree-42 torsion envelope changed")
    return degree34_envelope, degree42_envelope


def format_counter(counter: Counter) -> str:
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter, key=str))


def main() -> None:
    verify_tau_and_projectors()
    hidden_low = classify_hidden_low_degrees()
    (totals, a_zero, residual, distributions, orbit_data,
     associate_counts, nonzero_a) = enumerate_full_residual()
    boundary_resultant = verify_hidden_boundary_and_ratio()
    envelope34, envelope42 = verify_torsion_envelopes()

    print("W=0 tau-projector squeeze exact certificate")
    print("basis=[u,f,g,h] h=(v+omega^2*f+g)/2")
    print("tau=u:-omega; f:g; g:-omega*f; h:omega^2*h-omega*f")
    print("rosati_isometry=PASS tau_order=12 tau_trace_operator=0")
    print("projectors=P_u:m->a*u; P_v:m->d*v; P_hidden:m->2*ell PASS")
    print(f"full_theta degree34:{totals[34]} degree42:{totals[42]}")
    for degree in (6, 12, 24):
        count, orbits, determinant_norms = hidden_low[degree]
        print(
            f"hidden_degree{degree} vectors:{count} unit_tau_orbits:{orbits} "
            f"cyclic_det_norms:{format_counter(determinant_norms)}"
        )
    print("hidden_degree24=2*hidden_degree6 EXACT")
    print(f"hidden_two_torsion_boundary_resultant={boundary_resultant}")
    print("hidden_degree12_24_attachment_collapse=EXCLUDED_ON_GATE")
    print(f"nonzero_u_profiles degree34:{sorted(nonzero_a[34])}")
    print(f"nonzero_u_profiles degree42:{sorted(nonzero_a[42])}")
    print("nonzero_u_degree34_42_attachment_collapse=EXCLUDED")
    print(f"raw_a0 degree34:{len(a_zero[34])} degree42:{len(a_zero[42])}")
    print(
        f"spectral_residual degree34_vectors:{len(residual[34])} "
        f"orbits:{len(orbit_data[34][0])} orbit_sizes:{format_counter(orbit_data[34][1])}"
    )
    print(
        f"spectral_residual degree42_vectors:{len(residual[42])} "
        f"orbits:{len(orbit_data[42][0])} orbit_sizes:{format_counter(orbit_data[42][1])}"
    )
    print(f"degree34_Nd_K_distribution={dict(sorted(distributions[34].items()))}")
    print(f"degree42_Nd_K_distribution={dict(sorted(distributions[42].items()))}")
    print(f"degree34_d_ideal_classes={dict(sorted(associate_counts[34].items()))}")
    print(f"degree42_d_ideal_classes={dict(sorted(associate_counts[42].items()))}")
    print("remaining_torsion_condition=[d*(omega^2-1)]v(Q0)=O")
    print("marked_ratio_from_v_Xcube R=-1/(X^3+1)")
    print("degree42_Nd3_X0_ratio=-1 EXCLUDED_GATE")
    print("degree42_Nd3_ratio=1/3")
    print(f"torsion_ratio_envelopes degree34:{envelope34} degree42:{envelope42}")
    print("trace_zero_normalization: common_attachment_value_in_E0[12]")
    print("destroyed_data=mixed_attachment_cancellation; explicit ratio polynomials; walls; entry")
    print("verdict=S34_S42_OPEN residual_orbits=176+132")


if __name__ == "__main__":
    main()
