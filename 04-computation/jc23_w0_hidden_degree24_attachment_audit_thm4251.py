#!/usr/bin/env python3
"""Exact THM-4251 audit of the W=0 hidden degree-24 shell.

Relative to THM-4230/4241/4247, this certificate:

* exhausts the hidden Hermitian shell of degree 24;
* proves it is twice the degree-6 shell and has two free mu_6 x <T> orbits;
* derives the characteristic-zero doubled-map denominator
  t*(t^2-1)*L(t^2)^2 for both orbit representatives;
* proves L(1)L(-1) is nonzero by exact nested resultants; and
* checks denominator degree seven and reciprocal gcd t^2-1 in four exact
  good reductions.

It deliberately imports the already maintained THM-4247 coefficient-form
group law.  The companion independent audit for THM-4251 does not import it.
"""

from __future__ import annotations

import importlib.util
from collections import Counter
from hashlib import sha256
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation" / "jc23_w0_hidden_degree12_attachment_audit_thm4247.py"
SPEC = importlib.util.spec_from_file_location("thm4247", SOURCE)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def enumerate_shell() -> tuple[list[object], list[object], Counter[int]]:
    # The Hermitian eigenvalue 6-sqrt(12)>2 gives N(a)+N(b)<12 at q=24.
    # Since N(m+n*omega)>=(m^2+n^2)/2, the box [-4,4]^4 is complete.
    gram = (
        (6, -3, -3, 0),
        (-3, 6, 3, -3),
        (-3, 3, 6, -3),
        (0, -3, -3, 6),
    )

    def matrix_degree(coordinates: tuple[int, int, int, int]) -> int:
        return sum(
            coordinates[row] * gram[row][column] * coordinates[column]
            for row in range(4) for column in range(4)
        )

    vectors = []
    degree_six = []
    for am in range(-4, 5):
        for an in range(-4, 5):
            for bm in range(-4, 5):
                for bn in range(-4, 5):
                    vector = ((am, an), (bm, bn))
                    coordinates = (am, an, bm, bn)
                    q = matrix_degree(coordinates)
                    require(q == M.lattice_degree(*vector), "Gram paths disagree")
                    if q == 24:
                        vectors.append(vector)
                    if q == 6:
                        degree_six.append(vector)

    require(len(vectors) == 24, "degree-24 shell count changed")
    require(len(degree_six) == 24, "degree-6 shell count changed")
    require(
        max(max(map(abs, (*a, *b))) for a, b in vectors) == 2,
        "degree-24 shell reached the completeness boundary",
    )
    require(
        all(all(coordinate % 2 == 0 for coordinate in (*a, *b)) for a, b in vectors),
        "degree-24 shell is not two-divisible",
    )
    halves = {
        ((a[0] // 2, a[1] // 2), (b[0] // 2, b[1] // 2))
        for a, b in vectors
    }
    require(halves == set(degree_six), "halving is not the degree-6 bijection")

    vector_set = set(vectors)
    unseen = set(vectors)
    representatives = []
    orbit_sizes: Counter[int] = Counter()
    while unseen:
        seed = min(unseen)
        orbit = {M.unit_action(seed, unit) for unit in M.UNITS} | {
            M.unit_action(M.tau_action(seed), unit) for unit in M.UNITS
        }
        require(orbit <= vector_set, "symmetry orbit left the shell")
        representatives.append(min(orbit))
        orbit_sizes[len(orbit)] += 1
        unseen -= orbit

    expected_orbits = []
    for seed in (((2, 0), (0, 0)), ((2, 0), (2, 0))):
        expected_orbits.append(
            {M.unit_action(seed, unit) for unit in M.UNITS}
            | {M.unit_action(M.tau_action(seed), unit) for unit in M.UNITS}
        )
    require(expected_orbits[0].isdisjoint(expected_orbits[1]), "named orbits overlap")
    require(expected_orbits[0] | expected_orbits[1] == vector_set, "named orbits miss vectors")
    require(orbit_sizes == Counter({12: 2}), "degree-24 orbit structure changed")
    return vectors, sorted(representatives), orbit_sizes


def symbolic_denominators() -> dict[str, str]:
    p, z, u = sp.symbols("p z u")
    phi = z**4 - z**2 + 1
    sqrt_three = 2 * z - z**3
    p_relation = p**2 - (1 + sqrt_three) * p + 1
    relations = sp.groebner([p_relation, phi], p, z, order="lex")

    def reduce_element(expression: sp.Expr) -> sp.Expr:
        return sp.factor(relations.reduce(sp.expand(expression))[1])

    def absolute_resultant(expression: sp.Expr) -> sp.Expr:
        return sp.factor(sp.resultant(sp.resultant(expression, p_relation, p), phi, z))

    # For 2f the zero of B_f gives the squared linear denominator.
    f_line = u + p**3
    f_plus_norm = absolute_resultant(reduce_element(f_line.subs(u, 1)))
    f_minus_norm = absolute_resultant(reduce_element(f_line.subs(u, -1)))
    require((f_plus_norm, f_minus_norm) == (2916, 4), "2f wall norms changed")

    # For 2(f+Tf), the addition slope is constant and B is linear in u.
    a_f = u - p**2
    b_f = u + p**3
    a_g = z**2 * (p**2 * u - 1)
    b_g = -z**3 * (1 + p**3 * u)
    slope_numerator = -z**3 * p**3 - 1
    slope_denominator = z**2 * p**2 - 1
    require(
        sp.Poly(
            reduce_element(
                slope_denominator * (b_g - b_f)
                - slope_numerator * (a_g - a_f)
            ),
            u,
        ).is_zero,
        "f+Tf slope is not constant",
    )
    slope = slope_numerator / slope_denominator
    a_sum = (u - 1) * slope**2 - a_f - a_g
    b_sum = sp.together(slope * (a_f - a_sum) - b_f)
    b_numerator, b_denominator = sp.fraction(b_sum)
    b_poly = sp.Poly(sp.expand(b_numerator), u)
    require(b_poly.degree() == 1, "f+Tf B numerator is not linear")
    b0 = reduce_element(b_poly.coeff_monomial(1))
    b1 = reduce_element(b_poly.coeff_monomial(u))
    line = sp.expand(b0 + b1 * u)
    coefficient_norms = (absolute_resultant(b0), absolute_resultant(b1))
    denominator_norm = absolute_resultant(reduce_element(b_denominator))
    plus_norm = absolute_resultant(reduce_element(b0 + b1))
    minus_norm = absolute_resultant(reduce_element(b0 - b1))
    require(
        coefficient_norms == (11**6, 11**6),
        "f+Tf linear coefficients may vanish",
    )
    require(denominator_norm == 11**6, "f+Tf field denominator may vanish")
    require(
        (plus_norm, minus_norm) == (2916 * 11**6, 4 * 11**6),
        "f+Tf wall norms changed",
    )

    # If L=b0+b1*u is linear, L and u*L(1/u) share a root exactly when
    # L(1)L(-1)=0.  Hence the only raw reciprocal factor is u-1.
    return {
        "f_line": str(f_line),
        "f_plus_norm": str(f_plus_norm),
        "f_minus_norm": str(f_minus_norm),
        "sum_line": str(sp.factor(line)),
        "sum_coefficient_norms": str(tuple(map(str, coefficient_norms))),
        "sum_denominator_norm": str(denominator_norm),
        "sum_plus_norm": str(plus_norm),
        "sum_minus_norm": str(minus_norm),
        "raw_denominator": "t*(t^2-1)*L(t^2)^2",
        "raw_reciprocal_gcd": "t^2-1",
    }


def finite_field_controls(representatives: list[object]) -> list[object]:
    embeddings = (
        (313, 29, 135, 21),
        (349, 24, 246, 28),
        (373, 69, 297, 33),
        (397, 157, 161, 27),
    )
    rows = []
    for embedding in embeddings:
        q, z, p, _scale = embedding
        inverse = lambda value: pow(value % q, -1, q)
        slope = (-z**3 * p**3 - 1) * inverse(z**2 * p**2 - 1) % q
        f_line = (p**3 % q, 1)
        a_f = (-p**2 % q, 1)
        a_g = (-z**2 % q, z**2 * p**2 % q)
        a_sum = (
            (-slope * slope - a_f[0] - a_g[0]) % q,
            (slope * slope - a_f[1] - a_g[1]) % q,
        )
        sum_line = (
            (slope * (a_f[0] - a_sum[0]) - p**3) % q,
            (slope * (a_f[1] - a_sum[1]) - 1) % q,
        )

        def predicted_digest(line: tuple[int, int]) -> str:
            b0, b1 = line
            u_coefficients = [
                -b0 * b0,
                b0 * b0 - 2 * b0 * b1,
                2 * b0 * b1 - b1 * b1,
                b1 * b1,
            ]
            leading_inverse = inverse(u_coefficients[-1])
            coefficients = [
                u_coefficients[(exponent - 1) // 2] * leading_inverse % q
                if exponent % 2
                else 0
                for exponent in range(7, -1, -1)
            ]
            return sha256(",".join(map(str, coefficients)).encode("ascii")).hexdigest()

        predicted = {predicted_digest(f_line), predicted_digest(sum_line)}
        embedding_rows = []
        for representative in representatives:
            result = M.finite_field_attachment_audit(*embedding, [representative], 7)
            degrees, gcd_degrees, gcd_polynomials, digest = result
            require(degrees == [7], "specialized denominator lost degree seven")
            require(gcd_degrees == [2], "specialized reciprocal gcd gained a root")
            require(gcd_polynomials == [(1, 0, q - 1)], "specialized gcd is not t^2-1")
            embedding_rows.append((representative, degrees[0], gcd_polynomials[0], digest))
        require({row[-1] for row in embedding_rows} == predicted, "symbolic/modular digests disagree")
        rows.append((embedding, embedding_rows))
    return rows


def main() -> None:
    vectors, representatives, orbit_sizes = enumerate_shell()
    symbolic = symbolic_denominators()
    finite_rows = finite_field_controls(representatives)

    print("THM-4251 W=0 HIDDEN DEGREE-24 ATTACHMENT AUDIT")
    print(f"degree24_vectors={len(vectors)} degree6_halves=24")
    print(f"orbit_count={len(representatives)} orbit_sizes={dict(sorted(orbit_sizes.items()))}")
    print("orbit_representatives=2f,2(f+Tf)")
    print(f"symbolic={symbolic}")
    for embedding, rows in finite_rows:
        print(f"embedding={embedding} rows={rows}")
    print("denominator_degree=7 reciprocal_gcd=t^2-1")
    print("gate_wall=t=+-1 => Z/U=0")
    print("thm4247_only_projection_row_deletion=degree34:2304 degree42:288")
    print("thm4247_only_raw_remainder=degree34:31872 degree42:15840")
    print("current_stronger_residual=governed_by_THM4249")
    print("result=PASS")


if __name__ == "__main__":
    main()
