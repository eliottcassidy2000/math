#!/usr/bin/env python3
"""Scratch-only degree-24 hidden-shell probe built on THM-4247 machinery."""

from __future__ import annotations

import importlib.util
from collections import Counter
from hashlib import sha256
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "04-computation" / "jc23_w0_hidden_degree12_attachment_audit_thm4247.py"
SPEC = importlib.util.spec_from_file_location("thm4247", SOURCE)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def enumerate_shell():
    # The smallest Hermitian eigenvalue is 6-sqrt(12)>2, so q=24 gives
    # N(a)+N(b)<12.  Since N(m+n*omega)>=(m^2+n^2)/2, [-4,4]^4 is complete.
    gram = (
        (6, -3, -3, 0),
        (-3, 6, 3, -3),
        (-3, 3, 6, -3),
        (0, -3, -3, 6),
    )

    def matrix_degree(coordinates):
        return sum(
            coordinates[row]*gram[row][column]*coordinates[column]
            for row in range(4) for column in range(4)
        )

    vectors = []
    for am in range(-4, 5):
        for an in range(-4, 5):
            for bm in range(-4, 5):
                for bn in range(-4, 5):
                    vector = ((am, an), (bm, bn))
                    coordinates = (am, an, bm, bn)
                    assert matrix_degree(coordinates) == M.lattice_degree(*vector)
                    if matrix_degree(coordinates) == 24:
                        vectors.append(vector)
    vector_set = set(vectors)
    assert max(max(map(abs, (*a, *b))) for a, b in vectors) == 2
    unseen = set(vectors)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {
            M.unit_action(seed, unit) for unit in M.UNITS
        } | {
            M.unit_action(M.tau_action(seed), unit) for unit in M.UNITS
        }
        assert orbit <= vector_set
        orbits.append(orbit)
        unseen -= orbit
    seed_orbits = []
    for seed in (((2, 0), (0, 0)), ((2, 0), (2, 0))):
        seed_orbits.append(
            {M.unit_action(seed, unit) for unit in M.UNITS}
            | {M.unit_action(M.tau_action(seed), unit) for unit in M.UNITS}
        )
    assert seed_orbits[0].isdisjoint(seed_orbits[1])
    assert seed_orbits[0] | seed_orbits[1] == vector_set
    return vectors, sorted((min(orbit) for orbit in orbits)), Counter(map(len, orbits))


def symbolic_raw_denominators():
    """Certify that both doubled degree-six orbits have only wall reciprocity."""
    p, z, u, t = sp.symbols("p z u t")
    phi = z**4 - z**2 + 1
    sqrt_three = 2*z - z**3
    p_relation = p**2 - (1 + sqrt_three)*p + 1
    relations = sp.groebner([p_relation, phi], p, z, order="lex")

    def reduce_element(expression):
        return sp.factor(relations.reduce(sp.expand(expression))[1])

    def absolute_resultant(expression):
        return sp.factor(sp.resultant(
            sp.resultant(expression, p_relation, p), phi, z
        ))

    # Orbit f.  Before any cancellation, duplication has denominator
    # t*(u-1)*(u+p^3)^2.  Its reverse has the line 1+p^3*u.
    f_line = u + p**3
    f_at_plus = reduce_element(f_line.subs(u, 1))
    f_at_minus = reduce_element(f_line.subs(u, -1))
    f_plus_norm = absolute_resultant(f_at_plus)
    f_minus_norm = absolute_resultant(f_at_minus)
    assert f_plus_norm != 0 and f_minus_norm != 0

    # Orbit f+Tf.  The slope quotient is constant because the two linear
    # differences are proportional modulo the characteristic-zero relations.
    a_f = u - p**2
    b_f = u + p**3
    a_g = z**2*(p**2*u - 1)
    b_g = -z**3*(1 + p**3*u)
    slope_numerator = -z**3*p**3 - 1
    slope_denominator = z**2*p**2 - 1
    slope = slope_numerator / slope_denominator
    slope_check = sp.Poly(
        reduce_element(
            slope_denominator*(b_g-b_f)
            - slope_numerator*(a_g-a_f)
        ), u
    )
    assert slope_check.is_zero
    a_sum = (u-1)*slope**2 - a_f - a_g
    b_sum = sp.together(slope*(a_f-a_sum) - b_f)
    b_numerator, b_denominator = sp.fraction(b_sum)
    b_poly = sp.Poly(sp.expand(b_numerator), u)
    assert b_poly.degree() == 1
    b0 = reduce_element(b_poly.coeff_monomial(1))
    b1 = reduce_element(b_poly.coeff_monomial(u))
    line = sp.expand(b0+b1*u)
    line_plus = reduce_element(b0+b1)
    line_minus = reduce_element(b0-b1)
    denominator_scalar = reduce_element(b_denominator)
    denominator_norm = absolute_resultant(denominator_scalar)
    plus_norm = absolute_resultant(line_plus)
    minus_norm = absolute_resultant(line_minus)
    assert denominator_norm != 0 and plus_norm != 0 and minus_norm != 0

    # For D=(u-1)L(u)^2, the reverse is -(u-1)(u*L(1/u))^2.
    # A second common root exists iff the root of L is its own reciprocal,
    # equivalently L(1)L(-1)=0.  The nonzero norms above rule this out.
    return {
        "f_line": str(f_line),
        "f_plus_norm": str(f_plus_norm),
        "f_minus_norm": str(f_minus_norm),
        "sum_line": str(sp.factor(line)),
        "sum_line_plus": str(line_plus),
        "sum_line_minus": str(line_minus),
        "sum_denominator_scalar": str(denominator_scalar),
        "sum_denominator_norm": str(denominator_norm),
        "sum_plus_norm": str(plus_norm),
        "sum_minus_norm": str(minus_norm),
        "raw_reciprocal_gcd": "u-1 (equivalently t^2-1)",
    }


def denominator_probe(representatives):
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
        slope = (-z**3*p**3-1)*inverse(z**2*p**2-1) % q
        # Coefficient pairs are constant, u for the lines B_f and B_(f+Tf).
        f_line = (p**3 % q, 1)
        a_f = (-p**2 % q, 1)
        a_g = (-z**2 % q, z**2*p**2 % q)
        a_sum = (
            (-slope*slope-a_f[0]-a_g[0]) % q,
            (slope*slope-a_f[1]-a_g[1]) % q,
        )
        sum_line = (
            (slope*(a_f[0]-a_sum[0])-p**3) % q,
            (slope*(a_f[1]-a_sum[1])-1) % q,
        )

        def predicted_digest(line):
            # Multiply t*(u-1)*(line[0]+line[1]*u)^2 and normalize in t.
            b0, b1 = line
            u_coefficients = [
                -b0*b0,
                b0*b0-2*b0*b1,
                2*b0*b1-b1*b1,
                b1*b1,
            ]
            leading_inverse = inverse(u_coefficients[-1])
            # Descending t coefficients, with zero even slots.
            coefficients = []
            for exponent in range(7, -1, -1):
                coefficients.append(
                    u_coefficients[(exponent-1)//2]*leading_inverse % q
                    if exponent % 2 else 0
                )
            return sha256(",".join(map(str, coefficients)).encode("ascii")).hexdigest()

        predicted_digests = {predicted_digest(f_line), predicted_digest(sum_line)}
        embedding_rows = []
        for representative in representatives:
            result = None
            for degree in (1, 3, 5, 7):
                try:
                    result = M.finite_field_attachment_audit(
                        *embedding, [representative], degree
                    )
                except AssertionError as error:
                    if str(error) != "moving denominator degree changed":
                        raise
                else:
                    break
            assert result is not None
            den_degrees, gcd_degrees, gcd_polynomials, digest = result
            embedding_rows.append(
                (representative, den_degrees[0], gcd_degrees[0], gcd_polynomials[0], digest)
            )
        assert {row[-1] for row in embedding_rows} == predicted_digests
        rows.append((embedding, embedding_rows))
    return rows


def main():
    vectors, representatives, orbit_sizes = enumerate_shell()
    degree_six = {
        ((am, an), (bm, bn))
        for am in range(-4, 5)
        for an in range(-4, 5)
        for bm in range(-4, 5)
        for bn in range(-4, 5)
        if M.lattice_degree((am, an), (bm, bn)) == 6
    }
    assert all(all(coordinate % 2 == 0 for coordinate in (*a, *b)) for a, b in vectors)
    halves = {
        ((a[0]//2, a[1]//2), (b[0]//2, b[1]//2)) for a, b in vectors
    }
    assert halves == degree_six
    print("degree24_vectors", len(vectors))
    print("degree6_vectors", len(degree_six))
    print("degree24_equals_twice_degree6", halves == degree_six)
    print("orbit_count", len(representatives))
    print("orbit_sizes", sorted(orbit_sizes.items()))
    print("representatives")
    for representative in representatives:
        print(representative, "det_norm", M.e_norm(M.determinant_element(representative)))
    print("symbolic", symbolic_raw_denominators())
    for embedding, rows in denominator_probe(representatives):
        print("embedding", embedding)
        for row in rows:
            print(row)


if __name__ == "__main__":
    main()
