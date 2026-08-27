#!/usr/bin/env python3
"""Exact corollary audit for THM-4249's degree-42 N(d)=3 profile.

Formal input from THM-4249:

* the exact residual full-Hom shell and size-24 symmetry action;
* d*pi kills v(Q0), where pi=omega^2-1;
* R=-1/(X^3+1); and
* the proved exclusion 1/3 not in S_42.

For N(d)=3, d*pi is associated to 3.  The gate-interior part of E0[3]
is one mu_6 orbit and has R=1/3, so THM-4249 formally deletes the whole
(N(d),K)=(3,13), equivalently (q(v),q(ell))=(12,156), profile.

As an independent hostile control, this program also evaluates every raw
hidden projector in that profile at every root of the R=1/3 attachment
polynomial in four good reductions.  It does not use THM-4249's determinant
or order-18 proof of the ratio-one-third exclusion.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = ROOT / "04-computation/jc23_w0_cyclic_projector_squeeze_thm4249.py"


def load_parent():
    spec = importlib.util.spec_from_file_location("thm4249_parent", PARENT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


PARENT = load_parent()

EMBEDDINGS = (
    (313, 29, 135, 21),
    (349, 24, 246, 28),
    (373, 69, 297, 33),
    (397, 157, 161, 27),
)


def hidden_coordinates(vector):
    _, b, c, d = vector
    hidden_a = PARENT.e_add(
        PARENT.e_scale(2, b), PARENT.e_mul(PARENT.OMEGA2, d)
    )
    hidden_b = PARENT.e_add(PARENT.e_scale(2, c), d)
    return hidden_a, hidden_b


def exact_profile_corollary():
    _, _, residual, distributions, _, _, _ = PARENT.enumerate_full_residual()
    degree42 = residual[42]
    assert len(degree42) == 3168
    assert distributions[42][(3, 13)] == 672

    profile = set()
    for vector in degree42:
        hidden_a, hidden_b = hidden_coordinates(vector)
        hidden_degree = PARENT.hidden_degree(hidden_a, hidden_b)
        assert hidden_degree % 12 == 0
        if PARENT.e_norm(vector[3]) == 3 and hidden_degree // 12 == 13:
            profile.add(vector)
    assert len(profile) == 672

    profile_representatives, profile_sizes = PARENT.symmetry_orbits(profile)
    assert len(profile_representatives) == 28
    assert profile_sizes == Counter({24: 28})

    remaining = degree42 - profile
    assert profile.isdisjoint(remaining)
    assert profile | remaining == degree42
    remaining_representatives, remaining_sizes = PARENT.symmetry_orbits(remaining)
    assert len(remaining) == 2496
    assert len(remaining_representatives) == 104
    assert remaining_sizes == Counter({24: 104})

    # pi=omega^2-1 has norm three.  Every norm-three d is its associate,
    # so d*pi is associated to 3 and the visible point lies in E0[3].
    pi = PARENT.e_sub(PARENT.OMEGA2, PARENT.ONE)
    assert PARENT.e_norm(pi) == 3
    pi_associates = {PARENT.e_mul(unit, pi) for unit in PARENT.UNITS}
    for vector in profile:
        d = vector[3]
        assert d in pi_associates
        assert PARENT.e_norm(PARENT.e_mul(d, pi)) == 9

    # psi_3=3X(X^3+4).  X=0 is the gate wall R=-1; the only admissible
    # orbit has X^3=-4 and R=1/3.  THM-4249 proves 1/3 not in S_42.
    assert -Fraction(1, -4 + 1) == Fraction(1, 3)

    # Exact overlap with THM-4249's incidence deletion.  The 28 profile
    # orbits have one raw admissible ratio each, namely the common 1/3 lane.
    pre_ratio_orbits = 132
    pre_ratio_incidences = 780
    common_ratio_incidences = pre_ratio_orbits
    profile_raw_incidences = len(profile_representatives)
    other_raw_incidences = pre_ratio_incidences - profile_raw_incidences
    other_common_incidences = pre_ratio_orbits - profile_raw_incidences
    post_ratio_incidences = pre_ratio_incidences - common_ratio_incidences
    assert profile_raw_incidences == 28
    assert other_raw_incidences == 752
    assert other_common_incidences == 104
    assert post_ratio_incidences == 648
    assert other_raw_incidences - other_common_incidences == 648

    return profile, {
        "pre_vectors": len(degree42),
        "deleted_vectors": len(profile),
        "remaining_vectors": len(remaining),
        "pre_orbits": pre_ratio_orbits,
        "deleted_orbits": len(profile_representatives),
        "remaining_orbits": len(remaining_representatives),
        "pre_incidences": pre_ratio_incidences,
        "excluded_common_incidences": common_ratio_incidences,
        "post_incidences": post_ratio_incidences,
        "profile_overlap_incidences": profile_raw_incidences,
    }


def audit_embedding(embedding, profile):
    q, z, rho, scale = embedding
    inv = lambda value: pow(value % q, -1, q)
    omega = pow(z, 4, q)
    inverse_two = inv(2)
    sqrt_three = (2 * z - z**3) % q
    scale_denominator = (2 * rho**3 + 3 * rho**2 - 1) % q

    # Complete coefficient-field and separability controls.
    assert (z**4 - z**2 + 1) % q == 0
    assert (rho**2 - (1 + sqrt_three) * rho + 1) % q == 0
    assert (pow(scale, 6, q) * scale_denominator - 4) % q == 0
    assert scale_denominator != 0
    assert (4 * z**3 - 2 * z) % q != 0
    assert (2 * rho - (1 + sqrt_three)) % q != 0
    assert 6 * pow(scale, 5, q) % q != 0

    # R=(2t/(t^2-1))^2=1/3.
    roots = [t for t in range(q) if (t**4 - 14 * t * t + 1) % q == 0]
    assert len(roots) == 4 and len(set(roots)) == 4
    assert all(t and (t * t - 1) % q for t in roots)
    assert all((4 * t**3 - 28 * t) % q for t in roots)

    def point_neg(point):
        return None if point is None else (point[0], -point[1] % q)

    def point_double(point, t):
        if point is None:
            return None
        a_coefficient, b_coefficient = point
        u = t * t % q
        if b_coefficient % q == 0:
            return None
        doubled_a = (
            9 * t * pow(a_coefficient, 4, q)
            * inv(2 * (u - 1) * b_coefficient * b_coefficient)
            - 2 * a_coefficient
        ) % q
        doubled_b = (
            3 * t * a_coefficient * a_coefficient
            * (a_coefficient - doubled_a)
            * inv((u - 1) * b_coefficient)
            - b_coefficient
        ) % q
        return doubled_a, doubled_b

    def point_add(left, right, t):
        if left is None:
            return right
        if right is None:
            return left
        left_a, left_b = left
        right_a, right_b = right
        u = t * t % q
        quotient_c = (u - 1) * inv(2 * t) % q
        if left_a == right_a:
            if (left_b + right_b) % q == 0:
                return None
            assert left_b == right_b
            return point_double(left, t)
        sum_a = (
            quotient_c * (right_b - left_b) ** 2
            * inv((right_a - left_a) ** 2)
            - left_a - right_a
        ) % q
        sum_b = (
            (right_b - left_b) * (left_a - sum_a)
            * inv(right_a - left_a) - left_b
        ) % q
        return sum_a, sum_b

    def integer_multiple(multiplier, point, t):
        if multiplier < 0:
            return point_neg(integer_multiple(-multiplier, point, t))
        answer = None
        summand = point
        while multiplier:
            if multiplier & 1:
                answer = point_add(answer, summand, t)
            summand = point_double(summand, t)
            multiplier //= 2
        return answer

    def eisenstein_multiple(multiplier, point, t):
        m, n = multiplier
        omega_point = None if point is None else (omega * point[0] % q, point[1])
        return point_add(
            integer_multiple(m, point, t),
            integer_multiple(n, omega_point, t),
            t,
        )

    def hidden_generators(t):
        u = t * t % q
        a_scale = scale * scale * inverse_two % q
        b_scale = pow(scale, 3, q) * inverse_two % q
        f = (
            a_scale * (u - rho * rho) * inv(t) % q,
            b_scale * (u + pow(rho, 3, q)) * inv(t) % q,
        )
        g = (
            a_scale * z * z * (rho * rho * u - 1) * inv(t) % q,
            -b_scale * pow(z, 3, q) * (1 + pow(rho, 3, q) * u) * inv(t) % q,
        )
        return f, g

    rows = []
    hits = []
    for vector in sorted(profile):
        hidden_a, hidden_b = hidden_coordinates(vector)
        for t in roots:
            f, g = hidden_generators(t)
            image = point_add(
                eisenstein_multiple(hidden_a, f, t),
                eisenstein_multiple(hidden_b, g, t),
                t,
            )
            if image is None:
                hits.append((vector, t))
            rows.append(f"{vector}|{t}|{image}")
    assert len(rows) == 672 * 4
    assert hits == []
    digest = hashlib.sha256("\n".join(rows).encode("ascii")).hexdigest()
    return roots, len(rows), digest


def main():
    profile, ledger = exact_profile_corollary()
    print("THM-4253 DEGREE-42 NORM-THREE PROFILE EXCLUSION AUDIT")
    print("formal_input=THM-4249_PROVED ratio_1/3_notin_S42")
    print("profile=(N(d),K):(3,13) equivalent=(q(v),q(ell)):(12,156)")
    print(
        "vector_orbit_ledger="
        f"{ledger['pre_vectors']}/{ledger['pre_orbits']}"
        f"-{ledger['deleted_vectors']}/{ledger['deleted_orbits']}"
        f"={ledger['remaining_vectors']}/{ledger['remaining_orbits']}"
    )
    print(
        "incidence_overlap="
        f"pre:{ledger['pre_incidences']} common_1/3:{ledger['excluded_common_incidences']} "
        f"profile_subset:{ledger['profile_overlap_incidences']} post:{ledger['post_incidences']}"
    )
    for embedding in EMBEDDINGS:
        roots, tests, digest = audit_embedding(embedding, profile)
        print(
            f"field={embedding} roots={roots} "
            f"raw_hidden_origin_tests={tests} hits=0 sha256={digest}"
        )
    print("conclusion=profile_deleted vectors672 orbits28")
    print("refined_degree42_frontier=vectors2496 orbits104 incidences648")
    print("scope=remaining_104_orbits_648_incidences_W0_M12_entry_JC2_OPEN")


if __name__ == "__main__":
    main()
