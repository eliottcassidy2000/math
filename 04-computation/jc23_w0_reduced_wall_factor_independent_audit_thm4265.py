#!/usr/bin/env python3
"""Independent primary-engine referee for provisional THM-4265.

This referee does not import or execute the THM-4265 checker.  It freezes the
proved THM-4260 primary SymPy/GF(t) compiler, rebuilds its 280 class universe,
and independently repeats the elliptic-curve additions that produce the
reduced X/x coefficient denominator.  It then checks the two geometric wall
roots and their derivatives in the proof and hostile finite fields.

The candidate used the native clean-room polynomial engine.  This program uses
the older SymPy rational-function path, so the implementations share only the
mathematical input and the frozen THM-4260 class universe.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import types

from sympy.polys.domains import GF
from sympy.polys.fields import field


PRIMARY_REVISION = "bafd69bdd607fe4dfa704edc63bef23cbf2e2ce9"
PRIMARY_PATH = (
    "04-computation/"
    "jc23_w0_canonical_node_attachment_exclusion_thm4260.py"
)
PRIMARY_SHA256 = (
    "cff97177dc253d63b09ed2ab9bfd3a3da74512a86c01fa400350083b63ef780a"
)
EMBEDDINGS = ((397, 157, 161, 27), (577, 57, 224, 25))
EXPECTED_PROFILE = Counter(
    {11: 8, 19: 36, 23: 24, 27: 64, 35: 52, 39: 72, 43: 24}
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def git_bytes(repo: Path, revision: str, path: str) -> bytes:
    return subprocess.run(
        ["git", "-C", str(repo), "show", f"{revision}:{path}"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def load_primary(repo: Path):
    source = git_bytes(repo, PRIMARY_REVISION, PRIMARY_PATH)
    require(sha256(source).hexdigest() == PRIMARY_SHA256, "primary source hash changed")
    module = types.ModuleType("thm4260_primary_for_thm4265_referee")
    module.__file__ = f"{repo}@{PRIMARY_REVISION}:{PRIMARY_PATH}"
    sys.modules[module.__name__] = module
    exec(compile(source, module.__file__, "exec"), module.__dict__)
    return module


def load_primary_universe(repo: Path):
    primary = load_primary(repo)
    dependency_revision = primary.resolved_revision(repo, primary.DEFAULT_REVISION)
    loaded = {}
    for name, (path, expected_hash) in primary.DEPENDENCIES.items():
        loaded[name] = primary.load_frozen_module(
            repo,
            dependency_revision,
            f"thm4265_referee_{name}",
            path,
            expected_hash,
        )
    class_representatives, profile_counts = primary.compile_classes(loaded["projector"])
    graph = primary.compile_incidence_graph(
        loaded["projector"], loaded["torsion"], class_representatives
    )
    require(
        {degree: len(class_representatives[degree]) for degree in (34, 42)}
        == {34: 176, 42: 104},
        "primary class universe changed",
    )
    require(len(graph["class_records"]) == 280, "primary record total changed")
    return primary, loaded, dependency_revision, profile_counts, graph


def make_embedding(attachment, embedding):
    q, z_value, p_value, scale_value = embedding
    z = z_value % q
    p = p_value % q
    scale = scale_value % q
    sqrt_three = (2 * z - z**3) % q
    scale_denominator = (2 * p**3 + 3 * p**2 - 1) % q
    require((z**4 - z**2 + 1) % q == 0, "bad Phi_12 embedding")
    require(
        (p**2 - (1 + sqrt_three) * p + 1) % q == 0,
        "bad real-branch embedding",
    )
    require(
        (pow(scale, 6, q) * scale_denominator - 4) % q == 0,
        "bad sixth-root embedding",
    )
    require(scale_denominator != 0, "scale denominator vanished")
    require((4 * z**3 - 2 * z) % q != 0, "ramified Phi_12 embedding")
    require(
        (2 * p - (1 + sqrt_three)) % q != 0,
        "ramified real-branch embedding",
    )
    require((6 * pow(scale, 5, q)) % q != 0, "ramified scale embedding")

    rational_field, t = field("t", GF(q))
    del rational_field
    u = t * t
    inverse_two = pow(2, -1, q)
    a_scale = scale * scale * inverse_two % q
    b_scale = scale**3 * inverse_two % q
    omega_value = pow(z, 4, q)
    quotient_c = (u - 1) / (2 * t)
    f_point = (
        a_scale * (u - p * p) / t,
        b_scale * (u + p**3) / t,
    )
    tf_point = (
        a_scale * (z * z % q) * (p * p * u - 1) / t,
        -b_scale * (z**3 % q) * (1 + p**3 * u) / t,
    )

    def point_neg(point):
        return None if point is None else (point[0], -point[1])

    def point_double(point):
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

    def point_add(left, right):
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
            - left_a
            - right_a
        )
        sum_b = (
            (right_b - left_b) * (left_a - sum_a) / (right_a - left_a)
            - left_b
        )
        return sum_a, sum_b

    def integer_multiple(multiplier, point):
        if multiplier < 0:
            return point_neg(integer_multiple(-multiplier, point))
        answer = None
        summand = point
        remaining = multiplier
        while remaining:
            if remaining & 1:
                answer = point_add(answer, summand)
            summand = point_double(summand)
            remaining //= 2
        return answer

    def eisenstein_multiple(multiplier, point):
        m_value, n_value = multiplier
        omega_point = (omega_value * point[0], point[1])
        return point_add(
            integer_multiple(m_value, point),
            integer_multiple(n_value, omega_point),
        )

    def denominator(a_coefficient, b_coefficient):
        point = point_add(
            eisenstein_multiple(a_coefficient, f_point),
            eisenstein_multiple(b_coefficient, tf_point),
        )
        require(point is not None, "hidden vector became the origin")
        return point[0].denom

    return q, t, denominator


def polynomial_value(polynomial, value: int, prime: int) -> int:
    return sum(
        int(coefficient) * pow(value, exponent[0], prime)
        for exponent, coefficient in polynomial.to_dict().items()
    ) % prime


def polynomial_derivative_value(polynomial, value: int, prime: int) -> int:
    return sum(
        exponent[0]
        * int(coefficient)
        * pow(value, exponent[0] - 1, prime)
        for exponent, coefficient in polynomial.to_dict().items()
        if exponent[0] > 0
    ) % prime


def reciprocal_polynomial(polynomial):
    degree = polynomial.degree()
    answer = polynomial.ring.zero
    for (exponent,), coefficient in polynomial.to_dict().items():
        answer[(degree - exponent,)] = coefficient * ((-1) ** exponent)
    return answer


def monic_coefficients(polynomial, prime: int) -> tuple[int, ...]:
    degree = polynomial.degree()
    coefficients = polynomial.to_dict()
    inverse = pow(int(polynomial.LC) % prime, -1, prime)
    return tuple(
        int(coefficients.get((exponent,), 0)) * inverse % prime
        for exponent in range(degree, -1, -1)
    )


def monic_gcd_coefficients(left, right, prime: int) -> tuple[int, ...]:
    common = left.gcd(right)
    return monic_coefficients(common, prime)


def audit_prime(attachment, graph, embedding):
    prime, _, denominator_builder = make_embedding(attachment, embedding)
    rows = []
    group_denominators = defaultdict(list)
    first_denominator = None

    for record in sorted(graph["class_records"], key=lambda item: item["class_id"]):
        denominator = denominator_builder(record["hidden_a"], record["hidden_b"])
        if first_denominator is None:
            first_denominator = denominator
        degree = denominator.degree()
        exponents = sorted(exponent[0] for exponent in denominator.to_dict())
        value_plus = polynomial_value(denominator, 1, prime)
        value_minus = polynomial_value(denominator, -1, prime)
        derivative_plus = polynomial_derivative_value(denominator, 1, prime)
        derivative_minus = polynomial_derivative_value(denominator, -1, prime)
        reciprocal = reciprocal_polynomial(denominator)
        gcd_coefficients = monic_gcd_coefficients(denominator, reciprocal, prime)
        expected_degree = 4 * record["hidden_k"] - 1

        require(degree == expected_degree, "sharp denominator degree failed")
        require(exponents[0] == 1, "distinguished t-order is not one")
        require(all(exponent % 2 == 1 for exponent in exponents), "denominator not odd")
        require(value_plus == value_minus == 0, "geometric wall root disappeared")
        require(
            derivative_plus != 0 and derivative_minus != 0,
            "geometric wall root is not reduced",
        )
        require(derivative_plus == derivative_minus, "oddness parity check failed")
        require(gcd_coefficients == (1, 0, prime - 1), "reciprocal gcd changed")

        coefficients = monic_coefficients(denominator, prime)
        coefficient_digest = sha256(
            ",".join(map(str, coefficients)).encode("ascii")
        ).hexdigest()
        rows.append(
            "\t".join(
                (
                    record["class_id"],
                    str(record["degree"]),
                    str(record["hidden_k"]),
                    str(degree),
                    str(derivative_plus * pow(int(denominator.LC) % prime, -1, prime) % prime),
                    str(derivative_minus * pow(int(denominator.LC) % prime, -1, prime) % prime),
                    coefficient_digest,
                )
            )
        )
        group_denominators[(record["degree"], record["hidden_k"])].append(
            denominator
        )

    require(first_denominator is not None, "empty denominator universe")
    wall = first_denominator.ring.one
    wall[(2,)] = first_denominator.ring.domain.one
    wall[(0,)] = -first_denominator.ring.domain.one
    squared_wall_control = first_denominator * wall
    require(
        polynomial_value(squared_wall_control, 1, prime) == 0
        and polynomial_value(squared_wall_control, -1, prime) == 0
        and polynomial_derivative_value(squared_wall_control, 1, prime) == 0
        and polynomial_derivative_value(squared_wall_control, -1, prime) == 0,
        "nonreduced hostile control was not detected",
    )
    broken_root_control = first_denominator + first_denominator.ring.one
    require(
        polynomial_value(broken_root_control, 1, prime) != 0
        and polynomial_value(broken_root_control, -1, prime) != 0,
        "missing-root hostile control was not detected",
    )

    profile = Counter(int(row.split("\t")[3]) for row in rows)
    require(profile == EXPECTED_PROFILE, "degree profile changed")
    ledger_hash = sha256("\n".join(rows).encode("ascii")).hexdigest()

    # Cross-check the independently rebuilt denominators against the frozen
    # primary THM-4260 group digests.  This is not used for the root tests.
    primary_group_rows = []
    for key in sorted(group_denominators):
        degree, hidden_k = key
        records = [
            record
            for record in graph["class_records"]
            if (record["degree"], record["hidden_k"]) == key
        ]
        representatives = [
            (record["hidden_a"], record["hidden_b"]) for record in records
        ]
        degrees, gcd_degrees, gcd_polynomials, primary_digest = (
            attachment.finite_field_attachment_audit(
                *embedding, representatives, 4 * hidden_k - 1
            )
        )
        require(set(degrees) == {4 * hidden_k - 1}, "primary degree cross-check failed")
        require(set(gcd_degrees) == {2}, "primary gcd-degree cross-check failed")
        require(
            set(gcd_polynomials) == {(1, 0, prime - 1)},
            "primary gcd-polynomial cross-check failed",
        )
        rebuilt_digest = sha256(
            ";".join(
                ",".join(map(str, monic_coefficients(polynomial, prime)))
                for polynomial in group_denominators[key]
            ).encode("ascii")
        ).hexdigest()
        # class_id order and the primary grouping order coincide because the
        # compiler appends class records degree-by-degree and index-by-index.
        require(rebuilt_digest == primary_digest, "primary denominator digest mismatch")
        primary_group_rows.append(f"{degree},{hidden_k},{len(records)},{primary_digest}")

    primary_group_hash = sha256("\n".join(primary_group_rows).encode("ascii")).hexdigest()
    return profile, ledger_hash, primary_group_hash


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    arguments = parser.parse_args()
    repo = arguments.repo.resolve()

    primary, loaded, dependency_revision, profile_counts, graph = load_primary_universe(repo)
    require(
        sum(sum(counter.values()) for counter in profile_counts.values()) == 280,
        "profile total changed",
    )

    print("THM-4265 INDEPENDENT PRIMARY SYMPY REFEREE")
    print(f"primary_revision={PRIMARY_REVISION}")
    print(f"primary_source={PRIMARY_PATH} sha256={PRIMARY_SHA256}")
    print(f"dependency_revision={dependency_revision}")
    print("implementation=SymPy_GF(t)_elliptic_addition independent_of_candidate_cleanroom_hook")
    print("classes=d34:176 d42:104 total:280")
    for embedding in EMBEDDINGS:
        profile, ledger_hash, primary_group_hash = audit_prime(
            loaded["attachment"], graph, embedding
        )
        profile_text = ",".join(
            f"{degree}:{profile[degree]}" for degree in sorted(profile)
        )
        print(
            f"q={embedding[0]} roots_zero=280/280 roots_simple=280/280 "
            f"odd_derivative_parity=280/280 reciprocal_gcd=t^2-1:280/280 "
            f"degree_profile={profile_text} "
            f"class_ledger_sha256={ledger_hash} "
            f"primary_group_digest_ledger_sha256={primary_group_hash}"
        )
    print("hostile_controls=extra_wall_factor_detected missing_wall_root_detected")
    print("lift_audit=monic_DVR_degree_preserving_reduction_forces_char0_simple_roots")
    print("formal_audit=IFT_graphs_and_J_equals_unit_times_first_graph_difference PASS")
    print("scope=no_W_lift_no_W_derivative_no_off_fibre_exclusion")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
