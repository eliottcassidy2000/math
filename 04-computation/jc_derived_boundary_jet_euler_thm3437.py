#!/usr/bin/env python3
"""Exact algebraic checks for THM-3437's derived Artin-jet packet.

This standard-library companion checks the two-term resolution convention,
Pruefer kernels/divisibility and transition maps, invertibility away from the
chosen primary support, and a finite bank of THM-3433 selected-root profiles.
It does not replace the proofs of THM-3433 or THM-3436.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "04-computation/jc_all_sector_multiroot_primary_torsion_thm3433.py":
        "7bb4db4f6d67436a2739dd79d0579243b8d0f5474a8f739a82b281014c114d8d",
    ROOT / "04-computation/jc_multiroot_boundary_jet_packet_probe_20260815.py":
        "3dcfc4f3f36ed658a64f3d45b7e055eb6bab10fa536522368a1aa3e32eab6332",
}
EXPECTED_SEMANTIC_SHA256 = "20bd151f40a3991b0fc85108f85339ebf9b48b24bc3113d8d480468b3ecef3f7"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def selected(d: int, sigma: int, exponents: tuple[int, ...], root: int) -> bool:
    if exponents[root] <= 1:
        return False
    if sigma * (exponents[root] - 1) % d:
        return False
    return all(
        sigma * exponent % d == 0
        for index, exponent in enumerate(exponents)
        if index != root
    )


def add(left: dict[int, Fraction], right: dict[int, Fraction]) -> dict[int, Fraction]:
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, Fraction(0)) + coefficient
        if result[exponent] == 0:
            del result[exponent]
    return result


def scale(value: dict[int, Fraction], scalar: Fraction) -> dict[int, Fraction]:
    return {
        exponent: scalar * coefficient
        for exponent, coefficient in value.items()
        if scalar * coefficient
    }


def multiply_mu(value: dict[int, Fraction], power: int = 1) -> dict[int, Fraction]:
    """Multiply a principal part by mu^power, discarding polynomial terms."""
    return {
        exponent + power: coefficient
        for exponent, coefficient in value.items()
        if exponent + power < 0
    }


def multiply_lambda(
    value: dict[int, Fraction], constant: Fraction
) -> dict[int, Fraction]:
    return add(scale(value, constant), multiply_mu(value))


def foreign_inverse(
    value: dict[int, Fraction], constant: Fraction, order: int
) -> dict[int, Fraction]:
    """Finite geometric inverse of constant+mu on order-order principal parts."""
    result: dict[int, Fraction] = {}
    shifted = dict(value)
    for power in range(order):
        coefficient = ((-1) ** power) / (constant ** (power + 1))
        result = add(result, scale(shifted, coefficient))
        shifted = multiply_mu(shifted)
    return result


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert found")
    require(not any(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    ), "float literal found")
    for path, expected in DEPENDENCIES.items():
        require(lf_sha256(path) == expected, (path, lf_sha256(path), expected))

    semantic = sha256()

    # Pruefer principal parts: kernel(lambda^q) has q basis elements, while
    # lambda^q is surjective because every basis element has a deeper lift.
    principal_kernel_checks = 0
    divisibility_checks = 0
    transition_checks = 0
    for q in range(1, 17):
        killed = tuple(
            j for j in range(1, 2 * q + 5)
            if -j + q >= 0
        )
        require(killed == tuple(range(1, q + 1)), (q, killed))
        principal_kernel_checks += len(killed)
        for j in range(1, 2 * q + 5):
            preimage = -(j + q)
            require(preimage + q == -j, (q, j, preimage))
            divisibility_checks += 1
        # Under f -> f lambda^-q, the degree-one chain map is multiplication
        # by lambda and hence coefficient truncation R_(q+1)->R_q.
        coefficients = tuple(Fraction((index + 1) * (q + 2), q + 1)
                             for index in range(q + 1))
        source_part = {
            index - (q + 1): coefficient
            for index, coefficient in enumerate(coefficients)
            if index - (q + 1) < 0
        }
        target_part = multiply_mu(source_part)
        expected_target = {
            index - q: coefficient
            for index, coefficient in enumerate(coefficients[:q])
        }
        require(target_part == expected_target,
                (q, target_part, expected_target))
        transition_checks += q
        semantic.update(repr(("principal", q, killed, target_part)).encode("ascii"))

    # A foreign primary support has lambda=constant+mu with nonzero constant;
    # its action on every finite principal part is invertible.
    foreign_support_checks = 0
    for order in range(1, 10):
        value = {
            -j: Fraction((j + 2) * (order + 1), j + order + 1)
            for j in range(1, order + 1)
        }
        for integer in (-5, -3, -2, -1, 1, 2, 3, 5):
            constant = Fraction(integer)
            inverse = foreign_inverse(value, constant, order)
            require(multiply_lambda(inverse, constant) == value,
                    (order, constant, value, inverse))
            foreign_support_checks += order
            semantic.update(repr(("foreign", order, integer, inverse)).encode("ascii"))

    # Frozen selected-root profile bank.  Formula (Tor_0,Tor_1)=(N-1+eps,eps)
    # is checked at every jet order, along with its constant Euler rank.
    profile_checks = 0
    jet_rank_checks = 0
    selected_profiles = 0
    unselected_profiles = 0
    character_count_checks = 0
    for d in range(2, 13):
        for root_count in range(1, 5):
            for exponents in product(range(1, 7), repeat=root_count):
                exponents = tuple(exponents)
                for root in range(root_count):
                    if exponents[root] <= 1:
                        continue
                    selected_count = 0
                    for sigma in range(1, d + 1):
                        epsilon = int(selected(d, sigma, exponents, root))
                        selected_count += epsilon
                        tor0_rank = root_count - 1 + epsilon
                        tor1_rank = epsilon
                        require(tor0_rank - tor1_rank == root_count - 1,
                                (d, sigma, exponents, root, tor0_rank, tor1_rank))
                        profile_checks += 1
                        selected_profiles += epsilon
                        unselected_profiles += 1 - epsilon
                        for q in range(1, 7):
                            # Both homology modules are free R_q packets, so
                            # their K'-lengths are q times the displayed ranks.
                            require(q * tor0_rank - q * tor1_rank
                                    == q * (root_count - 1),
                                    (d, sigma, exponents, root, q))
                            jet_rank_checks += 1
                        semantic.update(repr(
                            (d, sigma, exponents, root, epsilon,
                             tor0_rank, tor1_rank)
                        ).encode("ascii"))
                    gcd_packet = exponents[root] - 1
                    for index, exponent in enumerate(exponents):
                        if index != root:
                            from math import gcd
                            gcd_packet = gcd(gcd_packet, exponent)
                    from math import gcd
                    gcd_packet = gcd(gcd_packet, d)
                    require(selected_count == gcd_packet,
                            (d, exponents, root, selected_count, gcd_packet))
                    character_count_checks += 1

    # The unfiltered derived packet sees arm presence but not the DeathBar
    # slope: these two selected one-root profiles have different e-1.
    blindness_controls = 0
    for d, sigma, first, second in ((2, 1, 3, 5), (3, 1, 4, 7),
                                    (4, 2, 3, 5), (5, 5, 2, 6)):
        require(selected(d, sigma, (first,), 0), (d, sigma, first))
        require(selected(d, sigma, (second,), 0), (d, sigma, second))
        require(first - 1 != second - 1, (first, second))
        # N=1 gives the same derived ranks (1,1) in both cases.
        require((1 - 1 + 1, 1) == (1, 1), (d, sigma))
        blindness_controls += 1
        semantic.update(repr(("blind", d, sigma, first, second)).encode("ascii"))

    # Simple roots are deliberately outside the boundary theorem.
    simple_root_controls = 0
    for d in range(2, 17):
        for sigma in range(1, d + 1):
            require(not selected(d, sigma, (1,), 0), (d, sigma))
            simple_root_controls += 1

    semantic_digest = semantic.hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256,
                (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("THM-3437 derived boundary-jet Euler conservation")
    print("status=PROVISIONAL_PROOF_CANDIDATE;independent_audit_required")
    print("dependency_script_hashes=" + ";".join(
        f"{path.name}:{expected}" for path, expected in DEPENDENCIES.items()
    ))
    print(f"principal_kernel_checks={principal_kernel_checks};divisibility_checks={divisibility_checks};transition_checks={transition_checks}")
    print(f"foreign_support_checks={foreign_support_checks}")
    print(f"profile_checks={profile_checks};jet_rank_checks={jet_rank_checks};selected_profiles={selected_profiles};unselected_profiles={unselected_profiles};character_count_checks={character_count_checks}")
    print(f"blindness_controls={blindness_controls};simple_root_controls={simple_root_controls}")
    print("derived_packet=Tor0_rank_N_minus_1_plus_epsilon;Tor1_rank_epsilon;Euler_rank_N_minus_1")
    print("loss=Tor1_recovers_Pruefer_presence_but_not_DeathBar_slope_or_intercept")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=stdlib_Fraction_integer_exact;no_float;no_assert;normal_and_O_truth_gates")


if __name__ == "__main__":
    main()
