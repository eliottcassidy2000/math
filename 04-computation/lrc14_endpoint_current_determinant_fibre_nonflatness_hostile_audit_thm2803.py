#!/usr/bin/env python3
"""Independent affine-line audit for THM-2803.

This companion reuses only the proved THM-2779 reconstruction of the
canonical THM-2625 endpoint factors.  It does not import the THM-2803
candidate, choose a transverse vector, or use its sequence/minor routines.

For each nonzero direction s it enumerates every affine determinant line
directly.  A possible equivariant identification of two lines is determined
by the image of one lexicographically chosen origin.  All thirteen images
are tried, and scalar proportionality is tested by anchored cross-products.
The reversed orientation is audited separately.

Script:
  04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_hostile_audit_thm2803.py
Output:
  05-knowledge/results/lrc14_endpoint_current_determinant_fibre_nonflatness_hostile_audit_thm2803.out
"""

import hashlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_bockstein_symplectic_endpoint_gate_thm2779.py"
)
SPEC = importlib.util.spec_from_file_location(
    "thm2779_endpoint_gate_2803_independent",
    SOURCE,
)
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)

P = 13
POINTS = GATE.POINTS
NONZERO = GATE.NONZERO
EXPECTED_PRIMES = (352341050142921841, 956354278959359281)
EXPECTED_DIGESTS = (
    (
        "d5b84c7702cead0b4fa13842d6207a0ba3acc7d68c9d0ac334518df9e6cc5cfe",
        "b3b55a94876cc5da043448e33bf5f1253e58e0950255df33b9d70253c2851d3f",
        "7efb1cb22fea8a123340e220de26b9595f92c843f2347d472500e5ad0d89bac3",
        "b3b55a94876cc5da043448e33bf5f1253e58e0950255df33b9d70253c2851d3f",
        "47323b9ce436dfaebaa6e2eda6ca7fe664a13229595778239da7d088e74524bc",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
    ),
    (
        "2b0fe21a1f03eccf176398892d47924c764a9dc6850e4e48704452bc4517f31b",
        "b3b55a94876cc5da043448e33bf5f1253e58e0950255df33b9d70253c2851d3f",
        "cfcde42a4e8b81420caa92cb8bf7f68a47a602bbcc26ba2868a6d7e9ee630e8f",
        "b3b55a94876cc5da043448e33bf5f1253e58e0950255df33b9d70253c2851d3f",
        "a25975329d4f8825069efd801754d3bd51cf8ec453b7cdc0b9bb3e7bfa615383",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
        "1a6a41f6e9fea7a2d86551e7369572dbc64c15cba20e392b07818d42299e29d5",
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(values):
    payload = b"".join(int(value).to_bytes(8, "big") for value in values)
    return hashlib.sha256(payload).hexdigest()


def add(point, step):
    return ((point[0] + step[0]) % P, (point[1] + step[1]) % P)


def negate(step):
    return (-step[0] % P, -step[1] % P)


def advance(point, step, times):
    result = point
    for _ in range(times):
        result = add(result, step)
    return result


def endpoint_value(field, step, origin):
    prime, left, right = field
    return left[add(origin, step)] * right[origin] % prime


def determinant_lines(step):
    lines = []
    for delta in range(P):
        line = tuple(
            point for point in POINTS
            if GATE.det(step, point) == delta
        )
        require(len(line) == P, "affine determinant-line size drift")
        lines.append(line)
    require(
        set().union(*(set(line) for line in lines)) == set(POINTS),
        "affine determinant lines do not cover the plane",
    )
    require(
        sum(len(line) for line in lines) == len(POINTS),
        "affine determinant lines overlap",
    )
    return tuple(lines)


def anchored_comparison(field, step, source_origin, target_origin, reverse_target):
    """Test one affine line map directly in endpoint coordinates."""
    prime = field[0]
    target_step = negate(step) if reverse_target else step
    source_values = tuple(
        endpoint_value(field, step, advance(source_origin, step, index))
        for index in range(P)
    )
    target_values = tuple(
        endpoint_value(field, step, advance(target_origin, target_step, index))
        for index in range(P)
    )
    require(all(source_values) and all(target_values), "line-map current vanished")
    failures = 0
    first_witness = 0
    for index in range(1, P):
        difference = (
            target_values[index] * source_values[0]
            - target_values[0] * source_values[index]
        ) % prime
        if difference:
            failures += 1
            if first_witness == 0:
                first_witness = difference
    pairwise_nonzero = 0
    for left_index in range(P):
        for right_index in range(left_index + 1, P):
            cross = (
                target_values[right_index] * source_values[left_index]
                - target_values[left_index] * source_values[right_index]
            ) % prime
            pairwise_nonzero += cross != 0
    return failures, first_witness, pairwise_nonzero


def verify_origin_controls(field, step, line):
    """The same line must be found projectively after every true origin move."""
    origin = min(line)
    positive = 0
    for target_origin in line:
        failures, witness, pairwise_nonzero = anchored_comparison(
            field,
            step,
            origin,
            target_origin,
            False,
        )
        if failures == 0:
            require(witness == 0, "zero-failure control acquired a witness")
            require(pairwise_nonzero == 0, "same-line scalar ratio control drift")
            positive += 1
    require(positive >= 1, "same-line origin positive control failed")

    reverse_positive = 0
    for target_origin in line:
        failures, witness, pairwise_nonzero = anchored_comparison(
            field,
            step,
            origin,
            target_origin,
            True,
        )
        if failures == 0:
            require(witness == 0, "reversal control acquired a witness")
            require(pairwise_nonzero == 0, "same-line reversal ratio control drift")
            reverse_positive += 1
    return positive, reverse_positive


def verify_boolean_marked_power_boundary(prime):
    """Independent exhaustive controls for the finite Mathieu analogue."""
    exponents = (1, 2, 13)
    coarse_zero = 0
    marked_formula = 0
    marked_nonzero = 0
    coarse_marker_zero = 0
    for mask in range(1 << P):
        low = tuple((mask >> residue) & 1 for residue in range(P))
        for exponent in exponents:
            for residue in range(P):
                value = pow(low[residue], exponent, prime)
                unmarked = (value - value) % prime
                high_marked = (0 * value - 1 * value) % prime
                expected = -value % prime
                low_marked = (low[residue] * value - low[residue] * value) % prime
                require(unmarked == 0, "Boolean coarse pullback did not cancel")
                require(
                    high_marked == expected,
                    "Boolean high-digit marked-power formula drift",
                )
                require(
                    low_marked == 0,
                    "coarse marker incorrectly detected a forgotten digit",
                )
                coarse_zero += 1
                marked_formula += 1
                marked_nonzero += high_marked != 0
                coarse_marker_zero += 1
    expected_total = (1 << P) * len(exponents) * P
    expected_nonzero = (1 << (P - 1)) * len(exponents) * P
    require(coarse_zero == expected_total, "Boolean coarse census drift")
    require(marked_formula == expected_total, "Boolean marked census drift")
    require(marked_nonzero == expected_nonzero, "Boolean support census drift")
    require(coarse_marker_zero == expected_total, "coarse-marker hostile drift")
    return expected_total, expected_nonzero


def analyze_field(field, field_index):
    prime = field[0]
    require(prime == EXPECTED_PRIMES[field_index], "certified field drift")
    pair_count = 0
    oriented_maps = 0
    reversed_maps = 0
    oriented_matches = 0
    reversed_matches = 0
    oriented_failure_bank = []
    oriented_witness_bank = []
    reversed_failure_bank = []
    reversed_witness_bank = []
    oriented_pairwise_cross_bank = []
    reversed_pairwise_cross_bank = []
    line_value_bank = []
    same_line_origin_controls = 0
    same_line_reversal_controls = 0

    for step in NONZERO:
        lines = determinant_lines(step)
        positive, reverse_positive = verify_origin_controls(
            field,
            step,
            lines[0],
        )
        same_line_origin_controls += positive
        same_line_reversal_controls += reverse_positive

        for line in lines:
            origin = min(line)
            orbit = tuple(advance(origin, step, index) for index in range(P))
            require(set(orbit) == set(line), "central orbit misses affine line")
            line_value_bank.extend(
                endpoint_value(field, step, point)
                for point in orbit
            )

        for left_delta in range(P):
            source_origin = min(lines[left_delta])
            for right_delta in range(left_delta + 1, P):
                pair_count += 1
                for target_origin in lines[right_delta]:
                    failures, witness, pairwise_nonzero = anchored_comparison(
                        field,
                        step,
                        source_origin,
                        target_origin,
                        False,
                    )
                    oriented_maps += 1
                    oriented_matches += failures == 0
                    oriented_failure_bank.append(failures)
                    oriented_witness_bank.append(witness)
                    oriented_pairwise_cross_bank.append(pairwise_nonzero)

                    reverse_failures, reverse_witness, reverse_pairwise = (
                        anchored_comparison(
                        field,
                        step,
                        source_origin,
                        target_origin,
                        True,
                        )
                    )
                    reversed_maps += 1
                    reversed_matches += reverse_failures == 0
                    reversed_failure_bank.append(reverse_failures)
                    reversed_witness_bank.append(reverse_witness)
                    reversed_pairwise_cross_bank.append(reverse_pairwise)

    expected_pairs = len(NONZERO) * P * (P - 1) // 2
    expected_maps = expected_pairs * P
    require(pair_count == expected_pairs, "independent pair universe drift")
    require(oriented_maps == expected_maps, "oriented affine-map census drift")
    require(reversed_maps == expected_maps, "reversed affine-map census drift")
    require(oriented_matches == 0, "oriented projective affine-line collision")
    require(reversed_matches == 0, "reversed projective affine-line collision")
    require(
        all(value == P - 1 for value in oriented_failure_bank),
        "oriented anchored-ratio collision",
    )
    require(all(oriented_witness_bank), "oriented first witness failure")
    require(
        all(value == P - 1 for value in reversed_failure_bank),
        "reversed anchored-ratio collision",
    )
    require(all(reversed_witness_bank), "reversed first witness failure")
    require(
        all(value == P * (P - 1) // 2 for value in oriented_pairwise_cross_bank),
        "oriented pointwise ratio collision",
    )
    require(
        all(value == P * (P - 1) // 2 for value in reversed_pairwise_cross_bank),
        "reversed pointwise ratio collision",
    )
    require(
        len(line_value_bank) == len(NONZERO) * P * P,
        "independent current profile universe drift",
    )
    require(all(line_value_bank), "independent current support hole")
    boolean_total, boolean_nonzero = verify_boolean_marked_power_boundary(prime)

    result = {
        "prime": prime,
        "pairs": pair_count,
        "oriented_maps": oriented_maps,
        "reversed_maps": reversed_maps,
        "oriented_matches": oriented_matches,
        "reversed_matches": reversed_matches,
        "oriented_failure_min": min(oriented_failure_bank),
        "oriented_failure_max": max(oriented_failure_bank),
        "reversed_failure_min": min(reversed_failure_bank),
        "reversed_failure_max": max(reversed_failure_bank),
        "oriented_cross_min": min(oriented_pairwise_cross_bank),
        "oriented_cross_max": max(oriented_pairwise_cross_bank),
        "reversed_cross_min": min(reversed_pairwise_cross_bank),
        "reversed_cross_max": max(reversed_pairwise_cross_bank),
        "same_line_origins": same_line_origin_controls,
        "same_line_reversals": same_line_reversal_controls,
        "boolean_total": boolean_total,
        "boolean_nonzero": boolean_nonzero,
        "line_value_digest": digest(line_value_bank),
        "oriented_failure_digest": digest(oriented_failure_bank),
        "oriented_witness_digest": digest(oriented_witness_bank),
        "reversed_failure_digest": digest(reversed_failure_bank),
        "reversed_witness_digest": digest(reversed_witness_bank),
        "oriented_cross_digest": digest(oriented_pairwise_cross_bank),
        "reversed_cross_digest": digest(reversed_pairwise_cross_bank),
    }
    require(
        (
            result["line_value_digest"],
            result["oriented_failure_digest"],
            result["oriented_witness_digest"],
            result["reversed_failure_digest"],
            result["reversed_witness_digest"],
            result["oriented_cross_digest"],
            result["reversed_cross_digest"],
        )
        == EXPECTED_DIGESTS[field_index],
        "independent deterministic field digest drift",
    )
    return result


def main():
    fields = GATE.build_endpoint_factors()
    require(
        tuple(field[0] for field in fields) == EXPECTED_PRIMES,
        "certified dual-field universe drift",
    )
    results = tuple(
        analyze_field(field, field_index)
        for field_index, field in enumerate(fields)
    )

    print("THM-2803 INDEPENDENT AFFINE-LINE HOSTILE AUDIT")
    print("status=VERIFIED-EXACT independent dual-field companion")
    print("universe=all_nonzero_s_all_determinant_line_pairs_all_equivariant_origins")
    print("method=direct_affine_lines_and_anchored_cross_products; no_transverse_frame")
    print("orientation_hostile=all_affine_antiequivariant_maps_also_tested")
    print("mathieu_scope=finite_pointwise_function_algebra; not_Weyl_SIC")
    for result in results:
        print(
            f"p={result['prime']} pairs={result['pairs']} "
            f"oriented_maps={result['oriented_maps']} "
            f"oriented_matches={result['oriented_matches']} "
            f"anchored_failures={result['oriented_failure_min']}.."
            f"{result['oriented_failure_max']} "
            f"nonzero_pairwise_crosses={result['oriented_cross_min']}.."
            f"{result['oriented_cross_max']}"
        )
        print(
            f"p={result['prime']} reversed_maps={result['reversed_maps']} "
            f"reversed_matches={result['reversed_matches']} "
            f"reversed_anchored_failures={result['reversed_failure_min']}.."
            f"{result['reversed_failure_max']} "
            f"reversed_nonzero_pairwise_crosses={result['reversed_cross_min']}.."
            f"{result['reversed_cross_max']} "
            f"same_line_origin_controls={result['same_line_origins']} "
            f"same_line_reversal_controls={result['same_line_reversals']}"
        )
        print(
            f"p={result['prime']} boolean_boundary_checks="
            f"{result['boolean_total']} "
            f"marked_nonzero={result['boolean_nonzero']} "
            f"coarse_marker_nonzero=0"
        )
        print(
            f"p={result['prime']} line_values_sha256="
            f"{result['line_value_digest']} "
            f"oriented_failure_sha256={result['oriented_failure_digest']} "
            f"oriented_witness_sha256={result['oriented_witness_digest']}"
        )
        print(
            f"p={result['prime']} reversed_failure_sha256="
            f"{result['reversed_failure_digest']} "
            f"reversed_witness_sha256={result['reversed_witness_digest']} "
            f"oriented_cross_count_sha256={result['oriented_cross_digest']} "
            f"reversed_cross_count_sha256={result['reversed_cross_digest']}"
        )
    print("scope=coefficient anti-quotient evidence only; no SIC implication, physical descent, positive floor, row exclusion, or LRC14")
    print("INDEPENDENT HOSTILE AUDIT PASSED")


if __name__ == "__main__":
    main()
