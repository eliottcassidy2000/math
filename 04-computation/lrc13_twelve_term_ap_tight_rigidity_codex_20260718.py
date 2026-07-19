#!/usr/bin/env python3
"""Exact replay for THM-1171, twelve-term AP tight rigidity.

Only the Python standard library is used.  Every phase calculation is exact.
The scan is a replay of the proof's arithmetic consumers, not a finite-box
substitute for the symbolic argument in the theorem file.
"""

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    """An optimization-stable assertion (unlike the assert statement)."""
    if not condition:
        raise AssertionError(message)


def circle_distance(x: Fraction) -> Fraction:
    residue = x - (x.numerator // x.denominator)
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(circle_distance(Fraction(speed) * t) for speed in speeds)


def normalized_parameters(a: int, d: int) -> tuple[int, int, int]:
    g = gcd(a, d)
    return g, a // g, d // g


def primitive_step_witness(a: int, d: int) -> tuple[Fraction, Fraction, int, int]:
    """Return (time, clearance, residue, inverse clock) for primitive D >= 2."""
    g, A, D = normalized_parameters(a, d)
    require(D >= 2, "primitive-step witness requires D >= 2")
    r = D // 2
    j = (pow(A, -1, D) * r) % D
    t = Fraction(j, d)
    speeds = tuple(a + k * d for k in range(12))
    expected = Fraction(r, D)
    phases = tuple(circle_distance(Fraction(speed) * t) for speed in speeds)
    require(all(phase == expected for phase in phases), "common-phase collapse failed")
    require(clearance(speeds, t) == expected, "primitive clearance mismatch")
    require(expected >= Fraction(1, 3), "primitive clearance floor failed")
    require((A * j - r) % D == 0, "modular inverse clock failed")
    require(g * D == d, "normalization reconstruction failed")
    return t, expected, r, j


def consecutive_witness(a: int, d: int) -> tuple[Fraction, Fraction, int]:
    """Return the symmetric-endpoint witness on the normalized D = 1 branch."""
    g, A, D = normalized_parameters(a, d)
    require(D == 1, "consecutive witness requires normalized D = 1")
    denominator = g * (2 * A + 11)
    require(denominator == 2 * a + 11 * d, "clock denominator identity failed")
    t = Fraction(1, denominator)
    expected = Fraction(A, 2 * A + 11)
    speeds = tuple(a + k * d for k in range(12))
    require(clearance(speeds, t) == expected, "consecutive clearance mismatch")
    require(Fraction(speeds[0]) * t + Fraction(speeds[-1]) * t == 1,
            "endpoint symmetry failed")
    require((expected > Fraction(1, 13)) == (A > 1),
            "strict consecutive threshold equivalence failed")
    return t, expected, A


def transitive_tournament_fingerprint(n: int) -> tuple[tuple[int, ...], int, tuple[int, ...], int]:
    scores = tuple(range(n))
    directed_three_cycles = 0
    scc_sizes = (1,) * n
    hamiltonian_paths = 1
    return scores, directed_three_cycles, scc_sizes, hamiltonian_paths


def format_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def main() -> None:
    limit = 200
    rows = 0
    spread_rows = 0
    consecutive_rows = 0
    homogeneous_rows = 0
    strict_certificates = 0
    minimum_spread = Fraction(1, 1)

    for a in range(1, limit + 1):
        for d in range(1, limit + 1):
            rows += 1
            _, A, D = normalized_parameters(a, d)
            if D >= 2:
                spread_rows += 1
                _, witness, _, _ = primitive_step_witness(a, d)
                minimum_spread = min(minimum_spread, witness)
                require(witness > Fraction(1, 13), "spread AP was not strictly dispatched")
                strict_certificates += 1
            else:
                consecutive_rows += 1
                _, witness, normalized_start = consecutive_witness(a, d)
                require(normalized_start == A, "normalized start mismatch")
                if A == 1:
                    require(a == d, "homogeneous reconstruction failed")
                    require(witness == Fraction(1, 13), "homogeneous lower witness mismatch")
                    homogeneous_rows += 1
                else:
                    require(witness > Fraction(1, 13),
                            "shifted consecutive AP was not strictly dispatched")
                    strict_certificates += 1

    require(rows == limit * limit, "scan cardinality mismatch")
    require(homogeneous_rows == limit, "homogeneous row count mismatch")
    require(strict_certificates + homogeneous_rows == rows, "classification count mismatch")
    require(minimum_spread == Fraction(1, 3), "sharp scanned primitive-step floor mismatch")

    representatives = (
        (2, 41),
        (6, 15),
        (14, 7),
        (5, 5),
    )

    print("THM-1171 EXACT REPLAY")
    print(f"parameter_box=1..{limit}")
    print(f"rows={rows}")
    print(f"primitive_step_D_ge_2={spread_rows}")
    print(f"primitive_step_D_eq_1={consecutive_rows}")
    print(f"homogeneous_a_eq_d={homogeneous_rows}")
    print(f"strict_gt_1_over_13_certificates={strict_certificates}")
    print(f"minimum_spread_clearance={format_fraction(minimum_spread)}")
    print("classification_check=PASS")

    print("\nREPRESENTATIVE CERTIFICATES")
    for a, d in representatives:
        g, A, D = normalized_parameters(a, d)
        if D >= 2:
            t, witness, r, j = primitive_step_witness(a, d)
            print(
                f"a={a} d={d} g={g} A={A} D={D} branch=common_phase "
                f"r={r} j={j} t={t} clearance={witness}"
            )
        else:
            t, witness, _ = consecutive_witness(a, d)
            label = "homogeneous" if A == 1 else "symmetric_endpoints"
            print(
                f"a={a} d={d} g={g} A={A} D={D} branch={label} "
                f"t={t} clearance={witness}"
            )

    scores, cycles, scc_sizes, hpaths = transitive_tournament_fingerprint(12)
    print("\nTOURNAMENT AND CARRIER AUDIT")
    print("runner_pair_observable=sign(j-i) under positive-step gauge")
    print(f"runner_score_histogram={','.join(map(str, scores))}")
    print(f"runner_directed_3cycles={cycles}")
    print(f"runner_scc_sizes={','.join(map(str, scc_sizes))}")
    print(f"runner_hamiltonian_paths={hpaths}")
    print("tie_hamiltonian_path=0->1->...->11 (no ties)")
    print("runner_tournament_preserves=speed order only")
    print("runner_tournament_destroys=primitive step, modular inverse clock, common phase")
    print("faithful_spread_carrier=(A mod D,floor(D/2),inverse clock j/D)")
    print("faithful_consecutive_carrier=(endpoint owners k=0,11,symmetric clock)")
    print("preserved_LRC_predicate=explicit common-phase or endpoint witness")
    print("challenged_assumption=AP indices and runners need not be proof vertices")
    print("ordinary_and_pythonoptimized_outputs_must_match=YES")


if __name__ == "__main__":
    main()
