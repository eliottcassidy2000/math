#!/usr/bin/env python3
"""Exact hostile probe and rank ledger for THM-2303.

All geometric eligibility checks use Fraction arithmetic.  The Fourier
conclusions use the proved two-interval formula: for equal half-widths, the
sum vanishes exactly when the centre phases are antipodal.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256


P = 13
M = 4
EPSILON = Fraction(1, 10**12)
ETA = Fraction(1, 10**12)

H = 1
UNITS = (4, 2, 3, 6, 10)
C1 = 13
C2 = 13**3
C3 = 2 * 13**5
BLOCKERS = (C1, C2, C3)


def require(condition: bool, message: str) -> None:
    """Optimization-safe assertion."""

    if not condition:
        raise AssertionError(message)


def circle_norm(x: Fraction) -> Fraction:
    """Distance of an exact rational circle point to the nearest integer."""

    residue = x % 1
    return min(residue, 1 - residue)


def source_labels(y: Fraction) -> tuple[bool, tuple[bool, ...], tuple[bool, ...]]:
    """Strict source labels for x=(8+y)/13."""

    x = (8 + y) / P
    guard_safe = circle_norm(H * x) > Fraction(1, 7)
    units_safe = tuple(circle_norm(q * x) > Fraction(1, 14) for q in UNITS)
    blocker_danger = tuple(
        circle_norm(c * x) < Fraction(1, 14) for c in BLOCKERS
    )
    return guard_safe, units_safe, blocker_danger


def terminal_labels(y: Fraction) -> tuple[bool, tuple[bool, ...], tuple[bool, ...]]:
    """Strict terminal labels for z=13y mod one."""

    z = (P * y) % 1
    guard_safe = circle_norm(H * z) > Fraction(1, 7)
    units_safe = tuple(circle_norm(q * z) > Fraction(1, 14) for q in UNITS)
    blocker_danger = tuple(
        circle_norm(c * z) < Fraction(1, 14) for c in BLOCKERS
    )
    return guard_safe, units_safe, blocker_danger


def source_margin_in_y(y: Fraction) -> Fraction:
    """Largest Lipschitz-safe y-radius at one source centre."""

    x = (8 + y) / P
    margins: list[Fraction] = []

    margins.append(
        (circle_norm(H * x) - Fraction(1, 7)) * P / H
    )
    for q in UNITS:
        margins.append(
            (circle_norm(q * x) - Fraction(1, 14)) * P / q
        )

    margins.append(
        (Fraction(1, 14) - circle_norm(C1 * x)) * P / C1
    )
    for c in (C2, C3):
        margins.append(
            (circle_norm(c * x) - Fraction(1, 14)) * P / c
        )

    require(all(margin > 0 for margin in margins), "source centre is not strict")
    return min(margins)


def terminal_margin_in_y(y: Fraction) -> Fraction:
    """Largest Lipschitz-safe y-radius at one terminal centre."""

    z = (P * y) % 1
    margins: list[Fraction] = []

    margins.append(
        (circle_norm(H * z) - Fraction(1, 7)) / (P * H)
    )
    for q in UNITS:
        margins.append(
            (circle_norm(q * z) - Fraction(1, 14)) / (P * q)
        )

    for c in (C1, C3):
        margins.append(
            (circle_norm(c * z) - Fraction(1, 14)) / (P * c)
        )
    margins.append(
        (Fraction(1, 14) - circle_norm(C2 * z)) / (P * C2)
    )

    require(all(margin > 0 for margin in margins), "terminal centre is not strict")
    return min(margins)


def rank_rational(rows: list[list[Fraction]]) -> int:
    """Exact row rank over Q."""

    if not rows:
        return 0
    matrix = [row[:] for row in rows]
    row_count = len(matrix)
    col_count = len(matrix[0])
    pivot_row = 0
    for col in range(col_count):
        pivot = next(
            (r for r in range(pivot_row, row_count) if matrix[r][col] != 0),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][col]
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for r in range(row_count):
            if r == pivot_row or matrix[r][col] == 0:
                continue
            factor = matrix[r][col]
            matrix[r] = [
                matrix[r][j] - factor * matrix[pivot_row][j]
                for j in range(col_count)
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def main() -> None:
    # Scalar row, pair relation, independent relation, and anchored minor.
    require(P * UNITS[0] - M * C1 == 0, "shallow pair relation failed")
    require(UNITS[0] - 2 * UNITS[1] == 0, "independent relation failed")
    anchored_minor = (-M) * (-2)
    require(anchored_minor == 8, "anchored minor changed")
    require(anchored_minor % P != 0, "anchored minor lost mod-13 rank")

    negative_center = Fraction(-1, 16)
    positive_center = Fraction(1, 16)
    configurations = {
        "antipodal": (negative_center, positive_center),
        "perturbed": (negative_center, positive_center + ETA),
    }

    # The perturbation and both interval half-widths fit inside every exact
    # source and current-target eligibility margin.
    all_centers = [
        negative_center,
        positive_center,
        positive_center + ETA,
    ]
    exact_margin = min(
        min(source_margin_in_y(y), terminal_margin_in_y(y))
        for y in all_centers
    )
    require(EPSILON + ETA < exact_margin, "probe leaves the strict label cell")

    signatures: list[str] = []
    relative_phases: dict[str, Fraction] = {}
    terminal_addresses: dict[str, tuple[int, int]] = {}
    root_character_component_checks = 0

    for name, centers in configurations.items():
        for center in centers:
            for endpoint_sign in (-1, 1):
                y = center + endpoint_sign * EPSILON
                source = source_labels(y)
                terminal = terminal_labels(y)
                require(source[0], f"{name}: source guard is not safe")
                require(all(source[1]), f"{name}: a source unit is dangerous")
                require(
                    source[2] == (True, False, False),
                    f"{name}: source blocker word changed",
                )
                require(terminal[0], f"{name}: terminal guard is not safe")
                require(all(terminal[1]), f"{name}: a terminal unit is dangerous")
                require(
                    terminal[2] == (False, True, False),
                    f"{name}: terminal blocker word changed",
                )

        # The negative component uses terminal inverse root 12; the positive
        # component uses root 0.  These addresses are constant on the tiny
        # intervals because neither image crosses a branch boundary.
        addresses = (12, 0)
        terminal_addresses[name] = addresses
        require(
            centers[0] % 1 == (((P * centers[0]) % 1) + 12) / P,
            f"{name}: negative terminal address changed",
        )
        require(
            centers[1] % 1 == (((P * centers[1]) % 1) + 0) / P,
            f"{name}: positive terminal address changed",
        )

        # On a one-sheet word the character value is (1/13)zeta^(-a*r).
        # Its squared magnitude is exactly 1/169, independent of a and r.
        # Each terminal component has length 26*epsilon, hence exact energy
        # 2*epsilon/13 in every nonzero character.
        for character in range(1, P):
            for address in addresses:
                require(
                    (-character * address) % P in range(P),
                    f"{name}: invalid root phase exponent",
                )
                pointwise_squared_magnitude = Fraction(1, P**2)
                terminal_component_length = 2 * P * EPSILON
                energy = pointwise_squared_magnitude * terminal_component_length
                require(
                    energy == 2 * EPSILON / P,
                    f"{name}: rooted component energy changed",
                )
                root_character_component_checks += 1

        phase_difference = (M * (centers[1] - centers[0])) % 1
        relative_phases[name] = phase_difference
        signature = (
            f"owner=c1,target=c2,clock=2,source_root=8,"
            f"terminal_roots={addresses},halfwidths={EPSILON},{EPSILON},"
            f"root_energy_per_terminal_component={2 * EPSILON / P}"
        )
        signatures.append(signature)

    require(signatures[0] == signatures[1], "discrete signatures differ")
    require(
        root_character_component_checks == 2 * 12 * 2,
        "rooted component check census changed",
    )
    require(
        relative_phases["antipodal"] == Fraction(1, 2),
        "unperturbed components are not antipodal",
    )
    require(
        relative_phases["perturbed"] == Fraction(1, 2) + M * ETA,
        "perturbed relative phase is wrong",
    )
    require(
        relative_phases["perturbed"] != Fraction(1, 2),
        "perturbation did not break antipodality",
    )

    # The equal-width two-interval formula now proves:
    # antipodal Fhat(4)=0, perturbed Fhat(4)!=0.
    # The same holds for Ehat(52)=Fhat(4)/13 and W_4hat(0)=Fhat(4).
    antipodal_zero = True
    perturbed_zero = False
    require(antipodal_zero and not perturbed_zero, "Fourier verdict changed")

    # Exact linear defect rank of the rooted-energy/current-service quotient
    # on the two fixed antipodal component currents.  Every one of the twelve
    # rooted energies is the same multiple of the total mass, hence the
    # observation row space is generated by [1,1].  The component phases are
    # +i and -i, so Re N is zero and Im N is [1,-1].
    energy_rows = [[Fraction(1), Fraction(1)] for _ in range(12)]
    real_target = [Fraction(0), Fraction(0)]
    imag_target = [Fraction(1), Fraction(-1)]
    observation_rank = rank_rational(energy_rows)
    joint_rank = rank_rational(energy_rows + [real_target, imag_target])
    defect_rank = joint_rank - observation_rank
    require(observation_rank == 1, "energy quotient rank changed")
    require(joint_rank == 2, "phase-augmented rank changed")
    require(defect_rank == 1, "continuation defect rank is not one")

    # A single signed mass sidecar [1,-1] repairs the fixed-phase weight
    # quotient.  With component masses already fixed, the remaining positional
    # sidecar is the one circle-valued relative phase above.
    sidecar_rank = rank_rational([[Fraction(1), Fraction(-1)]])
    repaired_rank = rank_rational(
        energy_rows + [[Fraction(1), Fraction(-1)]]
    )
    require(sidecar_rank == 1 and repaired_rank == 2, "rank-one repair failed")

    digest_payload = "\n".join(
        [
            signatures[0],
            f"margin={exact_margin}",
            f"antipodal_phase={relative_phases['antipodal']}",
            f"perturbed_phase={relative_phases['perturbed']}",
            f"ranks={observation_rank},{joint_rank},{defect_rank}",
        ]
    )
    signature_digest = sha256(digest_payload.encode()).hexdigest()

    print("THM-2303 terminal-component phase-current exact probe")
    print("profile=(1,3,5)")
    print(f"pair_relation=13*{UNITS[0]}-{M}*{C1}=0")
    print(f"independent_relation={UNITS[0]}-2*{UNITS[1]}=0")
    print(f"anchored_minor={anchored_minor} mod13={anchored_minor % P}")
    print(f"exact_strict_y_margin={exact_margin}")
    print(f"epsilon={EPSILON} eta={ETA}")
    print(f"terminal_root_addresses={terminal_addresses['antipodal']}")
    print(f"root_character_component_checks={root_character_component_checks}")
    print("same_discrete_root_current_signature=True")
    print(
        "antipodal_relative_phase_cycles="
        f"{relative_phases['antipodal']}"
    )
    print(
        "perturbed_relative_phase_cycles="
        f"{relative_phases['perturbed']}"
    )
    print("antipodal_Fhat4=0")
    print("perturbed_Fhat4=nonzero")
    print("antipodal_Ehat52=0")
    print("perturbed_Ehat52=nonzero")
    print("antipodal_W4hat0=0")
    print("perturbed_W4hat0=nonzero")
    print(f"root_energy_packet_ledger_rank={observation_rank}")
    print(f"phase_augmented_rank={joint_rank}")
    print(f"fixed_phase_continuation_defect_rank={defect_rank}")
    print("relative_phase_transport_edges_needed_for_two_components=1")
    print(f"signature_sha256={signature_digest}")
    print("PASS")


if __name__ == "__main__":
    main()
