#!/usr/bin/env python3
"""Exact Kummer/Galois clutch and Bockstein seam referee for THM-2874."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import sys


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
PINNED = {
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    COMP / "lrc14_endpoint_galois_carry_torsor_thm2857.py":
        "0bae59c9b1460f37e1879a81746154593cb0699ee13b3e5e800ba0af95ea5e4c",
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    COMP / "prime_power_convolution_physical_intertwiner_thm2870.py":
        "8db743b70bb3b7942abcdc77b4561741937435dd65792d0332855cd203fed823",
}
for path, expected in PINNED.items():
    actual = sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, f"pinned dependency changed: {path.name}")

sys.path.insert(0, str(COMP))
import lrc14_two_origin_endpoint_projective_kummer_thm2868 as frequency


P = 13
P2 = P * P
SEAM_SECTION = {3: 0, 11: 8, 7: 4}
EDGES = ((3, 8, 11), (11, 9, 7), (3, 4, 7))


def carry(q: int, h: int) -> int:
    return (q + h) // P


def lifted(state: tuple[int, int], h: int) -> tuple[int, int]:
    r, q = state
    return ((r + carry(q, h)) % P, (q + h) % P)


def encode(state: tuple[int, int]) -> int:
    r, q = state
    return (P * r + q) % P2


def all_safe_restriction(carrier, address):
    """Keep all nine literal endpoint factors safe."""
    ell = frequency.endpoint_base.REPS[address]
    intervals = [(left, right) for left, right, weight in carrier if weight]
    for index, speed in enumerate(frequency.endpoint.W):
        if index == 0:
            intervals = frequency.endpoint.subtract_comb(
                intervals, speed, 91,
                -13 - 7 * ell[index], 13 - 7 * ell[index],
            )
        else:
            intervals = frequency.endpoint.subtract_comb(
                intervals, speed, 182,
                -13 - 14 * ell[index], 13 - 14 * ell[index],
            )
    return tuple(intervals)


def matrix_rank_mod(matrix: list[list[int]], prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    pivot_row = 0
    for col in range(cols):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][col]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][col], -1, prime)
        work[pivot_row] = [
            entry * inverse % prime for entry in work[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row or not work[row][col]:
                continue
            scalar = work[row][col]
            work[row] = [
                (left - scalar * right) % prime
                for left, right in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
    return pivot_row


def main() -> None:
    require(
        tuple(SEAM_SECTION[q] for q in (3, 7, 11)) == (0, 4, 8)
        and tuple((q, h, (q + h) % P) for q, h, _ in EDGES) == EDGES,
        "seam labels changed",
    )
    require(449 % P == 7 and 546 == 3 * (2366 // P), "integer data changed")

    (
        _module, full_module, _details, _e3, _clocks, _q_pairs
    ) = frequency.allocation.build_geometry()
    unit = full_module.T // P
    target = tuple(
        coordinate + frequency.allocation.physical.SHIFT
        for coordinate in frequency.allocation.ATOM_INTERVAL
    )
    target_atom = ((*target, 1),)
    all_safe_rows = tuple(
        (
            address,
            q,
            all_safe_restriction(
                frequency.allocation.physical.overlap.shift_weighted(
                    target_atom, q * unit
                ),
                address,
            ),
        )
        for address in (frequency.ORIGIN, frequency.STEPPED)
        for q in (3, 7, 11)
    )
    require(
        all(not restricted for _address, _q, restricted in all_safe_rows),
        "literal all-nine-safe hostile ceased to be empty",
    )

    field_rows = []
    for prime, root in frequency.endpoint.MODS:
        xi = pow(root, frequency.endpoint.NN // 2366, prime)
        omega = pow(xi, 2366 // P, prime)
        require(
            pow(xi, 2366, prime) == 1
            and pow(xi, 1183, prime) == prime - 1
            and pow(omega, P, prime) == 1
            and omega != 1,
            "primitive cyclotomic roots changed",
        )

        lift_actions = 0
        cocycle_checks = 0
        for r in range(P):
            for q in range(P):
                state = (r, q)
                for h in range(P):
                    image = lifted(state, h)
                    require(
                        encode(image) == (encode(state) + h) % P2,
                        "C169 encoding failed",
                    )
                    require(
                        pow(omega, 3 * image[0], prime)
                        == pow(omega, 3 * carry(q, h), prime)
                        * pow(omega, 3 * r, prime) % prime,
                        "frequency carry response failed",
                    )
                    lift_actions += 1
                    for k in range(P):
                        direct = lifted(state, (h + k) % P)
                        expected = (
                            (direct[0] + (h + k) // P) % P,
                            direct[1],
                        )
                        require(
                            lifted(image, k) == expected,
                            "natural-lift cocycle failed",
                        )
                        cocycle_checks += 1
                iterate = state
                for _ in range(P):
                    iterate = lifted(iterate, 1)
                require(
                    iterate == ((r + 1) % P, q),
                    "L_1^13 ceased to be one carry",
                )
        require(
            lift_actions == P ** 3 and cocycle_checks == P ** 4,
            "C169 census changed",
        )

        centered_galois = tuple(
            -pow(xi, 1020, prime) * pow(omega, 3 * r, prime) % prime
            for r in range(P)
        )
        projective_frequency = tuple(
            pow(xi, (955 + 546 * r) % 2366, prime)
            for r in range(P)
        )
        clutch = -pow(xi, -65, prime) % prime
        require(
            clutch == pow(pow(xi, 13, prime), 86, prime)
            and all(
                projective_frequency[r]
                == clutch * centered_galois[r] % prime
                for r in range(P)
            ),
            "F-rational Kummer/Galois clutch failed",
        )

        seam_ratio = (
            projective_frequency[4]
            * pow(projective_frequency[8], -1, prime)
        ) % prime
        seam_derivative = 449 * (seam_ratio - 1) % prime
        require(
            seam_ratio == omega and seam_derivative != 0,
            "normalized q11-to-q7 cotangent shadow failed",
        )

        flat_edge = {}
        ancestry_edge = {}
        for source, h, target_q in EDGES:
            displacement = (
                SEAM_SECTION[target_q] - SEAM_SECTION[source]
            ) % P
            flat_edge[source, target_q] = pow(
                omega, 3 * displacement, prime
            )
            ancestry_edge[source, target_q] = pow(
                omega, 3 * carry(source, h), prime
            )
        flat_direct = flat_edge[3, 7]
        flat_via = flat_edge[3, 11] * flat_edge[11, 7] % prime
        ancestry_direct = ancestry_edge[3, 7]
        ancestry_via = (
            ancestry_edge[3, 11] * ancestry_edge[11, 7] % prime
        )
        require(
            flat_via == flat_direct
            and ancestry_via
            == pow(omega, 3, prime) * ancestry_direct % prime,
            "flat/Bockstein holonomy comparison failed",
        )

        phi3 = 1
        phi11 = (
            ancestry_edge[3, 11] * phi3
            * pow(flat_edge[3, 11], -1, prime) % prime
        )
        phi7_direct = (
            ancestry_edge[3, 7] * phi3
            * pow(flat_edge[3, 7], -1, prime) % prime
        )
        phi7_via = (
            ancestry_edge[11, 7] * phi11
            * pow(flat_edge[11, 7], -1, prime) % prime
        )
        cech_defect = (phi7_via - phi7_direct) % prime
        require(
            phi7_via
            == pow(omega, 3, prime) * phi7_direct % prime
            and cech_defect
            == (pow(omega, 3, prime) - 1) * phi7_direct % prime
            and cech_defect != 0,
            "Cech transgression vanished",
        )

        carry_derivative_chi3 = (
            449 * (pow(omega, -3, prime) - 1)
        ) % prime
        normalization = (
            carry_derivative_chi3 * pow(cech_defect, -1, prime)
        ) % prime
        require(
            carry_derivative_chi3 != 0 and normalization != 0,
            "oriented derivative comparison failed",
        )

        convolution = [
            [1 if row == (column + 3) % P else 0 for column in range(P)]
            for row in range(P)
        ]
        diagonal = [
            [1 if row == column == 3 else 0 for column in range(P)]
            for row in range(P)
        ]
        convolution_minus_one = [
            [
                (convolution[row][column] - (row == column)) % prime
                for column in range(P)
            ]
            for row in range(P)
        ]
        require(
            matrix_rank_mod(convolution, prime) == P
            and matrix_rank_mod(diagonal, prime) == 1
            and P - matrix_rank_mod(convolution_minus_one, prime) == 1,
            "same-mask rank hostile changed",
        )

        field_rows.append((
            prime,
            clutch,
            seam_ratio,
            seam_derivative,
            flat_edge[3, 11],
            flat_edge[11, 7],
            flat_edge[3, 7],
            ancestry_via,
            cech_defect,
            carry_derivative_chi3,
            normalization,
        ))

    print("THM-2874 ENDPOINT KUMMER--GALOIS BOCKSTEIN TRANSGRESSION")
    print("C169_lift_actions=2197; cocycle_compositions=28561_per_field")
    print("L1^13=one_carry; q3_frequency_fibre=13_nonzero_chi3_states")
    print(
        "projective_t_r=(-xi^-65)*(c_r-mean_c); "
        "-xi^-65=(xi^13)^86_in_Q(zeta91)"
    )
    print("seam=(q3:r0,q11:r8,q7:r4); q11_to_q7_ratio=omega")
    print("cotangent_shadow=449*(omega-1)=7*(omega-1)_mod_(13,(omega-1)^2)")
    for row in field_rows:
        print(
            f"field={row[0]}; clutch={row[1]}; seam_ratio={row[2]}; "
            f"seam_derivative={row[3]}; flat_edges={row[4:7]}; "
            f"ancestry_via={row[7]}; cech_defect={row[8]}; "
            f"carry_derivative_chi3={row[9]}; normalization={row[10]}"
        )
    print("flat_holonomy=1; ancestry_holonomy=omega^3; global_edge_gauge=0")
    print(
        "target_mask=delta3; convolution_rank=13; diagonal_rank=1; "
        "convolution_one_eigenspace_dim=1"
    )
    print(
        "deep_bit_flip=all_nine_safe; six_q_origin_restrictions_empty=1"
    )
    print(
        "scope=external_frequency_fibre_over_q3_only; "
        "q11_cancels_q7_E3_zero; no_physical_C169_extension_or_LRC14"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
