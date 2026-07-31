#!/usr/bin/env python3
"""Exact polarized E3/complement seam probe after THM-2868.

This constructs three separately physical endpoint charts on the 20-cell
THM-2847 horn.  It does not contract the two macro-truth blocks.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


PINNED = {
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
    COMP / "lrc14_prime_power_unit_mass_q11_response_thm2839.py":
        "68ae72b62b7974e4f2c2bf7d570615c8c524746978c57cf120f6372a7250ece4",
    RESULTS / "lrc14_prime_power_unit_mass_q11_response_thm2839.out":
        "495829603ea0c3944f83d7ae269cbbc5cbdec9fdc452395e96c78de8b2e7e24b",
}
for path, expected in PINNED.items():
    actual = sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, f"pinned dependency changed: {path.name}")

import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas


allocation = atlas.allocation
endpoint_base = atlas.endpoint_base
endpoint = atlas.endpoint
horn = atlas.horn
P = 13
SEAM = (
    # frequency section r, allocation q, truth block, local offset
    (0, 3, "E3", 0),
    (4, 7, "not-E3", 1),
    (8, 11, "E3", 2),
)


def complement(intervals, modulus):
    out = []
    cursor = 0
    for left, right in intervals:
        require(0 <= left < right <= modulus, "bad present interval")
        require(cursor <= left, "present intervals overlap")
        if cursor < left:
            out.append((cursor, left))
        cursor = right
    if cursor < modulus:
        out.append((cursor, modulus))
    return tuple(out)


def restricted_piece(target_atom, q, unit, block):
    ell = endpoint_base.REPS[allocation.RIGHT_ORIGIN]
    present = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
    truth = (
        present if block == "E3"
        else complement(present, endpoint_base.T)
    )
    shifted = allocation.physical.overlap.shift_weighted(
        target_atom, q * unit
    )
    return allocation.indexed_weighted_intersection(
        shifted, truth, tuple(left for left, _right in truth)
    )


def natural_lift(state, h):
    ancestry, q = state
    return (
        ancestry + (q + h) // P,
        (q + h) % P,
    )


def main() -> None:
    (
        _module, full_module, _details, e3, clocks, _q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    unit = period // P
    interval = allocation.ATOM_INTERVAL
    target = tuple(x + allocation.physical.SHIFT for x in interval)
    target_atom = ((*target, 1),)

    cells3 = atlas.full_cells(
        3, interval, target, unit, period, full_module, e3, clocks
    )
    cells11 = atlas.full_cells(
        11, interval, target, unit, period, full_module, e3, clocks
    )
    common = tuple(sorted(set(cells3) & set(cells11)))
    target7 = horn.circular_shift_interval(target, 7 * unit, period)
    horn20 = tuple(
        cell for cell in common
        if (
            not atlas.signature(
                target7, *cell, full_module, e3, clocks
            )[0]
            and all(atlas.signature(
                target7, *cell, full_module, e3, clocks
            )[1:])
        )
    )
    expected20 = tuple(
        (s, t, 1)
        for s in (0, 3, 8, 9, 12)
        for t in (5, 6, 9, 10)
    )
    require(horn20 == expected20, "20-cell E3-only horn changed")
    start = (0, 3)
    via_11 = natural_lift(start, 8)
    composed = natural_lift(via_11, 9)
    direct = natural_lift(start, 4)
    require(
        via_11 == (0, 11)
        and composed == (1, 7)
        and direct == (0, 7),
        "THM-2851 Bockstein triangle changed",
    )

    occupancy = {}
    pieces = {}
    for _r, q, block, _offset in SEAM:
        piece = restricted_piece(target_atom, q, unit, block)
        pieces[q, block] = piece
        occupancy[q, block] = int(bool(piece))
        require(
            len(piece) == 1
            and piece[0][2] == 1
            and piece[0][1] - piece[0][0] == 26444880,
            f"seam endpoint piece changed at q={q}, block={block}",
        )
        opposite = "not-E3" if block == "E3" else "E3"
        require(
            not restricted_piece(target_atom, q, unit, opposite),
            f"seam interval met both truth blocks at q={q}",
        )

    field_rows = []
    for embedding in endpoint.MODS:
        prime, root = embedding
        xi = pow(root, endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        source_value = atlas.COMMON_SOURCE[prime]

        base_piece = pieces[3, "E3"][0]
        left, right, _weight = base_piece
        alpha_left = pow(
            root, 12 * endpoint.RDIL * left % endpoint.NN, prime
        )
        alpha_right = pow(
            root, 12 * endpoint.RDIL * right % endpoint.NN, prime
        )
        lambda_left = pow(
            root, 26 * endpoint.RDIL * left % endpoint.NN, prime
        )
        lambda_right = pow(
            root, 26 * endpoint.RDIL * right % endpoint.NN, prime
        )
        inverse_difference = pow(
            (lambda_left - lambda_right) % prime, -1, prime
        )
        require(
            pow(lambda_left, 42, prime) == pow(omega, 3, prime)
            and pow(lambda_right, 42, prime) == 1,
            "seam Prony monodromy changed",
        )

        seam_u = []
        seam_v = []
        raw_pairs = []
        for r, q, block, offset in SEAM:
            formal_m = 1 + 42 * r
            actual_m = formal_m + offset
            values = []
            for m in (actual_m, actual_m + 1):
                frequency = -(12 + 26 * m)
                value = allocation.endpoint_sum(
                    pieces[q, block], frequency, embedding
                )
                require(value != 0, "seam endpoint coefficient vanished")
                values.append(source_value * value % prime)
            current, current_next = values
            split_left = (
                current_next - lambda_right * current
            ) * inverse_difference % prime
            split_minus_right = (
                lambda_left * current - current_next
            ) * inverse_difference % prime
            transported_left = (
                split_left
                * pow(pow(lambda_left, offset, prime), -1, prime)
            ) % prime
            transported_right = (
                split_minus_right
                * pow(pow(lambda_right, offset, prime), -1, prime)
            ) % prime
            seam_u.append(transported_left)
            seam_v.append(transported_right)
            raw_pairs.append(tuple(values))

        require(
            all(
                seam_u[index]
                == seam_u[0] * pow(omega, 3 * r, prime) % prime
                for index, (r, _q, _block, _offset) in enumerate(SEAM)
            )
            and seam_v == [seam_v[0]] * 3,
            "polarized seam split law failed",
        )

        compensated_u = tuple(
            seam_u[index] * pow(omega, 10 * q, prime) % prime
            for index, (_r, q, _block, _offset) in enumerate(SEAM)
        )
        compensated_v = tuple(
            seam_v[index] * pow(omega, 10 * q, prime) % prime
            for index, (_r, q, _block, _offset) in enumerate(SEAM)
        )
        require(
            len(set(compensated_u)) == 1
            and len(set(compensated_v)) == 3,
            "unique character-three seam compensation failed",
        )
        exponents = tuple((3 * r + 10 * q) % P for r, q, *_ in SEAM)
        require(
            exponents == (4, 4, 4)
            and all(q - r == 3 for r, q, *_ in SEAM),
            "diagonal seam equation changed",
        )
        rho = lambda h: pow(omega, 3 * h, prime)
        frequency_holonomy = (
            rho(9) * rho(8) * pow(rho(4), -1, prime)
        ) % prime
        carry_character = omega
        require(
            rho(9) == carry_character
            and frequency_holonomy == 1
            and carry_character != 1
            and all(
                rho(h) * rho(k) % prime == rho((h + k) % P)
                for h in range(P) for k in range(P)
            ),
            "flat-frequency versus Bockstein holonomy check failed",
        )
        field_rows.append((
            prime,
            tuple(raw_pairs),
            tuple(seam_u),
            tuple(seam_v),
            compensated_u[0],
            frequency_holonomy,
            carry_character,
        ))

    print("THM-2868 POLARIZED E3/COMPLEMENT FREQUENCY SEAM")
    print(
        f"horn_cells={len(horn20)}; first={horn20[0]}; last={horn20[-1]}"
    )
    print(
        "seam=(r,q,truth,offset)="
        f"{SEAM}; equation=q-r=3; invariant_exponent=3r+10q=4"
    )
    print(
        "endpoint_occupancy="
        f"{tuple((q, block, occupancy[q, block]) for _, q, block, _ in SEAM)}"
    )
    for (
        prime, raw, seam_u, seam_v, invariant,
        frequency_holonomy, carry_character,
    ) in field_rows:
        print(
            f"field={prime}; raw_pair_values={raw}; "
            f"U_seam={seam_u}; V_seam={seam_v}; "
            f"compensated_U={invariant}; "
            f"frequency_triangle_holonomy={frequency_holonomy}; "
            f"ancestry_carry_character={carry_character}"
        )
    print(
        "positive=the character-three branch has one exact phase reference "
        "across q3/E3, q7/not-E3, q11/E3; the trivial branch does not"
    )
    print(
        "bockstein_boundary=rho(9)=omega matches the carry leg in isolation, "
        "but rho(9)rho(8)/rho(4)=1 whereas chi(T)=omega; the full triangle "
        "proves the frequency seam factors through q and is carry-blind"
    )
    print(
        "scope=three separately physical graded charts on the 20-cell horn; "
        "no contraction between E3 and not-E3, no single positive current, "
        "no ancestry action, and no row closure"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
