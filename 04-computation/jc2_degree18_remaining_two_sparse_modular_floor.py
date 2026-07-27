#!/usr/bin/env python3
"""Finite-place genus-floor scout for the B--W, C--D, and D--W banks."""

from __future__ import annotations

import hashlib
from pathlib import Path

from jc2_degree18_bc_algebraic_modular_floor import witness_for


CORE_SHA256 = "262d99129916bb6904507322b824830e3c6ab140c6c996f76a0d0d31f12b2d4b"


def main() -> None:
    core_path = Path(__file__).with_name(
        "jc2_degree18_bc_algebraic_modular_floor.py"
    )
    core_hash = hashlib.sha256(
        core_path.read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    if core_hash != CORE_SHA256:
        raise RuntimeError("finite-field core hash changed")
    print(f"finite_field_core_sha256={core_hash}")

    packets = {
        ("BW", "quadratic"): [
            5313800000,
            4508659468656,
            -136738899331083,
        ],
        ("BW", "sextic"): [
            5511577600000000000000000000,
            4983290602536960000000000000000,
            -6564822237254419568640000000000,
            -3094052863483309848285092659200000,
            -81862566455344350924421142159812608,
            -744088924275617882256518828471658624,
            -2973811237322720333634598763466407943,
        ],
        ("CD", "linear"): [
            22143375,
            6397664,
        ],
        ("CD", "cubic"): [
            387420489,
            -8964338040,
            54880100352,
            16544432128,
        ],
        ("DW", "linear"): [
            935886848,
            430565625,
        ],
        ("DW", "cubic"): [
            36028797018963968,
            17932072576352256,
            -1448500838400000,
            56162900390625,
        ],
    }
    total_roots = 0
    for (plane, packet), polynomial in packets.items():
        (
            prime,
            field,
            branch_degree,
            gcd_degree,
            infinity_separable,
            no_root,
            specialization,
        ) = witness_for(polynomial, plane)
        degree = len(polynomial) - 1
        total_roots += degree
        print(f"plane={plane}")
        print(f"packet={packet}")
        print(f"parameter_degree={degree}")
        print(f"witness_prime={prime}")
        print(
            "monic_parameter_mod_prime="
            + ",".join(str(value) for value in field.modulus)
        )
        print(f"residue_field_order={field.order}")
        print("cubic_leading_unit=1")
        print(f"branch_degree={branch_degree}")
        print(f"branch_derivative_gcd_degree={gcd_degree}")
        print(f"infinity_separable={int(infinity_separable)}")
        print(f"irreducible_specialization_x={specialization}")
        print(f"specialized_cubic_no_residue_root={int(no_root)}")
    print(f"candidate_roots_covered={total_roots}")
    print("status=REMAINING_TWO_SPARSE_FINITE_PLACE_GENUS_FLOOR_EXACT")


if __name__ == "__main__":
    main()
