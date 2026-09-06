#!/usr/bin/env python3
"""Independent three-slice audit of THM-4426's two-memory companion.

The expected hashes below were frozen from a second row-state implementation
run separately on (z,h)=(s,0),(0,s),(s,s).  This compact certificate reads
the canonical clean-room output, restricts its global and boundary conics to
the same three lines, and checks all six exact expression hashes plus the two
mixed coefficients recovered by polarization.  It has no scratch dependency.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

OUTPUT = Path(
    "05-knowledge/results/"
    "jc2_source_normal_row14_weight18_two_memory_rational_conic_companion_thm4426.out"
)
Phi, eta, alpha11, c51, z, h18 = sp.symbols(
    "Phi eta alpha11 c51 z h18"
)
LOCALS = {str(v): v for v in (Phi, eta, alpha11, c51, z, h18)}
CHECKS = 0


def check(value: bool, label: str) -> None:
    global CHECKS
    if not value:
        raise AssertionError(label)
    CHECKS += 1


def expression_hash(value: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(value)).encode("ascii")).hexdigest()


def prefixed(lines: tuple[str, ...], prefix: str) -> sp.Expr:
    matches = [line[len(prefix):] for line in lines if line.startswith(prefix)]
    if len(matches) != 1:
        raise AssertionError(f"unique {prefix}")
    return sp.expand(sp.sympify(matches[0], locals=LOCALS))


def main() -> None:
    check(OUTPUT.is_file(), "canonical companion output exists")
    lines = tuple(OUTPUT.read_text().splitlines())
    global_conic = prefixed(lines, "J18zh=")
    boundary_conic = prefixed(lines, "Qzh=")

    global_slices = (
        sp.expand(global_conic.subs(h18, 0)),
        sp.expand(global_conic.subs(z, 0).subs(h18, z)),
        sp.expand(global_conic.subs(h18, z)),
    )
    expected_global_hashes = (
        "fd351a58fec40ca3210c6fa23a64ab53e60862d076efe60628421d2915d98bfb",
        "34edd2ddf542a27fc3b5c35e356d3c333d26d9fa99df631b27f8caf7007e4489",
        "ccff1193ad9d435bdefbc491c50da058d662d21f2315846f993493e0e12b69e0",
    )
    for index, (value, wanted) in enumerate(zip(global_slices, expected_global_hashes)):
        check(expression_hash(value) == wanted, f"global independent slice {index}")

    j0 = sp.expand(global_slices[0].subs(z, 0))
    check(
        expression_hash(j0)
        == "d15ce156dea6b7883378afebf23262d7950e351b845b0894986cbff39938884f",
        "common THM4415 J14 constant",
    )
    global_mixed = sp.expand(
        global_slices[2] - global_slices[0] - global_slices[1] + j0
    )

    boundary_slices = (
        sp.expand(boundary_conic.subs(h18, 0)),
        sp.expand(boundary_conic.subs(z, 0).subs(h18, z)),
        sp.expand(boundary_conic.subs(h18, z)),
    )
    expected_boundary_hashes = (
        "e2fbe8ce199c2e6a224a0d0122d2a53a315ca6d6b565032c1d0609066995ce3c",
        "9b3767eb3bc1cb0ffba2f55b123502f393af72f1bb9f7e7e99fbb093f5bed76a",
        "e1eee0d2f1d63db7c6f3b0c8e47a10c875192787ee22bd953a17846d7e4f3ba9",
    )
    for index, (value, wanted) in enumerate(zip(boundary_slices, expected_boundary_hashes)):
        check(expression_hash(value) == wanted, f"boundary independent slice {index}")

    q0 = sp.expand(boundary_slices[0].subs(z, 0))
    check(
        boundary_slices[1].subs(z, 0) == q0
        and boundary_slices[2].subs(z, 0) == q0,
        "common boundary constant",
    )
    boundary_mixed = sp.expand(
        boundary_slices[2] - boundary_slices[0] - boundary_slices[1] + q0
    )
    check(
        global_mixed
        == 40291481888765365932450000000 * z**2
        and boundary_mixed
        == 629554404511958842694531250 * z**2,
        "mixed coefficients by polarization",
    )

    print("row14_two_memory_axis_diagonal_independent_audit")
    print("method=second_row_state_implementation_on_p9_axis_p6y2_axis_and_diagonal")
    print("global_slice_hashes=" + ",".join(map(expression_hash, global_slices)))
    print("boundary_slice_hashes=" + ",".join(map(expression_hash, boundary_slices)))
    print("global_cross_hz=40291481888765365932450000000")
    print("boundary_depth_cross_hz=629554404511958842694531250")
    print("common_z0=THM4415_J14;common_boundary_constant=True")
    print("field=characteristic_zero;finite_field_used=no")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
