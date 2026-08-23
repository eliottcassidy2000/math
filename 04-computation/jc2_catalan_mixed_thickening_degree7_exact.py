#!/usr/bin/env python3
"""Exact closure companion for the N=3, D=7 mixed Catalan cell.

The proof replays the complete verified D=6 artifact, transfers its exact
degree ledger to cap seven, and checks every newly admitted intermediate cell
against the first uniquely highest recurrence term.  The cap-seven terminal
cells are exactly the already-closed D<=6 terminal cells.

The result extends THM-3557's width-three affine closure from cap six to
cap seven.  Width at least four and width-three cap at least eight remain
open.
"""

from __future__ import annotations

import ast
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys

import sympy as sp


failures: list[str] = []
gates = 0


def gate(name: str, condition: object) -> None:
    global gates
    gates += 1
    if not bool(condition):
        failures.append(name)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(lf_bytes(path)).hexdigest()


def tied_max(degrees: list[int]) -> bool:
    top = max(degrees)
    return sum(value == top for value in degrees) >= 2


def polynomial_degree(expression: sp.Expr, variable: sp.Symbol) -> int:
    return int(sp.Poly(sp.expand(expression), variable).degree())


def coefficient(expression: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expression), variable).coeff_monomial(variable**degree)


HERE = Path(__file__).resolve()
ROOT = next(
    parent
    for parent in HERE.parents
    if (parent / "04-computation").is_dir()
    and (parent / "05-knowledge" / "results").is_dir()
)
D6_SCRIPT = ROOT / "04-computation" / "jc2_catalan_mixed_thickening_degree6_exact.py"
D6_OUTPUT = ROOT / "05-knowledge" / "results" / "jc2_catalan_mixed_thickening_degree6_exact.out"
EXPECTED_D6_SCRIPT = "a3f94ab9eb4157bec7effc15ac1f1a8ab0842814c3a4f4ae4f6705a97fe4346f"
EXPECTED_D6_OUTPUT = "02d23d80e4c50edc1323105a9d6fcc2cf5c77d729731a4be354c415251c31616"

gate("D6 script hash", sha256_lf(D6_SCRIPT) == EXPECTED_D6_SCRIPT)
gate("D6 output hash", sha256_lf(D6_OUTPUT) == EXPECTED_D6_OUTPUT)
d6_normal = subprocess.run(
    [sys.executable, str(D6_SCRIPT)],
    cwd=ROOT,
    capture_output=True,
    check=False,
    timeout=180,
)
d6_optimized = subprocess.run(
    [sys.executable, "-O", str(D6_SCRIPT)],
    cwd=ROOT,
    capture_output=True,
    check=False,
    timeout=180,
)
d6_stored = lf_bytes(D6_OUTPUT)
d6_normal_stdout = d6_normal.stdout.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
d6_optimized_stdout = d6_optimized.stdout.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
gate("D6 normal exit", d6_normal.returncode == 0)
gate("D6 optimized exit", d6_optimized.returncode == 0)
gate("D6 normal equals frozen", d6_normal_stdout == d6_stored)
gate("D6 optimized equals frozen", d6_optimized_stdout == d6_stored)
gate("D6 normal equals optimized", d6_normal_stdout == d6_optimized_stdout)
gate("D6 64 gates active", b"gates=64\n" in d6_normal_stdout)
gate("D6 verdict inherited", b"verdict=N3_D6_AFFINE_EMPTY_EXACT\n" in d6_normal_stdout)
gate("D6 no failures", d6_normal_stdout.endswith(b"failures=none\n"))


ABSENT = -10**6


def both_nonzero_positive_pair_ledger(
    cap: int, s_degree: int, e_degree: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Necessary (deg(t),deg(c)) states after transformed E2 and E3."""
    after_e2: list[tuple[int, int]] = []
    after_e3: list[tuple[int, int]] = []
    for t_degree in range(1, cap + 1):
        for c_degree in range(-1, cap + 1):
            # E2 = 3eH + (2a's-as') + (c't-2ct').
            degree_e_h = e_degree + 2
            degree_a_s = t_degree + s_degree - 2
            if 2 * (t_degree - 1) == s_degree:
                degree_a_s -= 1
            degree_c_t = ABSENT if c_degree < 0 else c_degree + t_degree - 1
            if c_degree >= 0 and c_degree == 2 * t_degree:
                degree_c_t -= 1
            if not tied_max([degree_e_h, degree_a_s, degree_c_t]):
                continue
            after_e2.append((t_degree, c_degree))

            # E3 = 2(c's-cs') + e't-3et'.
            degree_e_t = e_degree + t_degree - 1
            if e_degree == 3 * t_degree:
                degree_e_t -= 1
            degree_c_s = ABSENT if c_degree < 0 else c_degree + s_degree - 1
            if c_degree >= 0 and c_degree == s_degree:
                degree_c_s -= 1
            if tied_max([degree_e_t, degree_c_s]):
                after_e3.append((t_degree, c_degree))
    return after_e2, after_e3


def e_zero_common_pair_ledger(
    cap: int, c_degree: int, f_degree: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int]], list[tuple[int, int]]]:
    """Necessary (deg(a),deg(d)) states after E1, E2, and E3."""
    after_e1: list[tuple[int, int]] = []
    after_e2: list[tuple[int, int]] = []
    after_e3: list[tuple[int, int]] = []
    for a_degree in range(cap):
        for d_degree in range(-1, cap + 1):
            # E1 = 2(p'd-cq') + a'b-ab'.
            degree_p_d = ABSENT if d_degree < 0 else d_degree + 1
            if not tied_max([degree_p_d, c_degree + 2, 2 * a_degree]):
                continue
            after_e1.append((a_degree, d_degree))

            # E2 = 3p'f + (2a'd-ad') + (c'b-2cb').
            degree_a_d = ABSENT if d_degree < 0 else a_degree + d_degree - 1
            if d_degree >= 0 and 2 * a_degree == d_degree:
                degree_a_d -= 1
            degree_c_b = c_degree + a_degree
            if c_degree == 2 * (a_degree + 1):
                degree_c_b -= 1
            if not tied_max([f_degree + 1, degree_a_d, degree_c_b]):
                continue
            after_e2.append((a_degree, d_degree))

            # E3 = 3a'f-af' + 2(c'd-cd').
            degree_a_f = a_degree + f_degree - 1
            if 3 * a_degree == f_degree:
                degree_a_f -= 1
            degree_c_d = ABSENT if d_degree < 0 else c_degree + d_degree - 1
            if d_degree >= 0 and c_degree == d_degree:
                degree_c_d -= 1
            if tied_max([degree_a_f, degree_c_d]):
                after_e3.append((a_degree, d_degree))
    return after_e1, after_e2, after_e3


def f_zero_common_pair_ledger(
    cap: int, d_degree: int, e_degree: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Necessary (deg(a),deg(c)) states after E1 and E2."""
    after_e1: list[tuple[int, int]] = []
    after_e2: list[tuple[int, int]] = []
    for a_degree in range(cap):
        for c_degree in range(-1, cap + 1):
            degree_c_q = ABSENT if c_degree < 0 else c_degree + 2
            if not tied_max([d_degree + 1, degree_c_q, 2 * a_degree]):
                continue
            after_e1.append((a_degree, c_degree))

            degree_a_d = a_degree + d_degree - 1
            if 2 * a_degree == d_degree:
                degree_a_d -= 1
            degree_c_b = ABSENT if c_degree < 0 else c_degree + a_degree
            if c_degree >= 0 and c_degree == 2 * (a_degree + 1):
                degree_c_b -= 1
            if tied_max([e_degree + 2, degree_a_d, degree_c_b]):
                after_e2.append((a_degree, c_degree))
    return after_e1, after_e2


def one_sided_ledgers(cap: int):
    e_zero_c_zero = []
    for a_degree in range(cap):
        d_degree = 2 * a_degree - 1
        f_degree = 3 * a_degree - 3
        if 0 <= d_degree <= cap and 0 <= f_degree <= cap:
            e_zero_c_zero.append(
                (a_degree, d_degree, f_degree, 3 * a_degree - f_degree)
            )

    f_zero_d_zero = []
    for a_degree in range(1, cap):
        c_degree = 2 * a_degree - 2
        e_degree = c_degree + a_degree - 2
        if c_degree <= cap and 0 <= e_degree <= cap:
            f_zero_d_zero.append(
                (a_degree, c_degree, e_degree, e_degree - 3 * (a_degree + 1))
            )
    return e_zero_c_zero, f_zero_d_zero


def cap_ledger(cap: int) -> dict[str, object]:
    allowed_pairs = [
        (left, right)
        for left in range(cap + 1)
        for right in range(cap + 1)
        if 2 * right == 3 * left
    ]
    both_23 = both_nonzero_positive_pair_ledger(cap, 2, 3)
    both_46 = both_nonzero_positive_pair_ledger(cap, 4, 6)
    e0_23 = e_zero_common_pair_ledger(cap, 2, 3)
    e0_46 = e_zero_common_pair_ledger(cap, 4, 6)
    f0_23 = f_zero_common_pair_ledger(cap, 2, 3)
    f0_46 = f_zero_common_pair_ledger(cap, 4, 6)
    e0_c0, f0_d0 = one_sided_ledgers(cap)
    s_zero_r = [degree for degree in range(1, cap + 1) if 3 * degree <= cap]
    s_zero_after_e1 = [
        (degree, 2 * degree - 4)
        for degree in s_zero_r
        if 0 <= 2 * degree - 4 <= cap
    ]
    s_zero_after_e2 = [
        state
        for state in s_zero_after_e1
        if tied_max([3 * state[0] + 2, state[1] + state[0] - 1])
    ]
    constant_cf = [
        degree
        for degree in range(cap)
        if tied_max([degree + 1, 2, 2 * degree])
    ]
    constant_de = [
        degree
        for degree in range(cap)
        if tied_max([1, degree + 3, 2 * degree])
    ]
    return {
        "allowed": allowed_pairs,
        "both_23": both_23,
        "both_46": both_46,
        "e0_23": e0_23,
        "e0_46": e0_46,
        "f0_23": f0_23,
        "f0_46": f0_46,
        "e0_c0": e0_c0,
        "f0_d0": f0_d0,
        "s_zero": (s_zero_r, s_zero_after_e1, s_zero_after_e2),
        "constant_cf": constant_cf,
        "constant_de": constant_de,
    }


d6 = cap_ledger(6)
d7 = cap_ledger(7)
gate("D7 has no new common-power type", d7["allowed"] == [(0, 0), (2, 3), (4, 6)])
gate("D6 inherited common-power types", d6["allowed"] == d7["allowed"])

gate(
    "D7 both (2,3) new E2 cell",
    sorted(set(d7["both_23"][0]) - set(d6["both_23"][0])) == [(7, 1)],
)
gate(
    "D7 both (4,6) new E2 cells",
    sorted(set(d7["both_46"][0]) - set(d6["both_46"][0])) == [(2, 7), (7, 3)],
)
gate("D7 both (2,3) still empty after E3", d7["both_23"][1] == [])
gate("D7 both (4,6) still empty after E3", d7["both_46"][1] == [])

gate(
    "D7 e=0 (2,3) sole new E1 cell",
    sorted(set(d7["e0_23"][0]) - set(d6["e0_23"][0])) == [(4, 7)],
)
gate(
    "D7 e=0 (4,6) sole new E1 cell",
    sorted(set(d7["e0_46"][0]) - set(d6["e0_46"][0])) == [(4, 7)],
)
gate("D7 e=0 (2,3) no new E2 cell", d7["e0_23"][1] == d6["e0_23"][1])
gate("D7 e=0 (4,6) no new E2 cell", d7["e0_46"][1] == d6["e0_46"][1])
gate("D7 e=0 (2,3) inherited terminal", d7["e0_23"][2] == [(2, 3)])
gate("D7 e=0 (4,6) inherited terminal", d7["e0_46"][2] == [(3, 5)])

gate("D7 f=0 (2,3) ledger unchanged", d7["f0_23"] == d6["f0_23"])
gate("D7 f=0 (4,6) ledger unchanged", d7["f0_46"] == d6["f0_46"])
gate("D7 e=0,c=0 ledger unchanged", d7["e0_c0"] == d6["e0_c0"])
gate("D7 f=0,d=0 ledger unchanged", d7["f0_d0"] == d6["f0_d0"])
gate("D7 s=0 ledger unchanged", d7["s_zero"] == d6["s_zero"])
gate("D7 constant c,f ledger unchanged", d7["constant_cf"] == d6["constant_cf"] == [1])
gate("D7 constant d,e ledger unchanged", d7["constant_de"] == d6["constant_de"] == [3])

# Extract the exact first failed coefficient for every genuinely new cell,
# independently from the degree comparison.  Lower terms cannot affect these
# uniquely highest coefficients.
v = sp.symbols("v")


def wronskian(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(left, v) * right - left * sp.diff(right, v))


def transformed_e3(s_degree: int, e_degree: int, t_degree: int, c_degree: int) -> sp.Expr:
    s_poly = v**s_degree
    e_poly = v**e_degree
    t_poly = v**t_degree
    c_poly = v**c_degree
    return sp.expand(2 * wronskian(c_poly, s_poly) + sp.diff(e_poly, v) * t_poly - 3 * e_poly * sp.diff(t_poly, v))


new_both_23 = transformed_e3(2, 3, 7, 1)
new_both_46_low_t = transformed_e3(4, 6, 2, 7)
new_both_46_high_t = transformed_e3(4, 6, 7, 3)
gate(
    "new both (2,3;r=7,c=1) exact obstruction",
    polynomial_degree(new_both_23, v) == 9 and coefficient(new_both_23, v, 9) == -18,
)
gate(
    "new both (4,6;r=2,c=7) exact obstruction",
    polynomial_degree(new_both_46_low_t, v) == 10
    and coefficient(new_both_46_low_t, v, 10) == 6,
)
gate(
    "new both (4,6;r=7,c=3) exact obstruction",
    polynomial_degree(new_both_46_high_t, v) == 12
    and coefficient(new_both_46_high_t, v, 12) == -15,
)

p = v**2
q = v**3 - v
a4 = v**4
b5 = sp.Rational(3, 2) * v**5  # E0's forced leading ratio.
d7_poly = v**7


def e_zero_e2(c_degree: int, f_degree: int) -> sp.Expr:
    c_poly = v**c_degree
    f_poly = v**f_degree
    return sp.expand(
        3 * sp.diff(p, v) * f_poly
        + 2 * sp.diff(a4, v) * d7_poly
        - a4 * sp.diff(d7_poly, v)
        + sp.diff(c_poly, v) * b5
        - 2 * c_poly * sp.diff(b5, v)
    )


new_e0_23 = e_zero_e2(2, 3)
new_e0_46 = e_zero_e2(4, 6)
gate(
    "new e=0 (2,3;A=4,D=7) exact obstruction",
    polynomial_degree(new_e0_23, v) == 10 and coefficient(new_e0_23, v, 10) == 1,
)
gate(
    "new e=0 (4,6;A=4,D=7) exact obstruction",
    polynomial_degree(new_e0_46, v) == 10 and coefficient(new_e0_46, v, 10) == 1,
)

# Positive controls: the two terminal cells are precisely the two cells that
# the inherited D6 replay closes by exact unit ideals/triangular contradictions.
terminal_cells = [
    ("e=0", 2, 3, 2, 3),
    ("e=0", 4, 6, 3, 5),
]
gate(
    "D7 terminal cells are exactly inherited",
    terminal_cells
    == [
        ("e=0", 2, 3, *d7["e0_23"][2][0]),
        ("e=0", 4, 6, *d7["e0_46"][2][0]),
    ],
)
gate("inherited (2,3) branch exact-empty transcript", b"D6_e0_CF23=" in d6_stored and b"old_exact=EMPTY" in d6_stored)
gate("inherited (4,6) branch exact-empty transcript", b"branch_A_groebner_QQ=GroebnerBasis([1]" in d6_stored and b"branch_B_groebner_QQ=GroebnerBasis([1]" in d6_stored)

source = Path(__file__).read_text(encoding="utf-8")
gate(
    "no Python assert",
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
)

semantic = {
    "hypothesis": "characteristic zero; N=3 mixed Catalan thickening; every coefficient degree <=7",
    "inheritance": "verified D6 artifact and its complete THM-3557 replay",
    "new_both_cells": "(S,E;r,C)=(2,3;7,1),(4,6;2,7),(4,6;7,3), killed in transformed E3 by -18,6,-15",
    "new_e0_cells": "(C,F;A,D)=(2,3;4,7),(4,6;4,7), killed in E2 by the unique a,d coefficient 2A-D=1 at degree 10",
    "terminal": "exactly inherited e=0 cells (C,F;A,D)=(2,3;2,3),(4,6;3,5)",
    "conclusion": "N=3,D<=7 affine coefficient space is empty; D>=8 and N>=4 remain open",
    "projective": "not computed and not needed for affine coefficient emptiness",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("JC2 Catalan mixed-thickening degree-seven exact audit")
print(f"python={platform.python_version()} sympy={sp.__version__}")
print("d6_script_sha256_lf=" + sha256_lf(D6_SCRIPT))
print("d6_output_sha256_lf=" + sha256_lf(D6_OUTPUT))
print("d6_replay=ordinary==optimized==stored;gates=64;verdict=AFFINE_EMPTY_EXACT")
print(f"D7_allowed_common_power_pairs={d7['allowed']}")
print("D7_new_both_E2=(2,3;7,1),(4,6;2,7),(4,6;7,3)")
print("D7_new_both_E3_obstructions=-18*v^9,6*v^10,-15*v^12")
print("D7_new_e0_E1=(2,3;4,7),(4,6;4,7)")
print("D7_new_e0_E2_obstruction=(2*A-D)*a0*d0*v^10=a0*d0*v^10")
print(f"D7_terminal_cells={terminal_cells}")
print("D7_terminal_transfer=all_terminal_cells_are_D6_exact_empty")
print("verdict=N3_D7_AFFINE_EMPTY_EXACT")
print("remaining_frontier=N3_D_ge_8;N_ge_4_OPEN")
print("projective_closure=NOT_COMPUTED_NOT_NEEDED_FOR_AFFINE_COEFFICIENT_EMPTINESS")
print(f"gates={gates}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("failures=" + ("none" if not failures else ",".join(failures)))

if failures:
    raise SystemExit(1)
