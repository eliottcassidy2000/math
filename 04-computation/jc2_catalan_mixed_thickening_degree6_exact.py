#!/usr/bin/env python3
"""Exact closure of the N=3, degree-six mixed Catalan thickening cell.

For

    P = v^2 + a(v) w + c(v) w^2 + e(v) w^3,
    Q = v^3-v + b(v) w + d(v) w^2 + f(v) w^3,

this script first replays the complete THM-3557 companion normally and under
``python -O``.  It then extends the degree ledger from coefficient cap five
to cap six.  The sole degree survivor is

    e=0,  deg(a,b,c,d,f)=(3,4,4,5,6),
    c=C*h^2,  f=F*h^3,  deg(h)=2.

The survivor is closed over QQ by an exact E0/E1 parametrization, saturation
of every prescribed nonzero leading coefficient, a triangular E2/E3 split,
and two tiny saturated Groebner certificates.  Every truth gate is an active
runtime check; no Python ``assert`` is used, so ``python -O`` changes nothing.
"""

from __future__ import annotations

import hashlib
import platform
from pathlib import Path
import subprocess
import sys

import sympy as sp


failures: list[str] = []
gates = 0


def gate(name: str, condition: bool) -> None:
    """Record an optimization-stable truth gate."""
    global gates
    gates += 1
    if not bool(condition):
        failures.append(name)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(lf_bytes(path)).hexdigest()


def wronskian(x: sp.Expr, y: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(x, variable) * y - x * sp.diff(y, variable))


def jacobian_rows(
    a_rows: list[sp.Expr], b_rows: list[sp.Expr], variable: sp.Symbol
) -> list[sp.Expr]:
    if len(a_rows) != len(b_rows):
        raise ValueError("row lists must have the same length")
    width = len(a_rows) - 1
    rows: list[sp.Expr] = []
    for k in range(2 * width):
        value = 0
        for i in range(width + 1):
            j = k + 1 - i
            if 0 <= j <= width:
                value += (
                    j * sp.diff(a_rows[i], variable) * b_rows[j]
                    - i * a_rows[i] * sp.diff(b_rows[j], variable)
                )
        rows.append(sp.expand(value))
    return rows


def tied_max(degrees: list[int]) -> bool:
    top = max(degrees)
    return sum(degree == top for degree in degrees) >= 2


def unit_basis(basis: sp.GroebnerBasis) -> bool:
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


def primitive_equations(
    expressions: list[sp.Expr], generators: tuple[sp.Symbol, ...]
) -> list[sp.Expr]:
    equations: list[sp.Expr] = []
    for expression in expressions:
        if expression == 0:
            continue
        polynomial = sp.Poly(expression, *generators, domain=sp.QQ)
        equations.append(polynomial.primitive()[1].as_expr())
    return equations


# ---------------------------------------------------------------------------
# 1. Replay the complete THM-3557 proof artifact.
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
OLD_SCRIPT = ROOT / "04-computation" / "jc2_catalan_mixed_thickening_recurrence_kps_s187.py"
OLD_OUTPUT = ROOT / "05-knowledge" / "results" / "jc2_catalan_mixed_thickening_recurrence_kps_s187.out"
EXPECTED_OLD_SCRIPT_SHA256 = "0444ad61a0bb2cd165243db1a97f0cb0b299eb19263378c97d1dee9ff39a7e1e"
EXPECTED_OLD_OUTPUT_SHA256 = "cd23937963962bc4f43b83fc2f3ab6b477970c76f3151b7247c0025145380832"

gate("THM-3557 script LF hash", sha256_lf(OLD_SCRIPT) == EXPECTED_OLD_SCRIPT_SHA256)
gate("THM-3557 output LF hash", sha256_lf(OLD_OUTPUT) == EXPECTED_OLD_OUTPUT_SHA256)

old_normal = subprocess.run(
    [sys.executable, str(OLD_SCRIPT)],
    cwd=ROOT,
    capture_output=True,
    check=False,
    timeout=60,
)
old_optimized = subprocess.run(
    [sys.executable, "-O", str(OLD_SCRIPT)],
    cwd=ROOT,
    capture_output=True,
    check=False,
    timeout=60,
)
old_stored = lf_bytes(OLD_OUTPUT)
old_normal_stdout = old_normal.stdout.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
old_optimized_stdout = old_optimized.stdout.replace(b"\r\n", b"\n").replace(b"\r", b"\n")

gate("THM-3557 ordinary replay exits zero", old_normal.returncode == 0)
gate("THM-3557 optimized replay exits zero", old_optimized.returncode == 0)
gate("THM-3557 ordinary replay matches stored output", old_normal_stdout == old_stored)
gate("THM-3557 optimized replay matches stored output", old_optimized_stdout == old_stored)
gate("THM-3557 ordinary/optimized byte equality", old_normal_stdout == old_optimized_stdout)
gate("THM-3557 all 86 gates active", b"gates=86\n" in old_normal_stdout)
gate("THM-3557 old branches close", old_normal_stdout.count(b"=AFFINE_EMPTY_EXACT") == 3)
gate("THM-3557 old exact ideals are unit", old_normal_stdout.count(b"GroebnerBasis([1]") == 4)
gate("THM-3557 replay has no failures", old_normal_stdout.endswith(b"failures=none\n"))


# ---------------------------------------------------------------------------
# 2. Recurrence controls and the complete cap-six degree ledger.
# ---------------------------------------------------------------------------

v, w = sp.symbols("v w")
p = v**2
q = v**3 - v

# Positive control: the recurrence recognizes the identity map.
identity_rows = jacobian_rows([v, sp.Integer(0)], [sp.Integer(0), sp.Integer(1)], v)
gate("identity positive control", identity_rows == [1, 0])

# Hostile near-solution: preserve the genuine three-row Catalan prefix and
# locate its first leak, rather than rejecting a valid prefix prematurely.
catalan_a = [p, sp.Integer(1), sp.Rational(3, 4), sp.Rational(9, 8)]
catalan_b = [q, sp.Rational(3, 2) * v, sp.Rational(9, 8) * v,
             sp.Rational(27, 16) * v]
catalan_rows = jacobian_rows(catalan_a, catalan_b, v)
catalan_expected = [
    sp.Integer(1), sp.Integer(0), sp.Integer(0),
    -sp.Rational(135, 16), -sp.Rational(405, 64), -sp.Rational(729, 128),
]
gate("Catalan prefix hostile control", catalan_rows == catalan_expected)

CAP = 6
ABSENT = -10**6  # degree marker for the zero polynomial in ledger output


def both_nonzero_positive_pair_ledger(
    s_degree: int, e_degree: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Necessary survivors after transformed E2 and E3.

    ``c_degree=-1`` means c=0.  The internal one-degree drops are the exact
    vanishing conditions of the displayed weighted leading coefficients.
    Every final exclusion below has an exact, nonzero uniquely highest term.
    """
    after_e2: list[tuple[int, int]] = []
    after_e3: list[tuple[int, int]] = []
    for r in range(1, CAP + 1):  # r=deg(t), deg(a)=r-1
        for c_degree in range(-1, CAP + 1):
            degree_e_h = e_degree + 2
            degree_a_s = r + s_degree - 2
            if 2 * (r - 1) == s_degree:
                degree_a_s -= 1
            degree_c_t = ABSENT if c_degree < 0 else c_degree + r - 1
            if c_degree >= 0 and c_degree == 2 * r:
                degree_c_t -= 1
            if not tied_max([degree_e_h, degree_a_s, degree_c_t]):
                continue
            after_e2.append((r, c_degree))

            degree_e_t = e_degree + r - 1
            if e_degree == 3 * r:
                degree_e_t -= 1
            degree_c_s = ABSENT if c_degree < 0 else c_degree + s_degree - 1
            if c_degree >= 0 and c_degree == s_degree:
                degree_c_s -= 1
            if tied_max([degree_e_t, degree_c_s]):
                after_e3.append((r, c_degree))
    return after_e2, after_e3


allowed_s_e = [
    (s_degree, e_degree)
    for s_degree in range(CAP + 1)
    for e_degree in range(CAP + 1)
    if 2 * e_degree == 3 * s_degree
]
gate("both-top common-power degree pairs", allowed_s_e == [(0, 0), (2, 3), (4, 6)])

# If s=0 then e=mu*t^3.  At cap six r=1 dies in transformed E1;
# r=2 forces c constant there, then 3eH has degree 8 in transformed E2
# while the only other term has degree 1.
s_zero_r = [r for r in range(1, CAP + 1) if 3 * r <= CAP]
s_zero_after_e1 = [
    (r, 2 * r - 4)
    for r in s_zero_r
    if 0 <= 2 * r - 4 <= CAP
]
s_zero_after_e2 = [
    state for state in s_zero_after_e1
    if tied_max([3 * state[0] + 2, state[1] + state[0] - 1])
]
gate("both-top s=0 cap-six degrees", s_zero_r == [1, 2])
gate("both-top s=0 transformed E1 survivor", s_zero_after_e1 == [(2, 0)])
gate("both-top s=0 transformed E2 emptiness", s_zero_after_e2 == [])

both_23_e2, both_23_e3 = both_nonzero_positive_pair_ledger(2, 3)
both_46_e2, both_46_e3 = both_nonzero_positive_pair_ledger(4, 6)
gate(
    "both-top (S,E)=(2,3) transformed E2 ledger",
    both_23_e2 == [(1, 5), (3, 3), (4, 2), (5, -1), (5, 0), (5, 1), (6, 1)],
)
gate("both-top (S,E)=(2,3) transformed E3 emptiness", both_23_e3 == [])
gate(
    "both-top (S,E)=(4,6) transformed E2 ledger",
    both_46_e2 == [(4, 5), (5, 4), (6, -1), (6, 0), (6, 1), (6, 2), (6, 3)],
)
gate("both-top (S,E)=(4,6) transformed E3 emptiness", both_46_e3 == [])

# Constant nonzero s,e has transformed E2 degrees (2, <=r-2, 2r-1),
# whose maximum is unique for every admissible r.
gate(
    "both-top constant s/e emptiness",
    all(not tied_max([2, r - 2, 2 * r - 1]) for r in range(1, CAP + 1)),
)


def e_zero_common_pair_ledger(
    c_degree: int, f_degree: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int]], list[tuple[int, int]]]:
    """Necessary (deg(a),deg(d)) states after E1, E2, and E3."""
    after_e1: list[tuple[int, int]] = []
    after_e2: list[tuple[int, int]] = []
    after_e3: list[tuple[int, int]] = []
    for a_degree in range(CAP):  # deg(b)=deg(a)+1 <= 6
        for d_degree in range(-1, CAP + 1):
            degree_p_d = ABSENT if d_degree < 0 else d_degree + 1
            if not tied_max([degree_p_d, c_degree + 2, 2 * a_degree]):
                continue
            after_e1.append((a_degree, d_degree))

            degree_a_d = ABSENT if d_degree < 0 else a_degree + d_degree - 1
            if d_degree >= 0 and 2 * a_degree == d_degree:
                degree_a_d -= 1
            degree_c_b = c_degree + a_degree
            if c_degree == 2 * (a_degree + 1):
                degree_c_b -= 1
            if not tied_max([f_degree + 1, degree_a_d, degree_c_b]):
                continue
            after_e2.append((a_degree, d_degree))

            degree_a_f = a_degree + f_degree - 1
            if 3 * a_degree == f_degree:
                degree_a_f -= 1
            degree_c_d = ABSENT if d_degree < 0 else c_degree + d_degree - 1
            if d_degree >= 0 and c_degree == d_degree:
                degree_c_d -= 1
            if tied_max([degree_a_f, degree_c_d]):
                after_e3.append((a_degree, d_degree))
    return after_e1, after_e2, after_e3


e_zero_c_zero = []
for a_degree in range(CAP):
    d_degree = 2 * a_degree - 1
    f_degree = 3 * a_degree - 3
    if 0 <= d_degree <= CAP and 0 <= f_degree <= CAP:
        e_zero_c_zero.append((a_degree, d_degree, f_degree, 3 * a_degree - f_degree))
gate(
    "e=0 c=0 ledger and terminal coefficient",
    e_zero_c_zero == [(1, 1, 0, 3), (2, 3, 3, 3), (3, 5, 6, 3)],
)

e0_23_e1, e0_23_e2, e0_23_e3 = e_zero_common_pair_ledger(2, 3)
e0_46_e1, e0_46_e2, e0_46_e3 = e_zero_common_pair_ledger(4, 6)
gate(
    "e=0 (C,F)=(2,3) E3 survivor",
    e0_23_e3 == [(2, 3)],
)
gate(
    "e=0 (C,F)=(4,6) E1 ledger",
    e0_46_e1 == [(0, 5), (1, 5), (2, 5), (3, -1), (3, 0), (3, 1),
                  (3, 2), (3, 3), (3, 4), (3, 5)],
)
gate(
    "e=0 (C,F)=(4,6) E2 ledger",
    e0_46_e2 == [(3, -1), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5)],
)
gate("e=0 (C,F)=(4,6) unique E3 survivor", e0_46_e3 == [(3, 5)])

constant_c_f_a = [
    a_degree for a_degree in range(CAP)
    if tied_max([a_degree + 1, 2, 2 * a_degree])
]
gate("e=0 constant c/f old branch", constant_c_f_a == [1])


def f_zero_common_pair_ledger(
    d_degree: int, e_degree: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Necessary (deg(a),deg(c)) states after E1 and E2."""
    after_e1: list[tuple[int, int]] = []
    after_e2: list[tuple[int, int]] = []
    for a_degree in range(CAP):
        for c_degree in range(-1, CAP + 1):
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


f_zero_d_zero = []
for a_degree in range(1, CAP):
    c_degree = 2 * a_degree - 2
    e_degree = c_degree + a_degree - 2
    if c_degree <= CAP and 0 <= e_degree <= CAP:
        f_zero_d_zero.append(
            (a_degree, c_degree, e_degree, e_degree - 3 * (a_degree + 1))
        )
gate(
    "f=0 d=0 ledger and terminal coefficient",
    f_zero_d_zero == [(2, 2, 2, -7), (3, 4, 5, -7)],
)

f0_23_e1, f0_23_e2 = f_zero_common_pair_ledger(2, 3)
f0_46_e1, f0_46_e2 = f_zero_common_pair_ledger(4, 6)
gate("f=0 (D,E)=(2,3) E2 emptiness", f0_23_e2 == [])
gate("f=0 (D,E)=(4,6) E2 emptiness", f0_46_e2 == [])

constant_d_e_a = [
    a_degree for a_degree in range(CAP)
    if tied_max([1, a_degree + 3, 2 * a_degree])
]
gate("f=0 constant d/e old branch", constant_d_e_a == [3])

# The old (2,3) square/cube survivor and constant branches are exactly closed
# by the replay above.  The only new degree survivor is therefore this cell.
degree_survivors = [
    ("e=0", a_degree, a_degree + 1, 4, d_degree, 6)
    for a_degree, d_degree in e0_46_e3
]
gate(
    "complete D=6 ledger unique survivor",
    degree_survivors == [("e=0", 3, 4, 4, 5, 6)],
)


# ---------------------------------------------------------------------------
# 3. Exact E0/E1 parametrization of the sole survivor.
# ---------------------------------------------------------------------------

t, A, B, u, r, C, F = sp.symbols("t A B u r C F")
g = t * v**2 + A * v + B
h = v**2 + u * v + r  # monic; C and F retain the two leading constants
a = 1 + 2 * v * g
b = sp.Rational(3, 2) * v + (3 * v**2 - 1) * g
c = C * h**2
f = F * h**3

E0 = sp.expand(sp.diff(p, v) * b - a * sp.diff(q, v))
gate("E0 exact parametrization", E0 == 1)
gate("E0 prescribed degrees", sp.Poly(a, v).degree() == 3 and sp.Poly(b, v).degree() == 4)

d_numerator = sp.expand(2 * c * sp.diff(q, v) - wronskian(a, b, v))
d_quotient, d_remainder_poly = sp.Poly(d_numerator, v).div(sp.Poly(v, v))
d = sp.expand(d_quotient.as_expr() / 4)
d_remainder = sp.expand(d_remainder_poly.as_expr())
expected_remainder = -(
    2 * A - 4 * B**2 + 4 * C * r**2 - 3
) / 2
gate("E1 exact divisibility condition", sp.expand(d_remainder - expected_remainder) == 0)

d_lead = sp.Poly(d, v).coeff_monomial(v**5)
gate("E1 degree-five leading coefficient", sp.expand(d_lead - sp.Rational(3, 2) * (C + t**2)) == 0)

E2 = sp.expand(
    3 * sp.diff(p, v) * f
    + 2 * sp.diff(a, v) * d - a * sp.diff(d, v)
    + sp.diff(c, v) * b - 2 * c * sp.diff(b, v)
)
E3 = sp.expand(
    3 * sp.diff(a, v) * f - a * sp.diff(f, v)
    + 2 * wronskian(c, d, v)
)
E4 = sp.expand(3 * sp.diff(c, v) * f - 2 * c * sp.diff(f, v))
gate("common-power parametrization solves E4", E4 == 0)

E2_top = sp.Poly(E2, v).coeff_monomial(v**7)
E3_top = sp.Poly(E3, v).coeff_monomial(v**8)
L2 = -3 * C * t + 2 * F + t**3
L3 = 2 * F * t - C**2 - C * t**2
gate("E2 leading equation", sp.expand(E2_top - 3 * L2) == 0)
gate("E3 leading equation", sp.expand(E3_top - 3 * L3) == 0)
gate("leading perfect-square eliminant", sp.expand(t * L2 - L3 - (C - t**2)**2) == 0)

# This is the full prescribed-degree saturation before using the top rows:
# t!=0 for deg(a,b), C!=0 for deg(c), F!=0 for deg(f), and C+t^2!=0
# for deg(d)=5.  The top rows force C=t^2,F=t^3, reducing it to t!=0.
saturation_product = sp.expand(t * C * F * (C + t**2))
gate("saturation product is nonzero polynomial", saturation_product != 0)
gate(
    "top rows reduce saturation to t nonzero",
    sp.expand(saturation_product.subs({C: t**2, F: t**3}) - 2 * t**8) == 0,
)


# ---------------------------------------------------------------------------
# 4. Triangular E2/E3 closure and QQ Nullstellensatz certificates.
# ---------------------------------------------------------------------------

# Since t is invertible, C=t^2 and F=t^3 let us absorb the monic h into
# k=t*h.  Write k=t*v^2+U*v+R.  The E1 remainder is now linear in A.
U, R = sp.symbols("U R")
k = t * v**2 + U * v + R
g_reduced = t * v**2 + A * v + B
a_reduced = 1 + 2 * v * g_reduced
b_reduced = sp.Rational(3, 2) * v + (3 * v**2 - 1) * g_reduced
c_reduced = k**2
f_reduced = k**3
d_reduced_numerator = sp.expand(
    2 * c_reduced * sp.diff(q, v) - wronskian(a_reduced, b_reduced, v)
)
d_reduced_quotient, d_reduced_remainder_poly = sp.Poly(
    d_reduced_numerator, v
).div(sp.Poly(v, v))
d_reduced = sp.expand(d_reduced_quotient.as_expr() / 4)
d_reduced_remainder = sp.expand(d_reduced_remainder_poly.as_expr())
A_forced = 2 * B**2 + sp.Rational(3, 2) - 2 * R**2
gate(
    "reduced E1 remainder",
    sp.expand(d_reduced_remainder + A - A_forced) == 0,
)

E2_reduced = sp.expand(
    3 * sp.diff(p, v) * f_reduced
    + 2 * sp.diff(a_reduced, v) * d_reduced
    - a_reduced * sp.diff(d_reduced, v)
    + sp.diff(c_reduced, v) * b_reduced
    - 2 * c_reduced * sp.diff(b_reduced, v)
)
E3_reduced = sp.expand(
    3 * sp.diff(a_reduced, v) * f_reduced
    - a_reduced * sp.diff(f_reduced, v)
    + 2 * wronskian(c_reduced, d_reduced, v)
)
E2_after_E1 = sp.Poly(sp.expand(E2_reduced.subs(A, A_forced)), v)
E3_after_E1 = sp.Poly(sp.expand(E3_reduced.subs(A, A_forced)), v)
gate(
    "E2 forces common linear coefficient",
    sp.expand(E2_after_E1.coeff_monomial(v**5) - 9 * t * (A_forced - U)**2) == 0,
)
gate(
    "E3 confirms common linear coefficient",
    sp.expand(E3_after_E1.coeff_monomial(v**6) - 3 * t**2 * (A_forced - U)**2) == 0,
)

E2_middle = sp.Poly(sp.expand(E2_after_E1.as_expr().subs(U, A_forced)), v)
E3_middle = sp.Poly(sp.expand(E3_after_E1.as_expr().subs(U, A_forced)), v)
delta = B - R
split_factor = sp.expand(delta * (3 * delta + 2 * t))
gate(
    "E2 branch split",
    sp.expand(E2_middle.coeff_monomial(v**3) - 3 * t * split_factor) == 0,
)
gate(
    "E3 branch split",
    sp.expand(E3_middle.coeff_monomial(v**4) - 3 * t**2 * split_factor) == 0,
)

# Branch A: B=R.  Then A=U=3/2, and E2 has the terminal coefficient -6t^2.
E2_branch_a = sp.Poly(sp.expand(E2_middle.as_expr().subs(R, B)), v)
E3_branch_a = sp.Poly(sp.expand(E3_middle.as_expr().subs(R, B)), v)
branch_a_terminal = sp.factor(E2_branch_a.coeff_monomial(v**2))
gate("branch A terminal contradiction", branch_a_terminal == -6 * t**2)

# Branch B: 3(B-R)+2t=0, so R=B+2t/3.  Two consecutive coefficients
# demand X=81 and 2X=135, an immediate characteristic-zero contradiction.
R_shift = B + sp.Rational(2, 3) * t
E2_branch_b = sp.Poly(sp.expand(E2_middle.as_expr().subs(R, R_shift)), v)
E3_branch_b = sp.Poly(sp.expand(E3_middle.as_expr().subs(R, R_shift)), v)
X = 48 * B * t + 16 * t**2
branch_b_e2 = sp.factor(E2_branch_b.coeff_monomial(v**2))
branch_b_e3 = sp.factor(E3_branch_b.coeff_monomial(v**3))
gate(
    "branch B E2 terminal equation",
    sp.expand(branch_b_e2 - sp.Rational(2, 9) * t**2 * (X - 81)) == 0,
)
gate(
    "branch B E3 terminal equation",
    sp.expand(branch_b_e3 - sp.Rational(4, 27) * t**3 * (2 * X - 135)) == 0,
)
gate("branch B scalar incompatibility", 2 * 81 != 135)

# Independent branchwise saturated ideals over QQ.  These use every E2/E3
# coefficient after the triangular substitutions, not only the displayed
# terminal coefficients.  Their unit bases certify empty affine varieties on
# t!=0; together with the exact split_factor equation they cover the survivor.
zeta = sp.symbols("zeta")
branch_generators = (zeta, B, t)
branch_a_equations = primitive_equations(
    list(E2_branch_a.all_coeffs()) + list(E3_branch_a.all_coeffs())
    + [zeta * t - 1],
    branch_generators,
)
branch_b_equations = primitive_equations(
    list(E2_branch_b.all_coeffs()) + list(E3_branch_b.all_coeffs())
    + [zeta * t - 1],
    branch_generators,
)
G_branch_a_QQ = sp.groebner(
    branch_a_equations, *branch_generators,
    domain=sp.QQ, order="grevlex", method="f5b",
)
G_branch_b_QQ = sp.groebner(
    branch_b_equations, *branch_generators,
    domain=sp.QQ, order="grevlex", method="f5b",
)
G_branch_a_5 = sp.groebner(
    branch_a_equations, *branch_generators,
    modulus=5, order="grevlex", method="f5b",
)
G_branch_b_5 = sp.groebner(
    branch_b_equations, *branch_generators,
    modulus=5, order="grevlex", method="f5b",
)
G_branch_a_7 = sp.groebner(
    branch_a_equations, *branch_generators,
    modulus=7, order="grevlex", method="f5b",
)
G_branch_b_7 = sp.groebner(
    branch_b_equations, *branch_generators,
    modulus=7, order="grevlex", method="f5b",
)
gate("branch A saturated Groebner QQ unit ideal", unit_basis(G_branch_a_QQ))
gate("branch B saturated Groebner QQ unit ideal", unit_basis(G_branch_b_QQ))
gate("branch A GF(5) control", unit_basis(G_branch_a_5))
gate("branch B GF(5) control", unit_basis(G_branch_b_5))
gate("branch A GF(7) control", unit_basis(G_branch_a_7))
gate("branch B GF(7) control", unit_basis(G_branch_b_7))

# A concrete hostile near-solution satisfies E0,E1,E4 and all degree
# requirements, lands in branch A, and leaks by exactly the terminal terms.
near_substitution = {t: 1, B: 0, R: 0, A: sp.Rational(3, 2), U: sp.Rational(3, 2)}
near_a = sp.expand(a_reduced.subs(near_substitution))
near_b = sp.expand(b_reduced.subs(near_substitution))
near_c = sp.expand(c_reduced.subs(near_substitution))
near_d = sp.expand(d_reduced.subs(near_substitution))
near_f = sp.expand(f_reduced.subs(near_substitution))
near_rows = jacobian_rows(
    [p, near_a, near_c, sp.Integer(0)],
    [q, near_b, near_d, near_f],
    v,
)
gate("survivor hostile satisfies E0/E1", near_rows[0] == 1 and near_rows[1] == 0)
gate("survivor hostile satisfies E4/E5", near_rows[4] == 0 and near_rows[5] == 0)
gate("survivor hostile has prescribed degrees", [sp.Poly(x, v).degree() for x in
     (near_a, near_b, near_c, near_d, near_f)] == [3, 4, 4, 5, 6])
gate("survivor hostile exposes E2 leak", sp.Poly(near_rows[2], v).coeff_monomial(v**2) == -6)
gate("survivor hostile exposes E3 leak", sp.Poly(near_rows[3], v).coeff_monomial(v**3) == -4)

# Dropping t!=0 leaves the prescribed degree cell.  This hostile confirms
# that saturation removes a genuine degree-drop boundary, not a candidate.
gate(
    "degree-drop saturation hostile",
    sp.Poly(a_reduced.subs(t, 0), v).degree() < 3
    and sp.Poly(c_reduced.subs({t: 0, U: A_forced}), v).degree() < 4,
)


# ---------------------------------------------------------------------------
# 5. Stable transcript.
# ---------------------------------------------------------------------------

print("JC2 Catalan mixed-thickening degree-six exact audit")
print(f"python={platform.python_version()} sympy={sp.__version__}")
print("old_script_sha256_lf=" + sha256_lf(OLD_SCRIPT))
print("old_output_sha256_lf=" + sha256_lf(OLD_OUTPUT))
print("old_replay=ordinary==optimized==stored;gates=86;failures=none")
print("identity_control=" + ",".join(sp.sstr(x) for x in identity_rows))
print("catalan_N3_rows=" + ",".join(sp.sstr(x) for x in catalan_rows))
print(f"D6_both_allowed_SE={allowed_s_e}")
print(f"D6_both_s0=raw:{s_zero_r};afterE1:{s_zero_after_e1};afterE2:{s_zero_after_e2}")
print(f"D6_both_SE23=afterE2:{both_23_e2};afterE3:{both_23_e3}")
print(f"D6_both_SE46=afterE2:{both_46_e2};afterE3:{both_46_e3}")
print(f"D6_e0_c0={e_zero_c_zero}")
print(f"D6_e0_CF23=afterE1:{e0_23_e1};afterE2:{e0_23_e2};afterE3:{e0_23_e3};old_exact=EMPTY")
print(f"D6_e0_CF46=afterE1:{e0_46_e1};afterE2:{e0_46_e2};afterE3:{e0_46_e3}")
print(f"D6_e0_constant_cf_A={constant_c_f_a};old_exact=EMPTY")
print(f"D6_f0_d0={f_zero_d_zero}")
print(f"D6_f0_DE23=afterE1:{f0_23_e1};afterE2:{f0_23_e2}")
print(f"D6_f0_DE46=afterE1:{f0_46_e1};afterE2:{f0_46_e2}")
print(f"D6_f0_constant_de_A={constant_d_e_a};old_exact=EMPTY")
print("D6_e_equals_f_equals_0=WIDTH_TWO_EMPTY_CAP_FREE")
print(f"D6_degree_survivors={degree_survivors}")
print("E0_parameter=a=1+2*v*g;b=3*v/2+(3*v**2-1)*g;g=t*v**2+A*v+B")
print("common_power=c=C*h**2;f=F*h**3;h=v**2+u*v+r")
print("E1_remainder=" + sp.sstr(sp.factor(d_remainder)))
print("E1_d_lead=" + sp.sstr(sp.factor(d_lead)))
print("saturation=t*C*F*(C+t**2)!=0")
print("E2_E3_top=" + sp.sstr(L2) + ";" + sp.sstr(L3))
print("top_eliminant=(C-t**2)**2;forced=C=t**2,F=t**3")
print("E1_reduced_forces=A=2*B**2+3/2-2*R**2")
print("next_square=(A-U)**2;forced=U=A")
print("branch_split=(B-R)*(3*(B-R)+2*t)")
print("branch_A_terminal=" + sp.sstr(branch_a_terminal))
print("branch_B_terminal=(2*t**2/9)*(X-81),(4*t**3/27)*(2*X-135);X=48*B*t+16*t**2")
print(f"branch_A_groebner_QQ={G_branch_a_QQ}")
print(f"branch_B_groebner_QQ={G_branch_b_QQ}")
print(f"branch_A_groebner_GF5={G_branch_a_5}")
print(f"branch_B_groebner_GF5={G_branch_b_5}")
print(f"branch_A_groebner_GF7={G_branch_a_7}")
print(f"branch_B_groebner_GF7={G_branch_b_7}")
print("hostile_survivor_prefix=E0,E1,E4,E5=(1,0,0,0);E2_v2=-6;E3_v3=-4")
print("verdict=N3_D6_AFFINE_EMPTY_EXACT")
print("projective_closure=NOT_COMPUTED_NOT_NEEDED_FOR_AFFINE_COEFFICIENT_EMPTINESS")
print("full_untriangularized_groebner=NOT_USED")
print(f"gates={gates}")
print("failures=" + ("none" if not failures else ",".join(failures)))

if failures:
    raise SystemExit(1)
