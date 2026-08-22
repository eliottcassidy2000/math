#!/usr/bin/env python3
"""Independent exact renewal-face probe for G=L^43 Norm(J).

The promoted THM-3506 proves three faces of G but leaves its z-top and
gamma=i-j-5k bottom faces open.  This companion closes exactly that fixed-map
gate without expanding G.  It reconstructs the frozen 66,146-term polynomial
J independently of THM-3506's face companion, extracts two hybrid Newton
weights, and evaluates their singleton initial form through two toric limits
of the inverse cubic.

The two hybrid weights are

    delta_6=i-j-6k=gamma-k,
    delta_8=i-j-8k=gamma-3k.

For J, both are minimized uniquely by x^66 z^76.  The c-top degeneration
reads delta_6; the target-gamma degeneration reads delta_8.  Exact Vieta
norms then give the complete faces

    in_max-z(G)       = C_z x^410 z^476,
    in_min-gamma(G)   = C_g z^271 (27 x^2 z+y^3)^205,

with explicit nonzero rational scalars.  ``require`` remains active under
``python -O``.

Scope: only the fixed polynomial G=L^43 Norm(J).  The calculation proves no
finite-sheet unit for the next norm, no fifth image equation or prime, no
all-level norm induction, and no general Jacobian-conjecture statement.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import runpy
from collections import Counter
from fractions import Fraction
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest_fraction(value: Fraction) -> str:
    payload = f"{value.numerator}/{value.denominator}".encode("ascii")
    return hashlib.sha256(payload).hexdigest()


ROOT = Path(__file__).resolve().parents[1]
GLOBAL_PROBE = ROOT / "04-computation/keller_level_three_global_norm_probe_20260816.py"
J_LEDGER_SHA256 = "9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe"

# Reconstruct J through its canonical exact route, but do not import the
# THM-3506 face script.  Its report is suppressed only to keep this output
# deterministic and focused.
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    namespace = runpy.run_path(str(GLOBAL_PROBE))
require(
    captured.getvalue().rstrip().endswith("all exact checks passed"),
    "the independent global J reconstruction did not finish cleanly",
)
J = namespace["J"]
sp = namespace["sp"]
require(namespace["coefficient_hash"](J) == J_LEDGER_SHA256, "global J ledger changed")
terms = {
    tuple(map(int, monomial)): int(coefficient)
    for monomial, coefficient in J.to_dict().items()
}
require(len(terms) == 66146, "J term count changed")


def weight(monomial: tuple[int, int, int], vector: tuple[int, int, int]) -> int:
    return sum(exponent * entry for exponent, entry in zip(monomial, vector))


def extremal_face(
    vector: tuple[int, int, int], *, take_minimum: bool
) -> tuple[int, dict[tuple[int, int, int], int], Counter[int]]:
    histogram = Counter(weight(monomial, vector) for monomial in terms)
    extreme = min(histogram) if take_minimum else max(histogram)
    face = {
        monomial: coefficient
        for monomial, coefficient in terms.items()
        if weight(monomial, vector) == extreme
    }
    return extreme, face, histogram


TOP_MONOMIAL = (66, 0, 76)
TOP_COEFFICIENT = 2**15 * 3**171
TOP_FACE = {TOP_MONOMIAL: TOP_COEFFICIENT}

z_max, z_face, _ = extremal_face((0, 0, 1), take_minimum=False)
gamma_min, gamma_face, _ = extremal_face((1, -1, -5), take_minimum=True)
delta6_min, delta6_face, delta6_histogram = extremal_face(
    (1, -1, -6), take_minimum=True
)
delta8_min, delta8_face, delta8_histogram = extremal_face(
    (1, -1, -8), take_minimum=True
)

require(z_max == 76 and z_face == TOP_FACE, "J z-top singleton changed")
require(gamma_min == -314 and len(gamma_face) == 34, "J gamma face changed")
require(
    delta6_min == gamma_min - z_max == -390,
    "delta_6 minimum misses the gamma/z bound",
)
require(
    delta8_min == gamma_min - 3 * z_max == -542,
    "delta_8 minimum misses the gamma/z bound",
)
require(delta6_face == TOP_FACE, "delta_6 has an equal-weight competitor")
require(delta8_face == TOP_FACE, "delta_8 has an equal-weight competitor")

# The inequality route is logically separate from the direct hybrid scans:
# delta_6=gamma-k and delta_8=gamma-3k.  Equality in either lower bound
# requires simultaneous equality at min(gamma) and max(k), whose intersection
# is the singleton top monomial.
gamma_support = set(gamma_face)
z_support = set(z_face)
require(gamma_support & z_support == {TOP_MONOMIAL}, "endpoint intersection changed")
require(
    all(
        weight(monomial, (1, -1, -6)) >= gamma_min - z_max
        and weight(monomial, (1, -1, -8)) >= gamma_min - 3 * z_max
        for monomial in terms
    ),
    "hybrid lower-bound hostile failed",
)

# Retain the first strict gaps as cancellation hostiles.  A zero gap would
# mean another complete equal-weight contribution had entered the resultant.
delta6_levels = sorted(delta6_histogram)
delta8_levels = sorted(delta8_histogram)
require(delta6_levels[1] > delta6_levels[0], "delta_6 strict gap vanished")
require(delta8_levels[1] > delta8_levels[0], "delta_8 strict gap vanished")


# --- Exact inverse-chart initial forms -----------------------------------

A, B, C, q = sp.symbols("A B C q")

L = 27 * A**2 * C**2 - 18 * A * B * C + 16 * A + B**3 * C - B**2
T = 4 - 3 * B * C
S = 27 * A * C**2 - 9 * B * C + 8
K = 9 * A * C - B
Y0 = 81 * A * B * C**2 - 72 * A * C - 15 * B**2 * C + 16 * B
A1 = 27 * A * B * C**2 + 54 * A * C - 9 * B**2 * C + 2 * B
A2 = (
    27 * A * B**2 * C**2
    + 18 * A * B * C
    - 48 * A
    - 9 * B**3 * C
    + 10 * B**2
)
Z0 = (
    -2916 * A**3 * C**4
    + 2916 * A**2 * B * C**3
    - 4536 * A**2 * C**2
    + 621 * A * B**3 * C**3
    - 1026 * A * B**2 * C**2
    + 504 * A * B * C
    + 64 * A
    - 207 * B**4 * C**2
    + 454 * B**3 * C
    - 256 * B**2
)


def polynomial_face(expression, vector: tuple[int, int, int], *, take_minimum: bool):
    polynomial = sp.Poly(expression, A, B, C, domain=sp.QQ)
    weighted_terms = [
        (weight(tuple(map(int, monomial)), vector), monomial, coefficient)
        for monomial, coefficient in polynomial.terms()
    ]
    extreme = (
        min(row[0] for row in weighted_terms)
        if take_minimum
        else max(row[0] for row in weighted_terms)
    )
    face = sum(
        coefficient * A**monomial[0] * B**monomial[1] * C**monomial[2]
        for row_weight, monomial, coefficient in weighted_terms
        if row_weight == extreme
    )
    return extreme, sp.expand(face)


def require_zero_mod_cubic(expression, cubic, label: str) -> None:
    numerator, _ = sp.together(expression).as_numer_denom()
    coefficient_field = sp.QQ.frac_field(A, B, C)
    remainder = sp.rem(
        sp.Poly(numerator, q, domain=coefficient_field),
        sp.Poly(cubic, q, domain=coefficient_field),
    )
    require(remainder.is_zero, f"{label}: residual cubic identity failed: {remainder}")


def require_face_ledger(actual, expected, label: str) -> None:
    require(actual.keys() == expected.keys(), f"{label}: invariant names changed")
    for name in actual:
        actual_weight, actual_face = actual[name]
        expected_weight, expected_face = expected[name]
        require(actual_weight == expected_weight, f"{label}: {name} weight changed")
        require(
            sp.Poly(sp.expand(actual_face - expected_face), A, B, C).is_zero,
            f"{label}: {name} initial form changed",
        )


invariants = {"L": L, "T": T, "S": S, "K": K, "Y0": Y0, "A1": A1, "A2": A2, "Z0": Z0}

# c-top scaling: (a,b,c)=(A,B,C*s^3), w=q/s, s -> infinity.
ctop_expected = {
    "L": (6, 27 * A**2 * C**2),
    "T": (3, -3 * B * C),
    "S": (6, 27 * A * C**2),
    "K": (3, 9 * A * C),
    "Y0": (6, 81 * A * B * C**2),
    "A1": (6, 27 * A * B * C**2),
    "A2": (6, 27 * A * B**2 * C**2),
    "Z0": (12, -2916 * A**3 * C**4),
}
ctop_faces = {
    name: polynomial_face(expression, (0, 0, 3), take_minimum=False)
    for name, expression in invariants.items()
}
require_face_ledger(ctop_faces, ctop_expected, "c-top")

ctop_cubic = 27 * A**2 * C * q**3 - 2
ctop_y = -3 * ctop_faces["K"][1] * ctop_faces["L"][1] * q**2 / (
    2 * ctop_faces["S"][1]
)
ctop_z = ctop_faces["Z0"][1] / (8 * ctop_faces["S"][1])
require_zero_mod_cubic(q * ctop_y + 1, ctop_cubic, "c-top y=-1/q")
require_zero_mod_cubic(q**3 * ctop_z + C, ctop_cubic, "c-top z=-C/q^3")

# target-gamma scaling: (a,b,c)=(A*t,B/t,C/t^5), w=q*t, t -> 0.
D = 27 * A**2 * C + B**3
gamma_expected = {
    "L": (-8, C * D),
    "T": (-6, -3 * B * C),
    "S": (-9, 27 * A * C**2),
    "K": (-4, 9 * A * C),
    "Y0": (-10, 81 * A * B * C**2),
    "A1": (-10, 27 * A * B * C**2),
    "A2": (-11, 27 * A * B**2 * C**2),
    "Z0": (-17, -2916 * A**3 * C**4 + 621 * A * B**3 * C**3),
}
gamma_faces = {
    name: polynomial_face(expression, (1, -1, -5), take_minimum=True)
    for name, expression in invariants.items()
}
require_face_ledger(gamma_faces, gamma_expected, "target-gamma")

gamma_cubic = D * q**3 - 3 * B * q - 2
gamma_y = (
    gamma_faces["Y0"][1]
    - 3 * gamma_faces["K"][1] * gamma_faces["L"][1] * q**2
) / (2 * gamma_faces["S"][1])
gamma_z = (
    gamma_faces["Z0"][1]
    + 6 * gamma_faces["L"][1] * gamma_faces["A1"][1] * q
    - 9 * gamma_faces["L"][1] * gamma_faces["A2"][1] * q**2
) / (8 * gamma_faces["S"][1])
require_zero_mod_cubic(q * gamma_y + 1, gamma_cubic, "target-gamma y=-1/q")
require_zero_mod_cubic(q**3 * gamma_z + C, gamma_cubic, "target-gamma z=-C/q^3")


# --- Sparse norm and the two complete output faces -----------------------

e_J, m_J = 43, 15
e_G, m_G = 7 * e_J - 2 * m_J, 3 * e_J - 2 * m_J
p_G = 2 * e_G - 4 * m_G // 3
d_G = 2 * e_G - 2 * m_G // 3
r_G = e_G - 2 * m_G // 3
gamma_G = -8 * e_G + 2 * m_G
require((e_G, m_G, p_G, d_G, r_G, gamma_G) == (271, 99, 410, 476, 205, -1970), "output exponents changed")

# On either toric limit the singleton J face evaluates as
#   c_J*C^76*q^(66-3*76)=c_J*C^76*q^-162.
# At the c-top, q^3=2/(27*A^2*C); at target-gamma, Vieta gives
# product(q)=2/D.  Multiplying the three J values and L^43 gives:
gamma_scalar = Fraction(TOP_COEFFICIENT**3, 2**162)
ctop_scalar = gamma_scalar * 27**205
require(gamma_scalar == Fraction(3**513, 2**117), "gamma scalar normalization changed")
require(ctop_scalar == Fraction(3**1128, 2**117), "c-top scalar normalization changed")
require(
    digest_fraction(gamma_scalar)
    == "47a8ec93d64f57847fb860261f99d4827cef1cdcef3f48efbf9c338411941676",
    "gamma scalar hash changed",
)
require(
    digest_fraction(ctop_scalar)
    == "0bf7a6d9f25117026b604cd1d9fe552cf5df804b216bf10513a0e695756b91e5",
    "c-top scalar hash changed",
)

# Weight and variable-exponent ledgers prove that these leading expressions
# are the *complete* output faces, rather than isolated surviving monomials.
require(6 * e_J + 3 * (-delta6_min) == 3 * d_G == 1428, "c-top order changed")
require(-8 * e_J + 3 * delta8_min == gamma_G == -1970, "gamma order changed")
require(2 * e_J + 2 * 162 == p_G, "c-top A exponent changed")
require(2 * e_J + 3 * 76 + 162 == d_G, "c-top C exponent changed")
require(e_J + 3 * 76 == e_G, "gamma C exponent changed")
require(e_J + 162 == r_G, "gamma D exponent changed")

semantic_payload = "\n".join(
    (
        J_LEDGER_SHA256,
        f"delta6={delta6_min}:{sorted(delta6_face.items())}",
        f"delta8={delta8_min}:{sorted(delta8_face.items())}",
        f"ctop={p_G},{d_G}:{ctop_scalar.numerator}/{ctop_scalar.denominator}",
        f"gamma={gamma_G},{e_G},{r_G}:{gamma_scalar.numerator}/{gamma_scalar.denominator}",
    )
).encode("ascii")
semantic_hash = hashlib.sha256(semantic_payload).hexdigest()

print("== independent fixed-G renewal-face probe ==")
print(f"J coefficient-ledger sha256={J_LEDGER_SHA256}; terms={len(terms)}")
print(
    f"J hybrid minima: delta_6={delta6_min}, terms={len(delta6_face)}, "
    f"next={delta6_levels[1]}; delta_8={delta8_min}, terms={len(delta8_face)}, "
    f"next={delta8_levels[1]}"
)
print(f"unique common endpoint=x^{TOP_MONOMIAL[0]}*z^{TOP_MONOMIAL[2]}; coefficient=2^15*3^171")
print("direct hybrid scan and gamma/z endpoint-intersection route: AGREE")
print("c-top inverse limit: q^3=2/(27*A^2*C), (x,y,z)~(q/s,-s/q,-C*s^6/q^3)")
print("target-gamma inverse limit: D*q^3-3*B*q-2=0, D=27*A^2*C+B^3")
print("  (x,y,z)~(q*t,-1/(q*t),-C/(q^3*t^8)); product(q)=2/D")
print(
    f"G z-top: max z={d_G}; complete face=C_z*x^{p_G}*z^{d_G}; "
    f"C_z sha256={digest_fraction(ctop_scalar)}"
)
print(
    f"G gamma-bottom: min gamma={gamma_G}; "
    f"complete face=C_g*z^{e_G}*(27*x^2*z+y^3)^{r_G}; "
    f"C_g sha256={digest_fraction(gamma_scalar)}"
)
print(f"semantic ledger sha256={semantic_hash}")
print("renewal-face gate for this fixed G: BOTH PASS")
print("scope: no next finite-sheet unit, fifth image prime, all-level induction, or general JC claim")
print("all exact checks passed")
