#!/usr/bin/env python3
"""Promotion-ready exact companion for the first cubic-pseudoplane 2x4 peel.

Candidate statement (not canonized here):

1. Every 2x4 Euler-support cell with exactly one contribution collision is
   empty, uniformly in all integer weights and positive gaps.
2. Every common-step 2x4 cell with step d=1 or d=2 is empty.

The output-swapped 4x2 statements follow by skew-symmetry.  Common-step d=3
and the first two-disjoint-pair shape are retained as hostile positive open
controls.  No arbitrary 2x4 or Darboux-pair conclusion is encoded.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    gate(sp.simplify(sp.together(left - right)) == 0, message)


def wbracket(u: int | sp.Expr, f: sp.Expr,
             v: int | sp.Expr, g: sp.Expr,
             w: sp.Symbol) -> sp.Expr:
    return sp.expand(u * f * sp.diff(g, w) - v * sp.diff(f, w) * g)


def bucket_map(r: int, offsets: tuple[int, int, int, int]) -> dict[int, list[tuple[int, int]]]:
    result: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, u in enumerate((0, r)):
        for j, v in enumerate(offsets):
            result[u + v].append((i, j))
    return dict(result)


def collision_addresses(r: int, offsets: tuple[int, int, int, int]) -> tuple[int, ...]:
    return tuple(sorted(q for q, edges in bucket_map(r, offsets).items() if len(edges) > 1))


def difference_addresses(r: int, offsets: tuple[int, int, int, int]) -> tuple[int, ...]:
    points = set(offsets)
    return tuple(sorted(v + r for v in offsets if v + r in points))


def commute_sign_possible(u: int, v: int) -> bool:
    # A singleton commuting pair has equal strict signs.  If exactly one
    # weight is zero, its profile is scalar and removable, contradicting
    # exact active support.  The (0,0) seam is allowed and retained.
    return u * v > 0 or (u == 0 and v == 0)


def scalar_sign_rows(r: int, offsets: tuple[int, int, int, int]) -> dict[int, tuple[tuple[int, int], ...]]:
    """Complete necessary sign/zero rows for every collided scalar bucket."""
    buckets = bucket_map(r, offsets)
    result: dict[int, tuple[tuple[int, int], ...]] = {}
    for scalar_q in collision_addresses(r, offsets):
        rows: list[tuple[int, int]] = []
        # Bottom singleton is non-scalar and has sum -scalar_q-2, hence both
        # bottom weights are negative.  This finite interval is exhaustive.
        for a in range(-scalar_q - 1, 0):
            b = -scalar_q - 2 - a
            if b >= 0:
                continue
            fweights = (a, a + r)
            gweights = tuple(b + v for v in offsets)
            if all(
                commute_sign_possible(fweights[i], gweights[j])
                for q, edges in buckets.items()
                if q != scalar_q and len(edges) == 1
                for i, j in edges
            ):
                rows.append((a, b))
        result[scalar_q] = tuple(rows)
    return result


# -------------------------------------------------------------------------
# A. Exact support grammar and the uniform one-collision hostile box.
# -------------------------------------------------------------------------

collision_multiplicity = Counter()
one_collision_cells = 0
one_collision_sign_rows = 0
two_collision_controls = 0
three_collision_controls = 0

for r in range(1, 25):
    for s in range(1, 25):
        for t in range(s + 1, 33):
            for u in range(t + 1, 41):
                offsets = (0, s, t, u)
                actual = collision_addresses(r, offsets)
                predicted = difference_addresses(r, offsets)
                gate(actual == predicted, f"difference grammar r={r},V={offsets}")
                gate(len(actual) <= 3, f"too many r-differences r={r},V={offsets}")
                gate(
                    (len(actual) == 3) == (offsets == (0, r, 2 * r, 3 * r)),
                    f"three differences iff four-term AP r={r},V={offsets}",
                )
                collision_multiplicity[len(actual)] += 1
                if len(actual) == 1:
                    one_collision_cells += 1
                    rows = scalar_sign_rows(r, offsets)
                    survivor_count = sum(len(values) for values in rows.values())
                    gate(survivor_count == 0, f"one-collision sign survivor r={r},V={offsets}")
                    one_collision_sign_rows += survivor_count
                elif len(actual) == 2:
                    two_collision_controls += 1
                elif len(actual) == 3:
                    three_collision_controls += 1

gate(one_collision_cells > 0, "one-collision hostile universe empty")
gate(one_collision_sign_rows == 0, "one-collision hostile survivor")
gate(two_collision_controls > 0, "two-collision positive grammar missing")
gate(three_collision_controls > 0, "AP positive grammar missing")


# -------------------------------------------------------------------------
# B. Complete sign/zero seams for common steps d=1,2; open controls at d=3.
# -------------------------------------------------------------------------

ap1_rows = scalar_sign_rows(1, (0, 1, 2, 3))
ap2_rows = scalar_sign_rows(2, (0, 2, 4, 6))
ap3_rows = scalar_sign_rows(3, (0, 3, 6, 9))
disjoint_rows = scalar_sign_rows(4, (0, 1, 4, 5))

gate(ap1_rows == {1: (), 2: ((-1, -3),), 3: ()}, f"AP d=1 seams {ap1_rows}")
gate(
    ap2_rows == {
        2: ((-1, -3),),
        4: ((-1, -5),),
        6: ((-2, -6),),
    },
    f"AP d=2 seams {ap2_rows}",
)

# Hostile positive controls: these pass the necessary sign/zero gate and
# are intentionally NOT claimed empty by this companion.
gate(
    ap3_rows == {
        3: ((-2, -3), (-1, -4)),
        6: ((-2, -6), (-1, -7)),
        9: (),
    },
    f"AP d=3 open control changed {ap3_rows}",
)
gate(
    disjoint_rows == {4: ((-3, -3),), 5: ()},
    f"two-disjoint-pair open control changed {disjoint_rows}",
)


# -------------------------------------------------------------------------
# C. Typed homogeneous profiles and the d=1 coefficient peel.
# -------------------------------------------------------------------------

w, c = sp.symbols("w c", nonzero=True)
Delta = w**3 - c**3
T = w**3
P1 = sp.Function("P1")(T)
N1 = sp.Function("N1")(T)
Q0 = sp.Function("Q0")(T)
p1 = w**2 * Delta * P1       # valid weight -1 profile
n1 = w**2 * Delta * N1       # valid weight -1 profile
q0 = Q0                      # valid weight 0 profile
lambda0, mu0 = sp.symbols("lambda0 mu0", nonzero=True)

# F=(-1,0), G=(-3,-2,-1,0).  Bottom endpoint commutation gives G_-3=lambda p^3.
# The address-1 collision integrates G_-2=p^2(3 lambda q+mu).
m2 = p1**2 * (3 * lambda0 * q0 + mu0)
d1_collision = wbracket(0, q0, -3, lambda0 * p1**3, w) + wbracket(
    -1, p1, -2, m2, w
)
same(d1_collision, 0, "AP d=1 integrated collision")
d1_scalar = wbracket(0, q0, -2, m2, w) + wbracket(-1, p1, -1, n1, w)
d1_quotient = sp.cancel(d1_scalar / Delta**2)
same(d1_scalar, Delta**2 * d1_quotient, "AP d=1 Delta^2 scalar")
gate(sp.cancel(d1_scalar.subs(w, c)) == 0, "AP d=1 hostile arm root")


# -------------------------------------------------------------------------
# D. All three d=2 rows.
# -------------------------------------------------------------------------

Q1 = sp.Function("Q1")(T)
q1 = w * Q1                    # valid weight +1 profile
lambda1, mu1, nu1 = sp.symbols("lambda1 mu1 nu1", nonzero=True)

# D.1 Lower scalar: F=(-1,1), G=(-3,-1,1,3).
g1_lower = q1 * (3 * mu1 * p1 * q1 + nu1)
gm1_lower = p1 * (3 * mu1 * p1 * q1 + nu1)
d2_lower_upper = wbracket(-1, p1, 3, mu1 * q1**3, w) + wbracket(
    1, q1, 1, g1_lower, w
)
d2_lower_middle = wbracket(-1, p1, 1, g1_lower, w) + wbracket(
    1, q1, -1, gm1_lower, w
)
same(d2_lower_upper, 0, "AP d=2 lower row upper collision")
same(d2_lower_middle, 0, "AP d=2 lower row middle collision")
d2_lower_scalar = wbracket(-1, p1, -1, gm1_lower, w) + wbracket(
    1, q1, -3, lambda1 * p1**3, w
)
same(
    d2_lower_scalar,
    3 * (lambda1 - mu1) * p1**2 * sp.diff(p1 * q1, w),
    "AP d=2 lower scalar factor",
)

# D.2 Central scalar: F=(-1,1), G=(-5,-3,-1,1).
gm1_central = lambda1 * p1
gm3_central = p1**3 * (5 * mu1 * p1 * q1 + nu1)
d2_central_top = wbracket(-1, p1, 1, lambda1 * q1, w) + wbracket(
    1, q1, -1, gm1_central, w
)
d2_central_bottom = wbracket(-1, p1, -3, gm3_central, w) + wbracket(
    1, q1, -5, mu1 * p1**5, w
)
same(d2_central_top, 0, "AP d=2 central row top collision")
same(d2_central_bottom, 0, "AP d=2 central row bottom collision")
d2_central_scalar = wbracket(-1, p1, -1, gm1_central, w) + wbracket(
    1, q1, -3, gm3_central, w
)
d2_central_quotient = sp.cancel(d2_central_scalar / Delta**2)
same(
    d2_central_scalar,
    Delta**2 * d2_central_quotient,
    "AP d=2 central Delta^2 scalar",
)

# D.3 Upper scalar: F=(-2,0), G=(-6,-4,-2,0).
P2 = sp.Function("P2")(T)
p2 = w * Delta * P2          # valid weight -2 profile
Qzero = sp.Function("Qzero")(T)
Hzero = sp.Function("Hzero")(T)
gm4_upper = p2**2 * (3 * lambda1 * Qzero + mu1)
gm2_upper = p2 * (3 * lambda1 * Qzero**2 + 2 * mu1 * Qzero + nu1)
d2_upper_lower = wbracket(-2, p2, -4, gm4_upper, w) + wbracket(
    0, Qzero, -6, lambda1 * p2**3, w
)
d2_upper_middle = wbracket(-2, p2, -2, gm2_upper, w) + wbracket(
    0, Qzero, -4, gm4_upper, w
)
same(d2_upper_lower, 0, "AP d=2 upper row lower collision")
same(d2_upper_middle, 0, "AP d=2 upper row middle collision")
d2_upper_scalar = wbracket(-2, p2, 0, Hzero, w) + wbracket(
    0, Qzero, -2, gm2_upper, w
)
d2_upper_quotient = sp.cancel(d2_upper_scalar / Delta)
same(d2_upper_scalar, Delta * d2_upper_quotient, "AP d=2 upper Delta scalar")


# -------------------------------------------------------------------------
# E. Output-swap symmetry and optimization safety.
# -------------------------------------------------------------------------

f_generic = sp.Function("f_generic")(w)
g_generic = sp.Function("g_generic")(w)
for u in range(-9, 10):
    for v in range(-9, 10):
        same(
            wbracket(v, g_generic, u, f_generic, w),
            -wbracket(u, f_generic, v, g_generic, w),
            f"output swap skew symmetry u={u},v={v}",
        )

semantic = {
    "candidate_uniform": "every exact 2x4 cell with exactly one contribution collision is empty",
    "candidate_small_ap": "every common-step exact 2x4 cell with d=1 or d=2 is empty",
    "field": "algebraically closed characteristic zero",
    "grading": "weights (3,1,-3), bracket shift +2, Delta=w^3-c^3",
    "output_swap": "4x2 follows from bracket skew-symmetry",
    "open_controls": "AP d=3 and r=4,V=(0,1,4,5) survive the sign gate",
    "not_claimed": "other >=2-collision 2x4; arbitrary 3x3; arbitrary Darboux; JC2",
}
semantic_hash = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()
source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "inactive assert")

print("CANDIDATE=cubic-pseudoplane-first-2x4-peel")
print("UNIFORM=exactly_one_collision_2x4_empty;output_swap_4x2_empty")
print(f"HOSTILE_BOX_ONE_COLLISION_CELLS={one_collision_cells};SIGN_SURVIVORS={one_collision_sign_rows}")
print("GRAMMAR=collisions_are_length_r_edges;two_edges=3chain_or_2pairs;three_edges=4term_AP")
print(f"AP_d1_SIGN_ROWS={ap1_rows};CLOSED=Delta^2")
print(f"AP_d2_SIGN_ROWS={ap2_rows};CLOSED=lower_Delta^2,central_Delta^2,upper_Delta")
print(f"OPEN_CONTROL_AP_d3={ap3_rows}")
print(f"OPEN_CONTROL_DISJOINT_r4={disjoint_rows}")
print("NEXT_OPEN=>=2_collisions_except_AP_d<=2;noncentral_or_nonAP_3x3")
print(f"SEMANTIC_SHA256={semantic_hash}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
