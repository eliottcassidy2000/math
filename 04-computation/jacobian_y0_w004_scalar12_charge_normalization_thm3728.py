#!/usr/bin/env python3
"""Exact companion for THM-3728's W004 scalar-12 nonentry theorem.

It imports the proved THM-3613 placement ledger, checks the parity boundary,
and verifies the normalized differential identities used in the theorem.
"""

from __future__ import annotations

from hashlib import sha256
import importlib.util
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
THM3613 = ROOT / "01-canon/theorems/THM-3613-three-by-four-size-seven-ray-parity-gate.md"
SCRIPT3613 = ROOT / "04-computation/jc2_three_by_four_size_seven_ray_parity_gate_thm3613.py"
OUTPUT3613 = ROOT / "05-knowledge/results/jc2_three_by_four_size_seven_ray_parity_gate_thm3613.out"
PINS = {
    THM3613: "5b679de822d6ff202208b355a0943a3b6ea41051a7efd7d443855f6f47f1a6db",
    SCRIPT3613: "2b4423460e89696f95b5a046affeeaf36920a59f529ee08762440cbd0260daed",
    OUTPUT3613: "0b6600e995c310c5279242e2d3b69c793989b488ed697eed17117936fc7de9b8",
}


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def zero(expression, payload):
    reduced = sp.powsimp(sp.factor(sp.cancel(expression)), force=True)
    if reduced != 0:
        raise RuntimeError((payload, reduced))


def lf_sha(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for path, expected in PINS.items():
    require(lf_sha(path) == expected, ("parent hash drift", path, lf_sha(path)))

spec = importlib.util.spec_from_file_location("thm3613_w004_scout", SCRIPT3613)
require(spec is not None and spec.loader is not None, "THM-3613 import")
thm3613 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = thm3613
spec.loader.exec_module(thm3613)


# ---------------------------------------------------------------------------
# Exact placement ledger and the non-overclaim controls.
# ---------------------------------------------------------------------------

ray = thm3613.EXPECTED_RAYS["W004"]
require(ray == (1, 2, 1, 1, 2), ray)
fibres = thm3613.P.sum_fibres(ray)
require(
    fibres
    == (
        ((0, 0),),
        ((0, 1), (1, 0)),
        ((0, 2), (1, 1)),
        ((1, 2), (2, 0)),
        ((0, 3), (2, 1)),
        ((1, 3), (2, 2)),
        ((2, 3),),
    ),
    fibres,
)

tail = thm3613.tail_records(ray)
CANDIDATE_AFFINE_KEY = (3, (-3, 1), (0, -2))
THM3722_AFFINE_KEY = (4, (-3, 1), (-1, -2))
candidate = next(row for row in tail if row["affine_key"] == CANDIDATE_AFFINE_KEY)
require(candidate["eligible"] == ((2, 0, "+-"),), candidate)
require(candidate["arm_exponents"] == (2, 1), candidate)
require(next(row for row in tail if row["affine_key"] == THM3722_AFFINE_KEY), "THM-3722 family")

n2_candidate = next(
    row
    for row in thm3613.direct_scale(ray, 2)
    if row["key"] == (3, -5, -2)
)
require(
    n2_candidate["eligible"] == ((2, 0, "+-", 2, ("P0", "Q0")),),
    n2_candidate,
)
require(n2_candidate["rejected"], ("even direct n=2", n2_candidate))

for scale in (4, 6, 100):
    require(thm3613.evaluated_tail_record(candidate, scale)["rejected"], ("even parity", scale))
for scale in (3, 5, 7, 101):
    # Scale three is checked by the direct evaluator because the symbolic tail
    # starts at four.
    if scale == 3:
        row = next(
            record
            for record in thm3613.direct_scale(ray, scale)
            if record["key"] == (3, 1 - 3 * scale, -2)
        )
        require(not row["rejected"], ("odd direct", row))
    else:
        require(not thm3613.evaluated_tail_record(candidate, scale)["rejected"], ("odd tail", scale))


def surviving_keys(scale):
    return tuple(row["key"] for row in thm3613.direct_scale(ray, scale) if not row["rejected"])


def known_closed_key(scale):
    # Incoming THM-3722, now proved upstream.
    return (4, 1 - 3 * scale, -scale - 2)


def new_closed_key(scale):
    return (3, 1 - 3 * scale, -2)


small_residual_counts = {}
for scale in (1, 2, 3):
    survivors = set(surviving_keys(scale))
    survivors.discard(known_closed_key(scale))
    survivors.discard(new_closed_key(scale))
    small_residual_counts[scale] = len(survivors)
require(small_residual_counts == {1: 0, 2: 6, 3: 7}, small_residual_counts)

tail_residual_counts = {}
for scale in (4, 5):
    survivors = {
        (row["affine_key"][0],) + tuple(
            slope * scale + intercept for slope, intercept in row["affine_key"][1:]
        )
        for row in tail
        if not thm3613.evaluated_tail_record(row, scale)["rejected"]
    }
    survivors.discard(known_closed_key(scale))
    survivors.discard(new_closed_key(scale))
    tail_residual_counts[scale % 2] = len(survivors)
require(tail_residual_counts == {0: 6, 1: 8}, tail_residual_counts)

# Hostile controls: the simultaneous support reversal is W005, not W004, and
# a generic target shear Q -> Q+P^2 leaves the four-weight chart.
simultaneous_reversal = (ray[1], ray[0], ray[4], ray[3], ray[2])
require(simultaneous_reversal == (2, 1, 2, 1, 1), simultaneous_reversal)
scale = 5
p_support = (1 - 3 * scale, 1 - 2 * scale, 1)
q_support = (-2, scale - 2, 2 * scale - 2, 4 * scale - 2)
p_square_support = {left + right for left in p_support for right in p_support}
require(len(set(q_support) | p_square_support) == 10, (q_support, p_square_support))


# ---------------------------------------------------------------------------
# K-weight gauge: wedge(K^r F,K^s G)=K^(r+s) wedge(F,G).
# ---------------------------------------------------------------------------

r, s = sp.symbols("r s", integer=True)
K, Kp, F, Fp, G, Gp = sp.symbols("K Kp F Fp G Gp", nonzero=True)
raw_left_derivative = r * K ** (r - 1) * Kp * F + K**r * Fp
raw_right_derivative = s * K ** (s - 1) * Kp * G + K**s * Gp
raw_wedge = s * raw_left_derivative * K**s * G - r * K**r * F * raw_right_derivative
zero(raw_wedge - K ** (r + s) * (s * Fp * G - r * F * Gp), "K-weight gauge")


def audit_raw_transports(odd_scale):
    """Derive, rather than assume, the two K-divisibility transports."""
    alpha_value = (3 * odd_scale - 1) // 2
    A_value = 3 * odd_scale - 1
    R_value = 2 * odd_scale - 1
    C1_value = odd_scale - 2
    C2_value = 2 * odd_scale - 2
    T_value = 4 * odd_scale - 2

    U0, P0, L0, N0 = sp.symbols(f"U0 P0 L0 N0", nonzero=True)
    U0p, P0p, L0p, N0p = sp.symbols(f"U0p P0p L0p N0p")
    a0, d0, t0 = sp.symbols("a0 d0 t0", nonzero=True)

    def raw_derivative(expression):
        return (
            sp.diff(expression, K) * Kp
            + sp.diff(expression, U0) * U0p
            + sp.diff(expression, P0) * P0p
            + sp.diff(expression, L0) * L0p
            + sp.diff(expression, N0) * N0p
        )

    def raw_row(left_weight, left, right_weight, right):
        return sp.expand(
            right_weight * raw_derivative(left) * right
            - left_weight * left * raw_derivative(right)
        )

    m0 = K ** (-R_value) * P0
    f00 = a0 * K ** (-A_value) * U0**alpha_value
    f20 = d0 * K
    g30 = t0 * K**T_value
    row_y = raw_row(-R_value, m0, T_value, g30)
    row_y += raw_row(1, f20, C2_value, N0)
    expected_y = K ** (C2_value + 1) * (
        t0 * T_value * P0p
        - d0 * raw_derivative(N0 / K**C2_value)
    )
    zero(row_y - expected_y, ("raw Y transport", odd_scale))

    row_x = raw_row(-A_value, f00, T_value, g30)
    row_x += raw_row(1, f20, C1_value, L0)
    expected_x = K ** (C1_value + 1) * (
        a0 * t0 * T_value * alpha_value * U0 ** (alpha_value - 1) * U0p
        - d0 * raw_derivative(L0 / K**C1_value)
    )
    zero(row_x - expected_x, ("raw X transport", odd_scale))


for raw_scale in (3, 5, 7, 11):
    audit_raw_transports(raw_scale)


# ---------------------------------------------------------------------------
# Odd tail n=2 ell+3: two transports, a pullback row, and a diagonal Euler
# operator.  All calculations below are in the normalized coefficient ring.
# ---------------------------------------------------------------------------

ell = sp.symbols("ell", integer=True, nonnegative=True)
n = 2 * ell + 3
alpha = (3 * n - 1) / 2
A = 3 * n - 1
R = 2 * n - 1
C1 = n - 2
C2 = 2 * n - 2
T = 4 * n - 2
require(sp.expand(2 * alpha - A) == 0, "alpha typing")

U, P, X, Y = sp.symbols("U P X Y", nonzero=True)
Up, Pp, Xp, Yp = sp.symbols("Up Pp Xp Yp")
a, c, d, t = sp.symbols("a c d t", nonzero=True)
lam, mu, nu = sp.symbols("lam mu nu")


def derivative(expression):
    return (
        sp.diff(expression, U) * Up
        + sp.diff(expression, P) * Pp
        + sp.diff(expression, X) * Xp
        + sp.diff(expression, Y) * Yp
    )


def wedge(left_weight, left, right_weight, right):
    return sp.expand(
        right_weight * derivative(left) * right
        - left_weight * left * derivative(right)
    )


f0, f1, f2 = a * U**alpha, P, d
g0, g1, g2, g3 = c * U, X, Y, t
rho = t * T / d
sigma = a * rho

upper_y = wedge(-R, f1, T, g3) + wedge(1, f2, C2, g2)
upper_x = wedge(-A, f0, T, g3) + wedge(1, f2, C1, g1)
zero(upper_y - (t * T * Pp - d * Yp), "upper Y transport")
zero(upper_x - (a * t * T * alpha * U ** (alpha - 1) * Up - d * Xp), "upper X transport")

transport = {
    Y: lam + rho * P,
    Yp: rho * Pp,
    X: mu + sigma * U**alpha,
    Xp: sigma * alpha * U ** (alpha - 1) * Up,
}

middle = wedge(-A, f0, C2, g2) + wedge(-R, f1, C1, g1)
lowest = wedge(-A, f0, C1, g1) + wedge(-R, f1, -2, g0)
scalar = wedge(-R, f1, C2, g2) + wedge(1, f2, -2, g0)

Pu = sp.symbols("Pu")
D = C1 * mu + sigma * (T - 1) * U**alpha
middle_expected = Up * (
    D * Pu
    + sp.diff(D, U) * P
    + C2 * a * alpha * lam * U ** (alpha - 1)
)
lowest_expected = Up * (
    -2 * U * Pu
    + R * P
    + a * alpha * U ** (alpha - 1) * (C1 * mu + sigma * (T - 1) * U**alpha) / c
) * c
scalar_expected = rho * (T - 1) * P * Pp + C2 * lam * Pp - c * d * Up

zero(middle.subs(transport).subs(Pp, Pu * Up) - middle_expected, "middle pullback row")
zero(lowest.subs(transport).subs(Pp, Pu * Up) - lowest_expected, "lowest Euler row")
zero(scalar.subs(transport) - scalar_expected, "scalar primitive")

# The middle row integrates to D(U)P+C2*a*lam*U^alpha=nu.  Since P is an
# original polynomial in b and U is nonconstant, the polynomial-pullback
# lemma gives P in C[U].  The odd R makes the lowest Euler operator invertible
# on every integer monomial.
P_forced = a * alpha * (
    mu * U ** (alpha - 1) + sigma * U ** (2 * alpha - 1)
) / c
zero(
    lowest_expected.subs({P: P_forced, Pu: sp.diff(P_forced, U)}),
    "unique polynomial solution of lowest row",
)

middle_primitive = sp.expand(D * P_forced + C2 * a * lam * U**alpha)
top_coefficient = a * alpha * sigma**2 * (T - 1) / c
top_exponent = 3 * alpha - 1
zero(
    middle_primitive
    - (
        a * alpha * C1 * mu**2 * U ** (alpha - 1) / c
        + a * alpha * sigma * mu * (C1 + T - 1) * U ** (2 * alpha - 1) / c
        + top_coefficient * U**top_exponent
        + C2 * a * lam * U**alpha
    ),
    "middle primitive expansion",
)
require(sp.expand(top_coefficient) != 0, "nonzero unique top coefficient")

# Bounded controls make the exponent ordering and nonzero top term concrete;
# the symbolic identities above are the all-scale proof.
for odd_scale in (3, 5, 7, 9, 21):
    alpha_value = (3 * odd_scale - 1) // 2
    exponents = (alpha_value - 1, alpha_value, 2 * alpha_value - 1, 3 * alpha_value - 1)
    require(tuple(sorted(exponents)) == exponents, (odd_scale, exponents))


# ---------------------------------------------------------------------------
# Exceptional n=1: the same two transports make one lower row an exact
# derivative of a product of two nonconstant polynomials.
# ---------------------------------------------------------------------------

U1, P1 = sp.symbols("U1 P1", nonzero=True)
U1p, P1p = sp.symbols("U1p P1p")
rho1 = 2 * t / d
sigma1 = a * rho1
X1 = mu + sigma1 * U1
Y1 = lam + rho1 * P1
row_n1 = 2 * a * U1 * (rho1 * P1p) - P1p * X1 + P1 * sigma1 * U1p
zero(row_n1 - ((sigma1 * U1 - mu) * P1p + sigma1 * U1p * P1), "n=1 product derivative")


print("W004 scalar-12 anchor-20 normalized scout")
print("family=P(1-3n,1-2n,1);Q(-2,n-2,2n-2,4n-2);scalar=12+20")
print("THM3613_even_tail_rejection=True")
print("odd_n_ge_3_middle_pullback_plus_euler_contradiction=True")
print("n1_product_derivative_contradiction=True")
print(f"post_THM3722_and_new_small_residual_counts={small_residual_counts}")
print(f"post_THM3722_and_new_tail_residual_counts_by_parity={tail_residual_counts}")
print("hostile_W005_reversal_not_transported=True")
print("hostile_generic_target_shear_weight_count=10")
print("scope=one_named_W004_family_only;W004_residuals_and_W005_W006_remain")
print("ALL CHECKS PASSED")
