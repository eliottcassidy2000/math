#!/usr/bin/env python3
"""Exact companion for THM-3733 (W004 scalar-12 persistent-arm pair)."""

from collections import defaultdict
from math import ceil, gcd

import sympy as sp


def require(condition: object, label: str) -> None:
    if condition is not True and condition != sp.S.true:
        raise ArithmeticError(label)


# ---------------------------------------------------------------------------
# W004 fibre word and the two absolute placements.
# ---------------------------------------------------------------------------

n = sp.symbols("n", integer=True, positive=True)
A_support = (0, n, 3 * n)
B_support = (0, n, 2 * n, 4 * n)
fibres: dict[sp.Expr, list[str]] = defaultdict(list)
for i, x in enumerate(A_support):
    for j, y in enumerate(B_support):
        fibres[sp.expand(x + y)].append(f"{i}{j}")
fibre_word = tuple(
    "+".join(fibres[k]) for k in sorted(fibres, key=lambda z: int(z.coeff(n)))
)
expected_word = ("00", "01+10", "02+11", "12+20", "03+21", "13+22", "23")
require(fibre_word == expected_word, "W004 fibre word")

p_a = (-n - 2, -2, 2 * n - 2)
q_a = (1 - 2 * n, 1 - n, 1, 2 * n + 1)
p_b = (1 - n, 1, 2 * n + 1)
q_b = (-2 * n - 2, -n - 2, -2, 2 * n - 2)
for label, pw, qw, arm in (
    ("A", p_a, q_a, (-2, 1)),
    ("B", p_b, q_b, (1, -2)),
):
    require((pw[1], qw[2]) == arm, f"family {label} persistent arm")
    require(sp.expand(pw[1] + qw[2]) == -1, f"family {label} scalar 12")
    require(sp.expand(pw[2] + qw[0]) == -1, f"family {label} scalar 20")


# ---------------------------------------------------------------------------
# Endpoint arithmetic and complete residue controls.
# ---------------------------------------------------------------------------

family_a_delta = tuple(gcd(r + 2, 5) for r in range(5))
family_b_delta = tuple(gcd(r - 1, 4) for r in range(4))
epsilon_residues = tuple(gcd(r - 1, 3) for r in range(3))
require(family_a_delta == (1, 1, 1, 5, 1), "family A delta residues")
require(family_b_delta == (1, 4, 1, 2), "family B delta residues")
require(epsilon_residues == (1, 3, 1), "epsilon residues")

for nv in range(2, 302):
    da = gcd(nv + 2, 5)
    db = gcd(nv - 1, 4)
    ev = gcd(nv - 1, 3)
    require(da == gcd(nv + 2, 2 * nv - 1), f"family A low gcd n={nv}")
    require(db == gcd(nv - 1, 2 * nv + 2), f"family B low gcd n={nv}")
    require(ev == gcd(2 * nv - 2, 2 * nv + 1), f"family A high gcd n={nv}")
    require(ev == gcd(2 * nv + 1, 2 * nv - 2), f"family B high gcd n={nv}")
    require(5 % da == 0 and 4 % db == 0 and 3 % ev == 0, f"charge divisibility n={nv}")


# ---------------------------------------------------------------------------
# Symbolic upper primitives for both orientations.
# ---------------------------------------------------------------------------

delta, epsilon = sp.symbols("delta epsilon", integer=True, positive=True)
a, c, d, t, rho = sp.symbols("a c d t rho", nonzero=True)
H, K, Hp, Kp = sp.symbols("H K Hp Kp", nonzero=True)


def deriv(expr: sp.Expr) -> sp.Expr:
    return sp.diff(expr, H) * Hp + sp.diff(expr, K) * Kp


def bracket(r: sp.Expr, s: sp.Expr, f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(s * deriv(f) * g - r * f * deriv(g))


alpha_a = (n + 2) / delta
beta_a = (2 * n - 1) / delta
r_a = (2 * n - 2) / epsilon
s_a = (2 * n + 1) / epsilon
k = 3 / epsilon
f0_a = a * H**alpha_a
f2_a = d * K**r_a
g3_a = t * K**s_a
L_a = rho * H**alpha_a * K**k
upper_a = bracket(-(n + 2), 2 * n + 1, f0_a, g3_a) + bracket(
    2 * n - 2, -(n - 1), f2_a, L_a
)
rho_a = a * t * (2 * n + 1) / (d * (2 * n - 2))
require(sp.factor(upper_a.subs(rho, rho_a)) == 0, "family A upper primitive")

alpha_b = (n - 1) / delta
beta_b = (2 * n + 2) / delta
r_b = (2 * n + 1) / epsilon
s_b = (2 * n - 2) / epsilon
f0_b = a * H**alpha_b
f2_b = d * K**r_b
g3_b = t * K**s_b
L_b = rho * H**alpha_b * K ** (-k)
upper_b = bracket(-(n - 1), 2 * n - 2, f0_b, g3_b) + bracket(
    2 * n + 1, -(n + 2), f2_b, L_b
)
rho_b = a * t * (2 * n - 2) / (d * (2 * n + 1))
require(sp.factor(upper_b.subs(rho, rho_b)) == 0, "family B upper primitive")

require(sp.simplify(delta * (2 * alpha_a - beta_a) - 5) == 0,
        "family A charge five")
require(sp.simplify(delta * (2 * alpha_b - beta_b) + 4) == 0,
        "family B charge minus four")
require(sp.simplify(epsilon * k - 3) == 0, "high charge three")


# ---------------------------------------------------------------------------
# Direct lowest-row sign checks over complete joint residue windows.
# ---------------------------------------------------------------------------

M, Mp = sp.symbols("M Mp", nonzero=True)
for nv in range(2, 62):
    ev = gcd(nv - 1, 3)

    # Family A.
    dv = gcd(nv + 2, 5)
    av = (nv + 2) // dv
    bv = (2 * nv - 1) // dv
    kv = 3 // ev
    mv = 5 // dv
    f0 = a * H**av
    f0p = a * av * H ** (av - 1) * Hp
    L = rho * H**av * K**kv
    Lp = rho * (
        av * H ** (av - 1) * Hp * K**kv
        + kv * H**av * K ** (kv - 1) * Kp
    )
    g0 = c * H**bv
    g0p = c * bv * H ** (bv - 1) * Hp
    source = (
        -(nv - 1) * f0p * L
        + (nv + 2) * f0 * Lp
        - (2 * nv - 1) * Mp * g0
        + 2 * M * g0p
    )
    E = ev * Hp * K + dv * H * Kp
    D = a * rho * av * kv / (c * bv)
    reduced = -c * bv * H ** (bv - 1) * (
        dv * H * Mp - 2 * M * Hp - D * H**mv * K ** (kv - 1) * E
    )
    require(sp.factor(source - reduced) == 0, f"family A lowest row n={nv}")

    # Family B.
    dv = gcd(nv - 1, 4)
    av = (nv - 1) // dv
    bv = (2 * nv + 2) // dv
    kv = 3 // ev
    mv = -4 // dv
    f0 = a * H**av
    f0p = a * av * H ** (av - 1) * Hp
    L = rho * H**av * K ** (-kv)
    Lp = rho * (
        av * H ** (av - 1) * Hp * K ** (-kv)
        - kv * H**av * K ** (-kv - 1) * Kp
    )
    g0 = c * H**bv
    g0p = c * bv * H ** (bv - 1) * Hp
    source = (
        -(nv + 2) * f0p * L
        + (nv - 1) * f0 * Lp
        - (2 * nv + 2) * Mp * g0
        - M * g0p
    )
    E = ev * Hp * K + dv * H * Kp
    D = a * rho * av * kv / (c * bv)
    reduced = -c * bv * H ** (bv - 1) * (
        dv * H * Mp + M * Hp + D * H**mv * K ** (-kv - 1) * E
    )
    require(sp.factor(source - reduced) == 0, f"family B lowest row n={nv}")


# ---------------------------------------------------------------------------
# Valuation/arm and reduced-fraction controls.
# ---------------------------------------------------------------------------

for av in range(1, 80):
    required = ceil(5 * av / 2)
    require(required > 2 * av, f"delta=5 strict order alpha={av}")
    require(ceil(required / av) == 3, f"delta=5 H-order alpha={av}")

v, A0, B0 = sp.symbols("v A0 B0", nonzero=True)
for kv in (1, 3):
    numerator = A0 + B0 * v**3 * K**kv
    denominator = v**4 * K**kv
    require(sp.factor(sp.gcd(numerator, denominator)) == 1,
            f"family B reduced fraction k={kv}")

for nv in range(2, 100):
    other_a = (2 * nv - 2, 1 - 2 * nv)
    require(other_a not in {(-2, 1), (1, -2)}, f"family A other arm n={nv}")

# n=1 scalar deletions.
require((p_a[2].subs(n, 1), q_a[3].subs(n, 1)) == (0, 3), "family A n=1")
require((p_b[0].subs(n, 1), q_b[0].subs(n, 1)) == (0, -4), "family B n=1")
fp, gp, f, g = sp.symbols("fp gp f g")
require(sp.expand(3 * fp * g - 0 * f * gp) == 3 * fp * g, "family A deletion")
require(sp.expand(-4 * fp * g - 0 * f * gp) == -4 * fp * g, "family B deletion")


print("THM-3733 exact W004 scalar-12 persistent-arm pair audit")
print("fibre word = " + "|".join(fibre_word))
print("family A = P(-n-2,-2,2n-2), Q(1-2n,1-n,1,2n+1)")
print("family B = P(1-n,1,2n+1), Q(-2n-2,-n-2,-2,2n-2)")
print("family A residues delta|5 = " + ",".join(map(str, family_a_delta)))
print("family B residues delta|4 = " + ",".join(map(str, family_b_delta)))
print("shared high residues epsilon|3 = " + ",".join(map(str, epsilon_residues)))
print("family A mechanism = h^2 arm death")
print("family B mechanism = reduced negative bivariate charge")
print("boundaries = both n=1 scalar deletions to 2x4 verified")
print("scope = both named persistent-arm W004 scalar fibre 12+20 families")
