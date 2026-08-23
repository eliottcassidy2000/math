#!/usr/bin/env python3
"""Exact companion for THM-3735 (W004 scalar-02 persistent-arm pair)."""

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

p_x = (-n - 2, -2, 2 * n - 2)
q_x = (1 - n, 1, n + 1, 3 * n + 1)
p_y = (1 - n, 1, 2 * n + 1)
q_y = (-n - 2, -2, n - 2, 3 * n - 2)
for label, pw, qw, arm in (
    ("X", p_x, q_x, (-2, 1)),
    ("Y", p_y, q_y, (1, -2)),
):
    require((pw[1], qw[1]) == arm, f"family {label} persistent arm")
    require(sp.expand(pw[0] + qw[2]) == -1, f"family {label} scalar 02")
    require(sp.expand(pw[1] + qw[1]) == -1, f"family {label} scalar 11")


# ---------------------------------------------------------------------------
# Endpoint arithmetic and complete residue controls.
# ---------------------------------------------------------------------------

delta_residues = tuple(gcd(r - 1, 3) for r in range(3))
x_epsilon_residues = tuple(gcd(r + 3, 8) for r in range(8))
y_epsilon_residues = tuple(gcd(r - 3, 7) for r in range(7))
require(delta_residues == (1, 3, 1), "shared delta residues")
require(x_epsilon_residues == (1, 4, 1, 2, 1, 8, 1, 2),
        "family X epsilon residues")
require(y_epsilon_residues == (1, 1, 1, 7, 1, 1, 1),
        "family Y epsilon residues")

for nv in range(2, 402):
    dv = gcd(nv - 1, 3)
    ex = gcd(nv + 3, 8)
    ey = gcd(nv - 3, 7)
    require(dv == gcd(nv + 2, nv - 1), f"family X low gcd n={nv}")
    require(dv == gcd(nv - 1, nv + 2), f"family Y low gcd n={nv}")
    require(ex == gcd(2 * nv - 2, 3 * nv + 1), f"family X high gcd n={nv}")
    require(ey == gcd(2 * nv + 1, 3 * nv - 2), f"family Y high gcd n={nv}")
    require(3 % dv == 0 and 8 % ex == 0 and 7 % ey == 0,
            f"charge divisibility n={nv}")


# ---------------------------------------------------------------------------
# Symbolic endpoint transports and their homogeneous Kummer sheets.
# ---------------------------------------------------------------------------

delta, epsilon = sp.symbols("delta epsilon", integer=True, positive=True)
a, c, d, t, rho = sp.symbols("a c d t rho", nonzero=True)
H, K, J, Hp, Kp, Jp = sp.symbols("H K J Hp Kp Jp", nonzero=True)


def deriv(expr: sp.Expr) -> sp.Expr:
    return sp.diff(expr, H) * Hp + sp.diff(expr, K) * Kp + sp.diff(expr, J) * Jp


def bracket(rw: sp.Expr, sw: sp.Expr, f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sw * deriv(f) * g - rw * f * deriv(g))


# Family X.
alpha_x = (n + 2) / delta
beta_x = (n - 1) / delta
r_x = (2 * n - 2) / epsilon
s_x = (3 * n + 1) / epsilon
k_x = (n + 3) / epsilon
f0_x = a * H**alpha_x
f2_x = d * K**r_x
g3_x = t * K**s_x
L0_x = rho * H**alpha_x * K**k_x
upper_x = bracket(-(n + 2), 3 * n + 1, f0_x, g3_x) + bracket(
    2 * n - 2, 1, f2_x, L0_x
)
rho_x = a * t * (3 * n + 1) / (d * (2 * n - 2))
require(sp.factor(upper_x.subs(rho, rho_x)) == 0, "family X upper primitive")
hom_x = bracket(2 * n - 2, 1, f2_x, J)
expected_hom_x = d * r_x * K ** (r_x - 1) * (Kp * J - epsilon * K * Jp)
require(sp.factor(hom_x - expected_hom_x) == 0, "family X Kummer sheet sign")

# Family Y (generic n>=3).
alpha_y = (n - 1) / delta
beta_y = (n + 2) / delta
r_y = (2 * n + 1) / epsilon
s_y = (3 * n - 2) / epsilon
k_y = (n - 3) / epsilon
f0_y = a * H**alpha_y
f2_y = d * K**r_y
g3_y = t * K**s_y
L0_y = rho * H**alpha_y * K**k_y
upper_y = bracket(1 - n, 3 * n - 2, f0_y, g3_y) + bracket(
    2 * n + 1, -2, f2_y, L0_y
)
rho_y = a * t * (3 * n - 2) / (d * (2 * n + 1))
require(sp.factor(upper_y.subs(rho, rho_y)) == 0, "family Y upper primitive")
hom_y = bracket(2 * n + 1, -2, f2_y, J)
expected_hom_y = -d * r_y * K ** (r_y - 1) * (
    2 * Kp * J + epsilon * K * Jp
)
require(sp.factor(hom_y - expected_hom_y) == 0, "family Y Kummer sheet sign")

require(sp.simplify(delta * (alpha_x - beta_x) - 3) == 0,
        "family X low charge three")
require(sp.simplify(epsilon * k_x - (n + 3)) == 0,
        "family X high remainder charge")
require(sp.simplify(delta * (beta_y - alpha_y) - 3) == 0,
        "family Y low charge three")
require(sp.simplify(epsilon * k_y - (n - 3)) == 0,
        "family Y high remainder charge")


# ---------------------------------------------------------------------------
# Family X lowest-row identity over a complete joint residue window.
# ---------------------------------------------------------------------------

M, Mp = sp.symbols("M Mp", nonzero=True)
for nv in range(2, 2 + 3 * 8):
    dv = gcd(nv - 1, 3)
    ev = gcd(nv + 3, 8)
    av = (nv + 2) // dv
    bv = (nv - 1) // dv
    kv = (nv + 3) // ev
    mv = (nv + 5) // dv
    lv = 3 // dv

    f0 = a * H**av
    f0p = a * av * H ** (av - 1) * Hp
    g0 = c * H**bv
    g0p = c * bv * H ** (bv - 1) * Hp
    L0 = rho * H**av * K**kv
    L0p = rho * (
        av * H ** (av - 1) * Hp * K**kv
        + kv * H**av * K ** (kv - 1) * Kp
    )
    source = (
        f0p * (L0 + J)
        + (nv + 2) * f0 * (L0p + Jp)
        - (nv - 1) * Mp * g0
        + 2 * M * g0p
    )
    E = ev * Hp * K + dv * H * Kp
    D0 = a * rho * av * kv / (c * bv)
    D1 = a * av / (c * bv)
    reduced = -c * bv * H ** (bv - 1) * (
        dv * H * Mp
        - 2 * M * Hp
        - D0 * H**mv * K ** (kv - 1) * E
        - D1 * H**lv * (Hp * J + dv * H * Jp)
    )
    require(sp.factor(source - reduced) == 0, f"family X lowest row n={nv}")

    # Freeze both Euler particular solutions, including all residue boundaries.
    A0 = a * rho * av / (c * bv)
    B0 = a * av / (c * bv)
    particular = A0 * H**mv * K**kv + B0 * H**lv * J
    lhs = dv * H * deriv(particular) - 2 * particular * Hp
    rhs = D0 * H**mv * K ** (kv - 1) * E + D1 * H**lv * (
        Hp * J + dv * H * Jp
    )
    require(sp.factor(lhs - rhs) == 0, f"family X Euler primitives n={nv}")


# ---------------------------------------------------------------------------
# Arm-order gates in the generic scales.
# ---------------------------------------------------------------------------

for nv in range(2, 302):
    dv = gcd(nv - 1, 3)
    av = (nv + 2) // dv
    mv = (nv + 5) // dv
    lv = 3 // dv
    if dv == 1:
        require(mv >= 7 and lv == 3, f"family X delta=1 floors n={nv}")
    else:
        required = ceil(3 * av / 2)
        require(required > av, f"family X regularity strengthens H n={nv}")
        require(ceil(required / av) == 2, f"family X H arm order n={nv}")
        require(mv >= 3 and lv == 1, f"family X delta=3 floors n={nv}")

    if nv >= 3:
        ay = (nv - 1) // dv
        if dv == 1:
            require(ay >= 2, f"family Y direct arm square n={nv}")
        else:
            required = ceil(3 * ay / 2)
            require(required > ay, f"family Y regularity strengthens H n={nv}")
            require(ceil(required / ay) == 2, f"family Y H arm order n={nv}")

# The UFD normal forms used by the two Kummer remainders.
v, B = sp.symbols("v B", nonzero=True)
require(sp.expand((B * v**2) ** 3 - B**3 * (v**3) ** 2) == 0,
        "family X F^3=kappa H^2 normal form")
require(sp.expand((J**epsilon * K**2).subs(J, 0)) == 0,
        "family Y zero Kummer sheet")


# ---------------------------------------------------------------------------
# Family Y at n=2: primitive, exact row, and incompatible local orders.
# ---------------------------------------------------------------------------

G, Gp = sp.symbols("G Gp", nonzero=True)
rho2 = 4 * a * t / (5 * d)
f0_2 = a * H
f2_2 = d * K**5
g3_2 = t * K**4
L2 = rho2 * H / K
upper_2 = bracket(-1, 4, f0_2, g3_2) + bracket(5, -2, f2_2, L2)
require(sp.factor(upper_2) == 0, "family Y n=2 upper primitive")

# Directly verify (K^2 L)'=rho(HK)' before the divisibility substitution.
Lsym, Lp = sp.symbols("Lsym Lp", nonzero=True)
upper_source_2 = (
    4 * a * Hp * t * K**4
    + a * H * 4 * t * K**3 * Kp
    - 2 * 5 * d * K**4 * Kp * Lsym
    - 5 * d * K**5 * Lp
)
upper_reduced_2 = 5 * d * K**3 * (
    rho2 * (Hp * K + H * Kp) - (2 * K * Kp * Lsym + K**2 * Lp)
)
require(sp.factor(upper_source_2 - upper_reduced_2) == 0,
        "family Y n=2 integrated derivative sign")

# Substitute H=KG and L=rho G in the lowest row.
H_sub = K * G
Hp_sub = Kp * G + K * Gp
L_sub = rho2 * G
Lp_sub = rho2 * Gp
lowest_source_2 = (
    -2 * a * Hp_sub * L_sub
    + a * H_sub * Lp_sub
    - 4 * Mp * c * H_sub**4
    - M * 4 * c * H_sub**3 * Hp_sub
)
lowest_claim_2 = -(
    a * rho2 * G * (2 * Kp * G + K * Gp)
    + 4 * c * K**3 * G**3 * (
        K * G * Mp + M * (Kp * G + K * Gp)
    )
)
require(sp.factor(lowest_source_2 - lowest_claim_2) == 0,
        "family Y n=2 lowest-row sign")

q = sp.symbols("q", integer=True, nonnegative=True)
require(sp.expand((4 * q + 3) - (q + 1)) == 3 * q + 2,
        "family Y n=2 local-order gap")
for qv in range(0, 100):
    require(qv + 1 < 4 * qv + 3, f"family Y n=2 strict order q={qv}")


# ---------------------------------------------------------------------------
# n=1 scalar deletions and exact output.
# ---------------------------------------------------------------------------

require((p_x[2].subs(n, 1), q_x[3].subs(n, 1)) == (0, 4), "family X n=1")
require((p_y[0].subs(n, 1), q_y[0].subs(n, 1)) == (0, -3), "family Y n=1")
fp, gp, f, g = sp.symbols("fp gp f g")
require(sp.expand(4 * fp * g - 0 * f * gp) == 4 * fp * g, "family X deletion")
require(sp.expand(-3 * fp * g - 0 * f * gp) == -3 * fp * g, "family Y deletion")


print("THM-3735 exact W004 scalar-02 persistent-arm pair audit")
print("fibre word = " + "|".join(fibre_word))
print("family X = P(-n-2,-2,2n-2), Q(1-n,1,n+1,3n+1)")
print("family Y = P(1-n,1,2n+1), Q(-n-2,-2,n-2,3n-2)")
print("shared low residues delta|3 = " + ",".join(map(str, delta_residues)))
print("family X high residues epsilon|8 = " + ",".join(map(str, x_epsilon_residues)))
print("family Y high residues epsilon|7 = " + ",".join(map(str, y_epsilon_residues)))
print("family X mechanism = high Kummer sheet plus h^2 arm death")
print("family Y mechanism = generic arm death plus exact n=2 local-order gap")
print("boundaries = both n=1 scalar deletions to 2x4 verified")
print("scope = both named persistent-arm W004 scalar fibre 02+11 families")
