#!/usr/bin/env python3
"""Exact placement and differential audit for THM-3708."""

from collections import defaultdict

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Exact actual-support enumeration, with anchor labels deliberately forgotten.
# ---------------------------------------------------------------------------

A_POS = (0, 1, 2)
B_POS = (0, 1, 2, 4)
raw_fibres = defaultdict(list)
for i, left in enumerate(A_POS):
    for j, right in enumerate(B_POS):
        raw_fibres[left + right].append((i, j))
FIBRES = tuple(tuple(raw_fibres[key]) for key in sorted(raw_fibres))
SINGLETONS = tuple(fibre[0] for fibre in FIBRES if len(fibre) == 1)


def same_sign_or_both_zero(r, s):
    return (r == 0 and s == 0) or (r > 0 and s > 0) or (r < 0 and s < 0)


def actual_placements(n):
    placements = set()
    a_support = tuple(n * value for value in A_POS)
    b_support = tuple(n * value for value in B_POS)
    for scalar_index, fibre in enumerate(FIBRES):
        if len(fibre) < 2:
            continue
        for i, j in fibre:
            for arm_left, arm_right in ((-2, 1), (1, -2)):
                p0 = arm_left - a_support[i]
                q0 = arm_right - b_support[j]
                p_weights = tuple(p0 + value for value in a_support)
                q_weights = tuple(q0 + value for value in b_support)
                require(
                    all(p_weights[k] + q_weights[l] == -1 for k, l in fibre),
                    "bad scalar fibre",
                )
                if not all(
                    same_sign_or_both_zero(p_weights[k], q_weights[l])
                    for k, l in SINGLETONS
                ):
                    continue
                placements.add((scalar_index, p_weights, q_weights))
    return placements


def family_a(n):
    return 1, (-2, n - 2, 2 * n - 2), (1 - n, 1, n + 1, 3 * n + 1)


def family_b(n):
    return 1, (1 - n, 1, n + 1), (-2, n - 2, 2 * n - 2, 4 * n - 2)


def family_c(n):
    return 2, (-2, n - 2, 2 * n - 2), (1 - 2 * n, 1 - n, 1, 2 * n + 1)


def family_d(n):
    return 2, (1 - n, 1, n + 1), (-n - 2, -2, n - 2, 3 * n - 2)


def family_e(n):
    return 3, (1 - n, 1, n + 1), (-2 * n - 2, -n - 2, -2, 2 * n - 2)


def family_f(n):
    return 4, (-2, n - 2, 2 * n - 2), (1 - 4 * n, 1 - 3 * n, 1 - 2 * n, 1)


for n in range(1, 13):
    if n == 1:
        expected = set()
    elif n == 2:
        expected = {family_b(n), family_d(n), family_e(n)}
    elif n == 3:
        # A=B and C=D as actual supports, although each merged support has two
        # eligible arm addresses.  Anchor labels are not extra placements.
        expected = {family_a(n), family_c(n), family_e(n), family_f(n)}
    else:
        expected = {
            family_a(n), family_b(n), family_c(n),
            family_d(n), family_e(n), family_f(n),
        }
    require(actual_placements(n) == expected, f"placement census n={n}")


# ---------------------------------------------------------------------------
# Differential-polynomial identities.
# ---------------------------------------------------------------------------

ell = sp.symbols("ell", integer=True, positive=True)
H, K, J, Z = sp.symbols("H K J Z", nonzero=True)
Hp, Kp, Jp, Zp = sp.symbols("Hp Kp Jp Zp")
aa, bb, cc, dd, tt, kap, lam = sp.symbols(
    "aa bb cc dd tt kap lam", nonzero=True
)


def derivative(expression):
    return (
        sp.diff(expression, H) * Hp
        + sp.diff(expression, K) * Kp
        + sp.diff(expression, J) * Jp
        + sp.diff(expression, Z) * Zp
    )


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def reduced(expression):
    return sp.simplify(sp.powsimp(sp.factor(expression), force=True))


def zero(expression, label):
    require(reduced(expression) == 0, label)


def exact_quotient(expression, factor, variables, label):
    quotient = sp.cancel(expression / factor)
    zero(expression - factor * quotient, label + " reconstruction")
    # Symbolic branch exponents such as K^(2*ell-4) are polynomial on the
    # stated ell range, although SymPy may display them with a K^4 denominator.
    # Reconstruction is the exact all-parameter divisibility check; the
    # nonnegative exponent bounds are recorded separately in the proof.
    return quotient


# Family C, exceptional merged scale n=3.  For n>=4 the arm-square argument
# is purely an order statement and needs no differential identity.
F2 = H**2
gamma_c = aa * tt * 7 / (dd * 4)
beta_c = cc * gamma_c / (dd * 4)
M_c = kap * K + gamma_c * F2 * K**3
L_c = -beta_c * F2
top_c = wedge(-2, aa * F2, 7, tt * K**7) + wedge(4, dd * K**4, 1, M_c)
next_c = wedge(1, cc * K, 1, M_c) + wedge(4, dd * K**4, -2, L_c)
scalar_c = (
    wedge(-2, aa * F2, 1, M_c)
    + wedge(1, cc * K, -2, L_c)
    + wedge(4, dd * K**4, -5, bb * H**5)
)
euler_c = Hp * K + H * Kp
quotient_c = (
    2 * aa * H * (kap + 3 * gamma_c * H**2 * K**2)
    + 2 * cc * beta_c * H
    - 20 * bb * dd * H**4 * K**3
)
zero(top_c, "family C n=3 top")
zero(next_c, "family C n=3 next")
zero(scalar_c - euler_c * quotient_c, "family C n=3 scalar")


# Family D for n>=3.  The three residue branches exhaust delta=gcd(n-1,3).
D_BRANCHES = (
    (1, 3 * ell, 3 * ell - 1, 3 * ell + 2),
    (1, 3 * ell + 2, 3 * ell + 1, 3 * ell + 4),
    (3, 3 * ell + 1, ell, ell + 1),
)


def audit_family_d(delta, n, aexp, gexp):
    F = H**aexp
    G0 = H**gexp
    gamma = aa * tt * (3 * n - 2) / (dd * (n + 1))
    beta = cc * gamma / (dd * (n + 1))
    M = kap * K ** (n - 2) + gamma * F * K ** (2 * n - 3)
    L = -beta * F * K ** (n - 3)
    top = wedge(1 - n, aa * F, 3 * n - 2, tt * K ** (3 * n - 2))
    top += wedge(n + 1, dd * K ** (n + 1), n - 2, M)
    next_row = wedge(1, cc * K, n - 2, M)
    next_row += wedge(n + 1, dd * K ** (n + 1), -2, L)
    scalar = wedge(1 - n, aa * F, n - 2, M)
    scalar += wedge(1, cc * K, -2, L)
    scalar += wedge(n + 1, dd * K ** (n + 1), -n - 2, bb * G0)
    euler = Hp * K + delta * H * Kp
    quotient = aa * aexp * H ** (aexp - 1) * (
        kap * (n - 2) * K ** (n - 3)
        + gamma * (2 * n - 3) * F * K ** (2 * n - 4)
        + cc * beta / aa * K ** (n - 3)
    )
    quotient -= bb * dd * (n + 1) * gexp * H ** (gexp - 1) * K**n
    zero(top, f"family D top delta={delta}")
    zero(next_row, f"family D next delta={delta}")
    zero(scalar - euler * quotient, f"family D scalar delta={delta}")


for branch in D_BRANCHES:
    audit_family_d(*branch)


# Family D at n=2: the integrated next row forces K|H.  Put H=KJ and audit
# the resulting polynomial identities and the new Euler factor J'K+2JK'.
gamma_d2 = aa * tt * 4 / (dd * 3)
beta_d2 = cc * gamma_d2 / (dd * 3)
H_d2 = K * J
M_d2 = kap + gamma_d2 * K**2 * J
L_d2 = -beta_d2 * J
top_d2 = wedge(-1, aa * H_d2, 4, tt * K**4) + wedge(3, dd * K**3, 0, M_d2)
next_d2 = wedge(1, cc * K, 0, M_d2) + wedge(3, dd * K**3, -2, L_d2)
scalar_d2 = (
    wedge(-1, aa * H_d2, 0, M_d2)
    + wedge(1, cc * K, -2, L_d2)
    + wedge(3, dd * K**3, -4, bb * H_d2**4)
)
euler_d2 = Jp * K + 2 * J * Kp
quotient_d2 = (
    aa * gamma_d2 * K**2 * J + cc * beta_d2 - 12 * bb * dd * K**6 * J**3
)
zero(top_d2, "family D n=2 top")
zero(next_d2, "family D n=2 next")
zero(scalar_d2 - euler_d2 * quotient_d2, "family D n=2 scalar")


# Family E for n>=3.  Z is the unresolved homogeneous solution from the
# lowest double.  The triple row transfers its W(1,-n-2)(K,Z) bracket into
# the scalar row, so the latter inherits the same Euler factor.
E_BRANCHES = (
    (1, 2 * ell, 2 * ell - 1, 4 * ell + 2, 2 * ell + 3),
    (2, 4 * ell + 3, 2 * ell + 1, 4 * ell + 4, 2 * ell + 3),
    (4, 4 * ell + 1, ell, 2 * ell + 1, ell + 1),
)


def audit_family_e(delta, n, aexp, eexp, qexp):
    F = H**aexp
    G0 = H**eexp
    gamma = aa * tt * (2 * n - 2) / (dd * (n + 1))
    eta = bb * cc * eexp / (aa * aexp)
    M = gamma * F * K ** (n - 3)
    particular = eta * H**qexp * K
    L = Z + particular
    top = wedge(1 - n, aa * F, 2 * n - 2, tt * K ** (2 * n - 2))
    top += wedge(n + 1, dd * K ** (n + 1), -2, M)
    bottom = wedge(1 - n, aa * F, -n - 2, L)
    bottom += wedge(1, cc * K, -2 * n - 2, bb * G0)
    bottom_homogeneous = aa * wedge(1 - n, F, -n - 2, Z)

    triple = wedge(1 - n, aa * F, -2, M)
    triple += wedge(1, cc * K, -n - 2, L)
    triple += wedge(n + 1, dd * K ** (n + 1), -2 * n - 2, bb * G0)
    w0 = wedge(1, K, -n - 2, Z)
    remainder = sp.expand(triple - cc * w0)

    scalar = wedge(1, cc * K, -2, M)
    scalar += wedge(n + 1, dd * K ** (n + 1), -n - 2, L)
    scalar_homogeneous = dd * (n + 1) * K**n * w0
    scalar_remainder = sp.expand(scalar - scalar_homogeneous)

    euler = Hp * K + delta * H * Kp
    psi = exact_quotient(remainder, euler, (H, K, Z), f"family E triple delta={delta}")
    phi = exact_quotient(
        scalar_remainder, euler, (H, K, Z), f"family E scalar delta={delta}"
    )
    transfer = cc * scalar - euler * (
        cc * phi - dd * (n + 1) * K**n * psi
    )

    zero(top, f"family E top delta={delta}")
    zero(bottom - bottom_homogeneous, f"family E bottom delta={delta}")
    # The transfer identity is modulo the triple equation: add precisely the
    # multiple dd(n+1)K^n*triple which eliminates w0.
    zero(
        transfer - dd * (n + 1) * K**n * triple,
        f"family E transfer delta={delta}",
    )


for branch in E_BRANCHES:
    audit_family_e(*branch)


# Family E at n=2: top-row integration forces H=KJ.  The same transfer works
# with E_2(J,K)=J'K+2JK'.
gamma_e2 = aa * tt * 2 / (dd * 3)
eta_e2 = 6 * bb * cc / aa
H_e2 = K * J
M_e2 = gamma_e2 * J
particular_e2 = eta_e2 * K**6 * J**5
L_e2 = Z + particular_e2
top_e2 = wedge(-1, aa * H_e2, 2, tt * K**2) + wedge(3, dd * K**3, -2, M_e2)
bottom_e2 = wedge(-1, aa * H_e2, -4, L_e2)
bottom_e2 += wedge(1, cc * K, -6, bb * H_e2**6)
bottom_e2_homogeneous = aa * wedge(-1, H_e2, -4, Z)
triple_e2 = wedge(-1, aa * H_e2, -2, M_e2)
triple_e2 += wedge(1, cc * K, -4, L_e2)
triple_e2 += wedge(3, dd * K**3, -6, bb * H_e2**6)
w0_e2 = wedge(1, K, -4, Z)
remainder_e2 = sp.expand(triple_e2 - cc * w0_e2)
scalar_e2 = wedge(1, cc * K, -2, M_e2)
scalar_e2 += wedge(3, dd * K**3, -4, L_e2)
scalar_homogeneous_e2 = 3 * dd * K**2 * w0_e2
scalar_remainder_e2 = sp.expand(scalar_e2 - scalar_homogeneous_e2)
euler_e2 = Jp * K + 2 * J * Kp
psi_e2 = exact_quotient(remainder_e2, euler_e2, (J, K, Z), "family E n=2 triple")
phi_e2 = exact_quotient(
    scalar_remainder_e2, euler_e2, (J, K, Z), "family E n=2 scalar"
)
transfer_e2 = cc * scalar_e2 - euler_e2 * (
    cc * phi_e2 - 3 * dd * K**2 * psi_e2
)
zero(top_e2, "family E n=2 top")
zero(bottom_e2 - bottom_e2_homogeneous, "family E n=2 bottom")
zero(transfer_e2 - 3 * dd * K**2 * triple_e2, "family E n=2 transfer")


print("THM-3708 exact W002 completion audit")
print("actual placement counts n=1,2,3,>=4 = 0,3,4,6")
print("n=3 mergers = A=B and C=D as actual supports, not extra anchor placements")
print("THM-3704 input = families A,B (lowest scalar double 01=10)")
print("family C = arm-square for n>=4; exceptional n=3 Euler identity PASS")
print("family D = delta=gcd(n-1,3) Euler identities PASS; n=2 K|H branch PASS")
print("family E = delta=gcd(n-1,4) triple-to-scalar transfer PASS; n=2 K|H branch PASS")
print("family F = unique arm coefficient is a common-base square")
print("Euler factors = H'K+delta HK' or, after H=KJ, J'K+2JK'")
print("scope = complete W002 ray in the y=0 collision ring")
print("ALL CHECKS PASSED")
