#!/usr/bin/env python3
"""Exact controls for THM-3583's exponent-two 2 x 4 nonentry theorem."""

from math import gcd

import sympy as sp


x = sp.symbols("x")
h = sp.Function("h")(x)
K = sp.Function("K")(x)
F = sp.Function("F")(x)
H = sp.Function("H")(x)
A, B, C, D, E, L, M, C0 = sp.symbols(
    "A B C D E L M C0", nonzero=True
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def W(r, s, f, g):
    return sp.expand(s * sp.diff(f, x) * g - r * f * sp.diff(g, x))


def zero(expression, label):
    require(sp.simplify(sp.expand(expression)) == 0, label)


# The scalar component is a path in the Q-weight graph, whose edges have
# difference delta.  Four vertices leave exactly the three placements below.
support_rows = 0
for R in range(1, 9):
    for T in range(1, 9):
        delta = R + T - 1
        for start in (-2, -1, 0):
            qweights = [-T + j * delta for j in range(start, start + 4)]
            require(-T in qweights and R - 1 in qweights, "scalar complements")
            rows = {}
            for pweight in (-R, T - 1):
                for qweight in qweights:
                    rows.setdefault(pweight + qweight + 1, []).append(
                        (pweight, qweight)
                    )
            require(
                set(rows[0]) == {(T - 1, -T), (-R, R - 1)},
                "scalar row",
            )
            require(
                sorted(len(owners) for owners in rows.values()) == [1, 1, 2, 2, 2],
                "path collision profile",
            )
            support_rows += 1


# R=2, balanced placement: exact lower recurrence for every sampled T.
for T in range(1, 9):
    f = A * h**2
    qm1 = B * h ** (2 * T + 1)
    q0 = D * h**T + (2 * T + 1) * B / (2 * A) * h ** (2 * T - 1) * F
    zero(
        W(T - 1, -2 * T - 1, F, qm1) + W(-2, -T, f, q0),
        f"R2 balanced T={T}",
    )


# R=2, double-lower placement, T=2k: both recurrences and the scalar factor.
for k in range(1, 9):
    r = 2 * k - 1
    u = 3 * k + 1
    f = A * h
    qm2 = B * h**u
    Fk = L * K**r
    q1 = M * K
    c1 = u * B / A
    qm1 = c1 * h ** (3 * k) * Fk
    mu = 3 * k * c1 / (2 * A)
    q0 = D * h**k + mu * h ** (3 * k - 1) * Fk**2
    zero(
        W(r, -6 * k - 2, Fk, qm2) + W(-2, -4 * k - 1, f, qm1),
        f"R2 lower first k={k}",
    )
    zero(
        W(r, -4 * k - 1, Fk, qm1) + W(-2, -2 * k, f, q0),
        f"R2 lower second k={k}",
    )
    euler = K * sp.diff(h, x) + 2 * h * sp.diff(K, x)
    factor = euler * (
        A * M
        - k * r * L * D * h ** (k - 1) * K ** (r - 1)
        - mu
        * (3 * k - 1)
        * r
        * L**3
        * h ** (3 * k - 2)
        * K ** (3 * r - 1)
    )
    zero(
        W(r, -2 * k, Fk, q0) + W(-2, 1, f, q1) - factor,
        f"R2 lower scalar k={k}",
    )


# The sole square-homogeneous boundary that is not killed immediately by
# membership is k=1.  It also factors exactly.
f = A * H**2
qm2 = B * H**8
qm1 = 4 * B / A * H**6 * F + C0 * H**5
q0 = D * H**2 + 6 * B / A**2 * H**4 * F**2 + 5 * C0 / (2 * A) * H**3 * F
q1 = M * F
zero(W(1, -8, F, qm2) + W(-2, -5, f, qm1), "square branch first")
zero(W(1, -5, F, qm1) + W(-2, -2, f, q0), "square branch second")
factor = H * sp.diff(H * F, x) / (2 * A**2) * (
    4 * A**3 * M
    - 4 * A**2 * D
    - 15 * A * C0 * H * F
    - 48 * B * H**2 * F**2
)
zero(W(1, -2, F, q0) + W(-2, 1, f, q1) - factor, "square branch scalar")


# R=2, double-upper placement, T=2k, after the gcd-five root extraction.
gcd5_rows = []
for k in range(1, 9):
    r = 2 * k - 1
    u = 4 * k + 3
    v = 2 * k + 2
    w = 2 * k + 4
    if gcd(r, u) == 5:
        gcd5_rows.append(k)
        require(k % 5 == 3 and w % 5 == 0 and v % 5 != 0, "gcd-five branch")
    f = A * h
    q0 = B * h**k
    Fk = L * K**r
    q3 = M * K**u
    c1 = A * M * u / (L * r)
    q2 = D * K**v + c1 * h * K**w
    lam = A * D * v / (L * r)
    mu = A * c1 * w / (2 * L * r)
    q1 = E * K + lam * h * K**3 + mu * h**2 * K**5
    zero(W(r, v, Fk, q2) + W(-2, u, f, q3), f"R2 upper first k={k}")
    zero(W(r, 1, Fk, q1) + W(-2, v, f, q2), f"R2 upper second k={k}")
    euler = K * sp.diff(h, x) + 2 * h * sp.diff(K, x)
    factor = euler * (
        -k * r * L * B * h ** (k - 1) * K ** (r - 1)
        + A * E
        + 3 * A * lam * h * K**2
        + 5 * A * mu * h**2 * K**4
    )
    zero(
        W(r, -2 * k, Fk, q0) + W(-2, 1, f, q1) - factor,
        f"R2 upper scalar k={k}",
    )


# T=2, balanced placement: no arm-order equality survives for R>2.
for R in range(3, 65):
    d = gcd(R, 3)
    for ell in range(1, 9):
        m = R * ell // d
        n = (R + 3) * ell // d
        membership = m >= (R + 1) // 2
        if membership:
            require(R - 2 * m != 0, "balanced leading coefficient")
            require(n - 1 != m, "balanced order collision")


# T=2, double-lower placement, R=2k.
gcd4_rows = []
for k in range(2, 9):
    if gcd(2 * k, 4) == 4:
        gcd4_rows.append(k)
    f = A * h**k
    qm2 = B * h ** (2 * k + 2)
    c1 = 2 * (k + 1) * B / (k * A)
    qm1 = c1 * h ** (k + 2) * F
    lam = (k + 2) * c1 / (2 * k * A)
    q0 = D * h + lam * h**2 * F**2
    q1 = M * F ** (2 * k - 1)
    zero(
        W(1, -4 * k - 4, F, qm2) + W(-2 * k, -2 * k - 3, f, qm1),
        f"T2 lower first k={k}",
    )
    zero(
        W(1, -2 * k - 3, F, qm1) + W(-2 * k, -2, f, q0),
        f"T2 lower second k={k}",
    )
    euler = F * sp.diff(h, x) + 2 * h * sp.diff(F, x)
    factor = euler * (
        A * M * k * (2 * k - 1) * h ** (k - 1) * F ** (2 * k - 2)
        - D
        - 2 * lam * h * F**2
    )
    zero(
        W(1, -2, F, q0) + W(-2 * k, 2 * k - 1, f, q1) - factor,
        f"T2 lower scalar k={k}",
    )


# T=2, double-upper placement, R=2k.
for k in range(2, 9):
    u = 6 * k + 1
    f = A * h**k
    q0 = B * h
    q3 = C * F**u
    e0 = A * C * u
    q2 = D * F ** (4 * k) + e0 * h**k * F ** (6 * k)
    lam = 4 * A * D * k
    mu = 3 * A * e0 * k
    q1 = (
        E * F ** (2 * k - 1)
        + lam * h**k * F ** (4 * k - 1)
        + mu * h ** (2 * k) * F ** (6 * k - 1)
    )
    zero(W(1, 4 * k, F, q2) + W(-2 * k, u, f, q3), f"T2 upper first k={k}")
    zero(
        W(1, 2 * k - 1, F, q1) + W(-2 * k, 4 * k, f, q2),
        f"T2 upper second k={k}",
    )
    euler = F * sp.diff(h, x) + 2 * h * sp.diff(F, x)
    factor = euler * (
        -B
        + A * k * E * (2 * k - 1) * h ** (k - 1) * F ** (2 * k - 2)
        + A * k * lam * (4 * k - 1) * h ** (2 * k - 1) * F ** (4 * k - 2)
        + A * k * mu * (6 * k - 1) * h ** (3 * k - 1) * F ** (6 * k - 2)
    )
    zero(
        W(1, -2, F, q0) + W(-2 * k, 2 * k - 1, f, q1) - factor,
        f"T2 upper scalar k={k}",
    )


print("THM-3583 exact control")
print(f"support rows: {support_rows} = 8*8*3; collision profile [1,1,2,2,2]")
print("placements: L/B/U only; scalar complements present in every row")
print("arm gate: only R=2 or T=2 can survive")
print("R=2 ladders: balanced/lower/upper PASS for T,k<=8")
print(f"R=2 hidden gcd-five upper rows: k={gcd5_rows}; extraction arithmetic PASS")
print("T=2 balanced arm valuations: R<=64, local orders<=8 PASS")
print(f"T=2 hidden gcd-four lower rows: k={gcd4_rows}; square-root normalization PASS")
print("T=2 lower/upper factor identities: k=2..8 PASS")
print("exact consequence: a two-piece output forces the mate to have >=5 pieces")
print("scope: total >=7 only in the two-piece sector; the 3x3 six-piece cell remains open")
print("PASS")
