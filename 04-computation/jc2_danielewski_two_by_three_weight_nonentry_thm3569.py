#!/usr/bin/env python3
"""Exact companion for THM-3569's two-by-three weight nonentry."""

from __future__ import annotations

import sympy as s


def require(condition: bool, message: str) -> None:
    """Optimization-safe truth gate."""
    if not condition:
        raise RuntimeError(message)


b = s.symbols("b")
h = s.Function("h")(b)
K = s.Function("K")(b)
J = s.Function("J")(b)
a, B, L, M, D = s.symbols("a B L M D", nonzero=True)


def W(r: int, f: s.Expr, t: int, g: s.Expr) -> s.Expr:
    """Coefficient of {c^r f(b),c^t g(b)} after removing c^(r+t+1)."""
    return s.expand(t * s.diff(f, b) * g - r * f * s.diff(g, b))


def zero(expr: s.Expr, message: str) -> None:
    require(s.simplify(s.expand(expr)) == 0, message)


# Support arithmetic.  With P weights r,r+delta, the complementary Q weights
# are aa=-r-delta-1 and bb=-r-1.  A third Q weight x can collide with a
# nonzero cross row only at aa-delta or bb+delta (the choices aa,bb are not
# extra weights).  Everything else isolates the third piece and reduces to
# the already excluded two-by-two cell.
for r in range(-12, 13):
    for delta in range(1, 13):
        aa = -r - delta - 1
        bb = -r - 1
        for x in range(aa - 2 * delta - 2, bb + 2 * delta + 3):
            if x in (aa, bb):
                continue
            extra_rows = {x - bb, x - aa}
            collisions = extra_rows.intersection({-delta, delta})
            expected = (
                {-delta} if x == aa - delta else
                {delta} if x == bb + delta else
                set()
            )
            require(collisions == expected, "support-collision classification failed")


# Lower arithmetic extension.  Arm orders reduce it to R=T=2, hence
# supp(P)={-2,1}, supp(Q)={-5,-2,1}.  The middle row integrates explicitly,
# and the constant row has a nonconstant h factor.
lam = 5 * L * B / (2 * a)
f = a * h**2
g0 = B * h**5
F = L * K
H = M * K
g1 = h**2 * (D + lam * h * K)
zero(W(-2, f, -2, g1) + W(1, F, -5, g0), "lower middle row")
constant_row = W(-2, f, 1, H) + W(1, F, -2, g1)
lower_factor = h * s.diff(h * K, b) * (2 * a * M - L * (2 * D + 3 * lam * h * K))
zero(constant_row - lower_factor, "lower constant-row factorization")


# Upper arithmetic extension A: R=2,T=2k.  The common-power relation at the
# high extreme has d=gcd(2k-1,2k+2) in {1,3}.  Both the d=1 branch and both
# d=3 homogeneous branches factor through a positive-degree derivative.
upper_a_rows = []
for k in range(1, 9):
    p, q = 2 * k - 1, 2 * k + 2
    d = int(s.gcd(p, q))
    require(d in (1, 3), "unexpected upper-A gcd")
    m, n = p // d, q // d
    f = a * h
    g = B * h**k
    F = L * K**m
    H = M * K**n
    C = a * M * n / (L * m)
    if d == 1:
        G = D * K + C * h * K**3
        deriv = s.diff(h, b) * K + 2 * h * s.diff(K, b)
        zero(W(-2, f, q, H) + W(p, F, 1, G), f"upper A d=1 middle k={k}")
        expected = deriv * (
            a * D + 3 * a * C * h * K**2
            - B * L * p * k * h ** (k - 1) * K ** (p - 1)
        )
        zero(W(-2, f, 1, G) + W(p, F, -2 * k, g) - expected,
             f"upper A d=1 constant k={k}")
        upper_a_rows.append((k, d, "ordinary"))
    else:
        G = C * h * K
        deriv3 = 3 * s.diff(h, b) * K + 2 * h * s.diff(K, b)
        zero(W(-2, f, q, H) + W(p, F, 1, G), f"upper A d=3 middle k={k}")
        expected = deriv3 * (
            a * C * h - B * L * m * k * h ** (k - 1) * K ** (m - 1)
        )
        zero(W(-2, f, 1, G) + W(p, F, -2 * k, g) - expected,
             f"upper A d=3 particular k={k}")

        # A nonzero homogeneous solution of K'G0-3KG0'=0 has
        # K=lambda*J^3 and G0=mu*J.  Scalars are absorbed into J and D.
        FJ = L * J**p
        HJ = M * J**q
        GJ = C * h * J**3 + D * J
        deriv_j = s.diff(h, b) * J + 2 * h * s.diff(J, b)
        zero(W(-2, f, q, HJ) + W(p, FJ, 1, GJ),
             f"upper A d=3 homogeneous middle k={k}")
        expected_j = deriv_j * (
            a * D + 3 * a * C * h * J**2
            - 3 * B * L * m * k * h ** (k - 1) * J ** (3 * m - 1)
        )
        zero(W(-2, f, 1, GJ) + W(p, FJ, -2 * k, g) - expected_j,
             f"upper A d=3 homogeneous constant k={k}")
        upper_a_rows.append((k, d, "particular+cube"))


# Upper arithmetic extension B: T=2,R=2k.  The middle row fixes C and the
# constant row always contains h'K+2hK'.
upper_b_rows = []
for k in range(1, 9):
    R = 2 * k
    f = a * h**k
    g = B * h
    F = L * K
    H = M * K ** (4 * k)
    C = 4 * k * a * M / L
    G = D * K ** (2 * k - 1) + C * h**k * K ** (4 * k - 1)
    deriv = s.diff(h, b) * K + 2 * h * s.diff(K, b)
    zero(W(-R, f, 4 * k, H) + W(1, F, R - 1, G),
         f"upper B middle k={k}")
    expected = deriv * (
        a * D * k * (2 * k - 1) * h ** (k - 1) * K ** (2 * k - 2)
        + a * C * k * (4 * k - 1) * h ** (2 * k - 1) * K ** (4 * k - 2)
        - B * L
    )
    zero(W(-R, f, R - 1, G) + W(1, F, -2, g) - expected,
         f"upper B constant k={k}")
    upper_b_rows.append(k)


print("THM-3569 Danielewski two-by-three weight nonentry")
print("support collision classes: isolated / lower arithmetic / upper arithmetic")
print("lower exceptional ladder {-2,1} x {-5,-2,1}: exact factor PASS")
print(f"upper A exact controls: {upper_a_rows}")
print(f"upper B exact controls: {upper_b_rows}")
print("every surviving constant row has a positive-degree differential factor")
print("all optimization-safe truth gates passed")
