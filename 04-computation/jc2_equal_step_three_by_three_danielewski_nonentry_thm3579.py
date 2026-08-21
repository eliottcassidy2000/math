#!/usr/bin/env python3
"""Exact companion for THM-3579's equal-step 3x3 Darboux nonentry."""

from __future__ import annotations

from math import gcd

import sympy as s


def require(condition: bool, message: str) -> None:
    """Optimization-safe truth gate."""
    if not condition:
        raise RuntimeError(message)


def zero(expr: s.Expr, message: str) -> None:
    """Require an exact symbolic identity."""
    require(s.simplify(s.expand(expr)) == 0, message)


b, t, z = s.symbols("b t z")
h = s.Function("h")(b)
K = s.Function("K")(b)
a = s.Function("a")(b)
q = s.Function("q")(b)
w = s.Function("w")(b)
A, B, L, M = s.symbols("A B L M", nonzero=True)


def W(r: int, f: s.Expr, u: int, g: s.Expr) -> s.Expr:
    """Coefficient of {c^r f(b),c^u g(b)} after removing c^(r+u+1)."""
    return s.expand(u * s.diff(f, b) * g - r * f * s.diff(g, b))


def E(H: s.Expr, v: s.Expr) -> s.Expr:
    """The degree-rigid operator H'v+2Hv'."""
    return s.expand(s.diff(H, b) * v + 2 * H * s.diff(v, b))


def J(H: s.Expr, v: s.Expr) -> s.Expr:
    """The reflected degree-rigid operator Hv'+2H'v."""
    return s.expand(H * s.diff(v, b) + 2 * s.diff(H, b) * v)


def C13(H: s.Expr, v: s.Expr) -> s.Expr:
    """The central companion operator H'v+3Hv'."""
    return s.expand(s.diff(H, b) * v + 3 * H * s.diff(v, b))


def D32(H: s.Expr, v: s.Expr) -> s.Expr:
    """The central bridge operator 3Hv'-2H'v."""
    return s.expand(3 * H * s.diff(v, b) - 2 * s.diff(H, b) * v)


# Laurent-chart compiler.  If P=c^r F(b,c^d), Q=c^s G(b,c^d), then
# {P,Q}=c^(r+s+1) times the expression below.  The coefficient rows are the
# five additive convolution rows of two length-three progressions.
r0, s0, d0 = s.symbols("r0 s0 d0")
f = [s.Function(f"f{i}")(b) for i in range(3)]
g = [s.Function(f"g{i}")(b) for i in range(3)]
Fpoly = sum(f[i] * t**i for i in range(3))
Gpoly = sum(g[j] * t**j for j in range(3))
compiler = s.expand(
    s0 * s.diff(Fpoly, b) * Gpoly
    - r0 * Fpoly * s.diff(Gpoly, b)
    + d0 * t * (
        s.diff(Fpoly, b) * s.diff(Gpoly, t)
        - s.diff(Fpoly, t) * s.diff(Gpoly, b)
    )
)
rows = s.expand(
    sum(
        (
            (s0 + j * d0) * s.diff(f[i], b) * g[j]
            - (r0 + i * d0) * f[i] * s.diff(g[j], b)
        )
        * t ** (i + j)
        for i in range(3)
        for j in range(3)
    )
)
zero(compiler - rows, "Laurent-chart compiler")
require(
    [sum(1 for i in range(3) for j in range(3) if i + j == k) for k in range(5)]
    == [1, 2, 3, 2, 1],
    "convolution multiplicities",
)


# Squarefree arm controls and the sharp degree-two threshold used by every
# leading-degree obstruction.  The proof itself uses an arbitrary Sigma.
squarefree_degrees = []
for degree in range(2, 8):
    Sigma = s.prod(b - root for root in range(degree))
    require(s.degree(Sigma, b) == degree, "squarefree degree")
    require(s.gcd(Sigma, s.diff(Sigma, b)) == 1, "squarefree gcd")
    for root in range(degree):
        require(Sigma.subs(b, root) == 0, "sample arm missing")
        require(s.diff(Sigma, b).subs(b, root) != 0, "sample arm not simple")
    squarefree_degrees.append(degree)


# A complementary pair has weights (-N,N-1).  At a simple arm, a unit can
# survive only for N=2, negative coefficient order one, positive coefficient
# order zero.  N=1 is the essential zero-multiplier boundary.
arm_survivors = []
for N in range(1, 17):
    for negative_order in range((N + 1) // 2, (N + 1) // 2 + 5):
        for positive_order in range(4):
            order = negative_order + positive_order - 1
            multiplier = (N - 1) * negative_order + N * positive_order
            if order == 0 and multiplier != 0:
                arm_survivors.append((N, negative_order, positive_order))
require(arm_survivors == [(2, 1, 0)], "simple-arm scalar gate")


# The operator H'v+2Hv' cannot equal one if deg(H)>=2.  This exact range is a
# hostile/positive control for the all-degree leading-coefficient proof.
degree_controls = []
for degree_H in range(2, 8):
    for degree_v in range(8):
        H = b**degree_H + 2 * b + 3
        v = 3 if degree_v == 0 else 2 * b**degree_v + b + 1
        Ev = s.Poly(E(H, v), b)
        expected_degree = degree_H + degree_v - 1
        expected_lc = (
            degree_H + 2 * degree_v
        ) * s.Poly(H, b).LC() * s.Poly(v, b).LC()
        require(Ev.degree() == expected_degree, "degree-rigid operator degree")
        require(Ev.LC() == expected_lc, "degree-rigid operator leading coefficient")
    degree_controls.append(degree_H)


# Scalar row kappa=1.  The arm and lowest-extreme equations leave, after
# swapping, supports {-2,2k-1,4k} and {-2k,1,2k+2}, with f0=A*h,
# g0=B*h^k and Sigma|h.  The scalar row is exactly E_h(u)=1.
row1_controls = []
for k in range(1, 13):
    f0 = A * h
    g0 = B * h**k
    scalar = W(-2, f0, 1, q) + W(2 * k - 1, a, -2 * k, g0)
    u = A * q - k * B * h ** (k - 1) * a
    zero(scalar - E(h, u), f"row-1 scalar identity k={k}")
    row1_controls.append(k)


# Scalar row kappa=3.  After swapping, the supports are
# {-n-4,-2,n} and {-2n-3,-n-1,1}.  The top extreme gives
# F=L*K^n,G=M*K, and the scalar row is E_K(v)=1.  Both middle coefficients
# have negative weight, so Sigma divides v.
row3_controls = []
for n in range(1, 13):
    Ftop = L * K**n
    Gtop = M * K
    scalar = W(-2, a, 1, Gtop) + W(n, Ftop, -n - 1, q)
    v = M * a - n * L * K ** (n - 1) * q
    zero(scalar - J(K, v), f"row-3 scalar identity n={n}")

    Sigma = b * (b - 1) * (b + 2)
    aa = Sigma * (b + n + 1)
    qq = Sigma ** ((n + 2) // 2) * (b**2 + n + 1)
    vv = s.expand(M * aa - n * L * (b + 3) ** (n - 1) * qq)
    require(s.rem(vv, Sigma, domain=s.QQ.frac_field(L, M)) == 0,
            f"row-3 arm divisibility n={n}")
    row3_controls.append(n)


# At the n=0 boundary, the high extreme has weights (0,1), so its zero
# Wronskian is F'G=0 and the weight-zero coefficient is removable.
F0 = s.Function("F0")(b)
G1 = s.Function("G1")(b)
zero(W(0, F0, 1, G1) - s.diff(F0, b) * G1,
     "row-3 n=0 removable boundary")


# Central row support census.  With low weights -R,-T, the step is
# d=(R+T-1)/2 and the middle P weight is p=(T-R-1)/2.  Endpoint scalar
# channels are killed by parity plus the low common-power exponents; only
# p=1 or p=-2 survives, i.e. |R-T|=3.
central_profiles = []
for R in range(1, 33):
    for T in range(1, 33):
        if (R + T) % 2 == 0:
            continue
        d = (R + T - 1) // 2
        p = (T - R - 1) // 2
        delta = gcd(R, T)
        endpoint_R = R == 2 and R // delta == 1
        endpoint_T = T == 2 and T // delta == 1
        middle = p in (1, -2)
        if endpoint_R or endpoint_T or middle:
            central_profiles.append((R, T, d, p))
require(
    all(abs(R - T) == 3 for R, T, _, _ in central_profiles),
    "central support classification",
)


# In the orientation T=R+3, the lower bridge can cancel at an arm only when
# gcd(R,3)=3 and the common base h is simple there.  This finite census is a
# control for the exact valuation equation (3/gcd(R,3))*ord(h)=1.
central_arm_rows = []
for R in range(1, 49):
    delta = gcd(R, 3)
    alpha = R // delta
    beta = (R + 3) // delta
    for arm_order in range(1, 9):
        if (beta - alpha) * arm_order == 1:
            central_arm_rows.append((R, delta, arm_order))
require(
    central_arm_rows == [(R, 3, 1) for R in range(3, 49, 3)],
    "central arm valuation gate",
)


# Exact central bridge identities.  Put R=3m,T=3m+3.  The lower bridge
# forces q=lambda*h*a; the upper bridge becomes a derivative and integrates
# to a(A0-A1*h*K^3)=C*K.
central_controls = []
for m in range(1, 9):
    R = 3 * m
    lam = s.Rational(m + 1, m) * B / A
    f0 = A * h**m
    g0 = B * h ** (m + 1)
    lower = W(-R, f0, -2, q) + W(1, a, -R - 3, g0)
    expected_lower = h ** (m - 1) * (
        m * A * D32(h, q) - (m + 1) * B * h * C13(h, a)
    )
    zero(lower - expected_lower, f"central lower bridge m={m}")
    zero(D32(h, h * a) - h * C13(h, a),
         f"central particular solution m={m}")

    q_particular = lam * h * a
    lower_particular = (
        W(-R, f0, -2, q_particular) + W(1, a, -R - 3, g0)
    )
    zero(lower_particular, f"central forced lower solution m={m}")

    Ftop = L * K ** (R + 2)
    Gtop = M * K ** (R - 1)
    upper = W(1, a, R - 1, Gtop) + W(R + 2, Ftop, -2, q_particular)
    derivative_form = K**R * s.diff(
        M * (R - 1) * a / K - lam * L * (R + 2) * h * a * K**2,
        b,
    )
    zero(upper - derivative_form, f"central upper derivative m={m}")
    central_controls.append(m)


# Rational first integral for the homogeneous central bridge.
zero(
    s.diff(w**3 / h**2, b) - w**2 * D32(h, w) / h**3,
    "central cube first integral",
)


# Coprimality/nonconstancy controls for the terminal equation.  Modulo K,
# A0-A1*h*K^3 is the nonzero constant A0.
coprime_controls = []
for m in range(1, 9):
    hs = b * (b - 1) * (b + 2)
    Ks = b**2 + b + m + 1
    terminal = (3 * m - 1) - (3 * m + 2) * hs * Ks**3
    require(s.degree(terminal, b) > 0, "terminal factor nonconstant")
    require(s.gcd(s.Poly(Ks, b), s.Poly(terminal, b)).degree() == 0,
            "terminal factor coprime to K")
    coprime_controls.append(m)


# Sharp hostiles and near misses.
# 1. The homogeneous differential equation has a cube-resonant solution if
#    one drops the simple squarefree arm.
zero(3 * z**3 * s.diff(z**2, z) - 2 * s.diff(z**3, z) * z**2,
     "cube-resonance hostile")

# 2. Degree one defeats the E_H degree gate; the corresponding target is A2.
zero(s.diff(b, b) * 1 + 2 * b * s.diff(s.Integer(1), b) - 1,
     "degree-one hostile")

# 3. The rational pole pair is an actual Darboux pair outside A_Sigma.
zero(W(0, b, -1, s.Integer(-1)) - 1, "rational-pole hostile")

# 4. THM-3572's selected symmetric central profile {-3,-1,1} against
#    {-2,0,2} already has scalar row divisible by the arm parameter once the
#    two extreme common-power equations are imposed.
Az, Bz, Lz, Mz = s.symbols("Az Bz Lz Mz", nonzero=True)
f_local = Az * z**3
g_local = Bz * z**2
a_local = z
q_local = s.Integer(1)
F_local = Lz
G_local = Mz


def Wz(r: int, fz: s.Expr, u: int, gz: s.Expr) -> s.Expr:
    return s.expand(u * s.diff(fz, z) * gz - r * fz * s.diff(gz, z))


selected_scalar = s.factor(
    Wz(-3, f_local, 2, G_local)
    + Wz(-1, a_local, 0, q_local)
    + Wz(1, F_local, -2, g_local)
)
require(s.rem(selected_scalar, z, domain=s.QQ.frac_field(Az, Bz, Lz, Mz)) == 0,
        "selected symmetric profile arm factor")


print("THM-3579 equal-step three-by-three Danielewski Darboux nonentry")
print("Laurent compiler and convolution multiplicities [1,2,3,2,1]: PASS")
print(f"squarefree degree controls: {squarefree_degrees}")
print(f"simple-arm scalar survivor: {arm_survivors}")
print(f"degree-rigid operator controls: deg(H)={degree_controls}, deg(v)=0..7")
print(f"row-1 exact ladder controls: k={row1_controls}")
print(f"row-3 exact ladder controls: n={row3_controls}")
print("central post-extreme support law: |R-T|=3")
print(f"central exact bridge controls: m={central_controls}")
print(f"terminal coprimality controls: m={coprime_controls}")
print("sharp hostiles: cube resonance / degree one / rational pole / selected profile")
print("all optimization-safe exact truth gates passed")
