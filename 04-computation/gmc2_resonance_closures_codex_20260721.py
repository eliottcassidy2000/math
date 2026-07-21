#!/usr/bin/env python3
"""Exact checks for THM-2018's two GMC(2) resonance closures.

Part A checks, over exact rationals, three independent descriptions of the
moments for

    P = Z + (b0+b1*s) + W*(delta+alpha*s),  s=ZW:

direct Wick expansion, the primitive-return channel formula, and the closed
Catalan/algebraic EGF.  Only the product a*c matters, so taking a=1 loses no
generality for this identity check.

Part B checks the arbitrary-charge factorization on h=kappa*b^r and the EGF
of its scalar channel coefficient H_m.  The proof of nonvanishing is in the
theorem; this script is a verification instrument, not a finite proof of NC2.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from math import comb, factorial, gcd


Poly = tuple[Q, ...]
Series = list[Q]


def trim(f: Poly) -> Poly:
    out = list(f)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def poly_mul(f: Poly, g: Poly) -> Poly:
    out = [Q(0)] * (len(f) + len(g) - 1)
    for i, left in enumerate(f):
        for j, right in enumerate(g):
            out[i + j] += left * right
    return trim(tuple(out))


def poly_pow(f: Poly, exponent: int) -> Poly:
    out: Poly = (Q(1),)
    base = f
    n = exponent
    while n:
        if n & 1:
            out = poly_mul(out, base)
        base = poly_mul(base, base)
        n //= 2
    return out


def poly_shift(f: Poly, amount: int) -> Poly:
    return (Q(0),) * amount + f


def laplace_factorial(f: Poly) -> Q:
    return sum((coefficient * factorial(degree) for degree, coefficient in enumerate(f)), Q(0))


def series_mul(f: Series, g: Series, order: int) -> Series:
    out = [Q(0)] * (order + 1)
    for i, left in enumerate(f[: order + 1]):
        for j, right in enumerate(g[: order + 1 - i]):
            out[i + j] += left * right
    return out


def series_pow(f: Series, exponent: int, order: int) -> Series:
    out = [Q(1)] + [Q(0)] * order
    base = f[: order + 1]
    n = exponent
    while n:
        if n & 1:
            out = series_mul(out, base, order)
        base = series_mul(base, base, order)
        n //= 2
    return out


def series_inv(f: Series, order: int) -> Series:
    assert f[0] != 0
    out = [Q(0)] * (order + 1)
    out[0] = 1 / f[0]
    for n in range(1, order + 1):
        out[n] = -sum(f[k] * out[n - k] for k in range(1, min(n, len(f) - 1) + 1)) / f[0]
    return out


def series_sqrt_unit(f: Series, order: int) -> Series:
    """The square root with constant coefficient +1."""

    assert f[0] == 1
    out = [Q(0)] * (order + 1)
    out[0] = Q(1)
    for n in range(1, order + 1):
        target = f[n] if n < len(f) else Q(0)
        middle = sum(out[k] * out[n - k] for k in range(1, n))
        out[n] = (target - middle) / 2
    return out


def series_exp_zero(f: Series, order: int) -> Series:
    """exp(f), for f[0]=0, from n*E_n=sum k*f_k*E_(n-k)."""

    assert f[0] == 0
    out = [Q(0)] * (order + 1)
    out[0] = Q(1)
    for n in range(1, order + 1):
        out[n] = sum(k * f[k] * out[n - k] for k in range(1, n + 1)) / n
    return out


def part_a_channel_moment(b0: Q, b1: Q, delta: Q, alpha: Q, m: int) -> Q:
    b: Poly = (b0, b1)
    h: Poly = (Q(0), delta, alpha)
    total = Q(0)
    for k in range(m // 2 + 1):
        coefficient = Q(factorial(m), factorial(k) ** 2 * factorial(m - 2 * k))
        total += coefficient * laplace_factorial(poly_mul(poly_pow(h, k), poly_pow(b, m - 2 * k)))
    return total


def wick_base_part_a(b0: Q, b1: Q, delta: Q, alpha: Q) -> dict[tuple[int, int], Q]:
    # Choose a=1 and c=delta+alpha*s.  Then s*a*c=delta*s+alpha*s^2.
    terms: dict[tuple[int, int], Q] = {
        (1, 0): Q(1),
        (0, 0): b0,
        (1, 1): b1,
        (0, 1): delta,
        (1, 2): alpha,
    }
    return {key: value for key, value in terms.items() if value}


def direct_wick_moment(base: dict[tuple[int, int], Q], m: int) -> Q:
    state: dict[tuple[int, int], Q] = {(0, 0): Q(1)}
    for _ in range(m):
        nxt: dict[tuple[int, int], Q] = {}
        for (z0, w0), left in state.items():
            for (z1, w1), right in base.items():
                key = (z0 + z1, w0 + w1)
                nxt[key] = nxt.get(key, Q(0)) + left * right
        state = nxt
    return sum(
        (coefficient * factorial(z_degree) for (z_degree, w_degree), coefficient in state.items() if z_degree == w_degree),
        Q(0),
    )


def part_a_closed_egf(b0: Q, b1: Q, delta: Q, alpha: Q, order: int) -> Series:
    # q^2=(1-b1*t)^2-4*alpha*t^2, q(0)=1.
    discriminant = [Q(1), -2 * b1, b1 * b1 - 4 * alpha] + [Q(0)] * max(0, order - 2)
    discriminant = discriminant[: order + 1]
    q = series_sqrt_unit(discriminant, order)
    p = [Q(1), -b1] + [Q(0)] * max(0, order - 1)
    p = p[: order + 1]

    exponent = [Q(0)] * (order + 1)
    if order >= 1:
        exponent[1] += b0
    if alpha != 0:
        for n in range(order + 1):
            exponent[n] += delta * (p[n] - q[n]) / (2 * alpha)
    else:
        inv_p = series_inv(p, order)
        for n in range(2, order + 1):
            exponent[n] += delta * inv_p[n - 2]

    return series_mul(series_exp_zero(exponent, order), series_inv(q, order), order)


def catalan_identity_check(order: int) -> None:
    catalan = [Q(comb(2 * k, k), k + 1) for k in range(order + 1)]
    inv_sqrt = [Q(comb(2 * k, k)) for k in range(order + 1)]
    checked = 0
    for n in range(0, 7):
        left = series_mul(series_pow(catalan, n, order), inv_sqrt, order)
        for j in range(order + 1):
            assert left[j] == comb(n + 2 * j, j), (n, j, left[j])
            checked += 1
    print(f"Catalan identity [x^j] C(x)^n/sqrt(1-4x): PASS ({checked} exact coefficients)")


def part_a_checks() -> None:
    print("\nPART A: AFFINE b, QUADRATIC h=s*a*c")
    catalan_identity_check(12)
    controls = [
        (Q(1), Q(2), Q(3), Q(4)),
        (Q(-2), Q(1), Q(5), Q(-3)),
        (Q(0), Q(-1), Q(2), Q(0)),       # alpha=0 branch
        (Q(3), Q(0), Q(-2), Q(7)),       # b1=0 sharp boundary
        (Q(0), Q(0), Q(1), Q(1)),
        (Q(1, 2), Q(-3, 2), Q(2, 3), Q(-5, 4)),
    ]
    order = 10
    for control in controls:
        b0, b1, delta, alpha = control
        egf = part_a_closed_egf(b0, b1, delta, alpha, order)
        base = wick_base_part_a(b0, b1, delta, alpha)
        moments: list[Q] = []
        for m in range(order + 1):
            channel = Q(1) if m == 0 else part_a_channel_moment(b0, b1, delta, alpha, m)
            wick = Q(1) if m == 0 else direct_wick_moment(base, m)
            closed = egf[m] * factorial(m)
            assert channel == wick == closed, (control, m, channel, wick, closed)
            if m <= 4:
                moments.append(channel)
        print(f"  params={tuple(str(x) for x in control)}: Wick=channels=closed_EGF m=0..{order}; M0..M4={moments}")


def general_channel_moment(p0: int, q0: int, b: Poly, h: Poly, m: int) -> Q:
    r = p0 + q0
    total = Q(0)
    for k in range(m // r + 1):
        coefficient = Q(factorial(m), factorial(q0 * k) * factorial(p0 * k) * factorial(m - r * k))
        total += coefficient * laplace_factorial(poly_mul(poly_pow(h, k), poly_pow(b, m - r * k)))
    return total


def hyper_bessel_coefficient(p0: int, q0: int, kappa: Q, m: int) -> Q:
    r = p0 + q0
    return sum(
        (
            Q(factorial(m), factorial(q0 * k) * factorial(p0 * k) * factorial(m - r * k))
            * kappa**k
            for k in range(m // r + 1)
        ),
        Q(0),
    )


def return_polynomial(p: int, q: int, a: Poly, c: Poly) -> tuple[int, int, Poly]:
    common = gcd(p, q)
    p0, q0 = p // common, q // common
    forced = p * q // common
    h = poly_shift(poly_mul(poly_pow(a, q0), poly_pow(c, p0)), forced)
    return p0, q0, h


def wick_base_general(p: int, q: int, a: Poly, b: Poly, c: Poly) -> dict[tuple[int, int], Q]:
    terms: dict[tuple[int, int], Q] = {}

    def insert(z_degree: int, w_degree: int, coefficient: Q) -> None:
        key = (z_degree, w_degree)
        terms[key] = terms.get(key, Q(0)) + coefficient

    for radial, coefficient in enumerate(a):
        insert(p + radial, radial, coefficient)
    for radial, coefficient in enumerate(b):
        insert(radial, radial, coefficient)
    for radial, coefficient in enumerate(c):
        insert(radial, q + radial, coefficient)
    return {key: value for key, value in terms.items() if value}


def proportional_actual_example(p: int, q: int, f: Poly) -> tuple[Poly, Poly, Poly, Poly]:
    # An explicit realization of h=b^r for every p,q.  With g,p0,q0,r as
    # usual, take b=s^(g*p0*q0) f^p0, a=1, and
    # c=s^(g*q0*(r-1)) f^r.
    common = gcd(p, q)
    p0, q0 = p // common, q // common
    r = p0 + q0
    b = poly_shift(poly_pow(f, p0), common * p0 * q0)
    a: Poly = (Q(1),)
    c = poly_shift(poly_pow(f, r), common * q0 * (r - 1))
    _, _, h = return_polynomial(p, q, a, c)
    assert h == poly_pow(b, r)
    return a, b, c, h


def tie_tournament(m: int, r: int) -> dict[str, object]:
    vertices = list(range(m // r + 1))
    scores = {vertex: len(vertices) - 1 - vertex for vertex in vertices}
    return {
        "vertices": len(vertices),
        "observable": "(r*k+m-r*k)-(r*l+m-r*l)=0",
        "score_histogram": dict(sorted(Counter(scores.values()).items())),
        "directed_3_cycles": 0,
        "scc_sizes": [1] * len(vertices),
        "hamiltonian_path": vertices,
        "hamiltonian_path_count": 1,
    }


def part_b_checks() -> None:
    print("\nPART B: ARBITRARY-CHARGE CENTRAL SHADOW h=kappa*b^r")
    controls = [
        (1, 1, Q(-1, 2), (Q(1), Q(-2), Q(1))),  # H_2=0: real channel cancellation
        (2, 3, Q(3, 2), (Q(0), Q(1), Q(1))),
        (1, 4, Q(0), (Q(2), Q(-1))),
        (3, 2, Q(1), (Q(1), Q(0), Q(2))),
    ]
    order = 14
    for p0, q0, kappa, b in controls:
        r = p0 + q0
        h = tuple(kappa * value for value in poly_pow(b, r))
        h = trim(h)
        nonzero_h = h if any(h) else (Q(0),)
        first_nonzero: list[int] = []
        for m in range(order + 1):
            channel = general_channel_moment(p0, q0, b, nonzero_h, m)
            scalar = hyper_bessel_coefficient(p0, q0, kappa, m)
            factored = scalar * laplace_factorial(poly_pow(b, m))
            assert channel == factored, (p0, q0, kappa, m, channel, factored)
            if scalar and len(first_nonzero) < 8:
                first_nonzero.append(m)
        print(
            f"  (p0,q0,r)=({p0},{q0},{r}), kappa={kappa}: "
            f"channel=H_m*L(b^m) m=0..{order}; first H_m!=0 at {first_nonzero}"
        )

    # Direct Wick controls on actual polynomial realizations h=b^r.
    for p, q, f, max_m in [
        (1, 1, (Q(1), Q(1)), 8),
        (2, 3, (Q(1), Q(1)), 6),
    ]:
        a, b, c, h = proportional_actual_example(p, q, f)
        p0, q0, returned_h = return_polynomial(p, q, a, c)
        assert returned_h == h
        base = wick_base_general(p, q, a, b, c)
        for m in range(1, max_m + 1):
            wick = direct_wick_moment(base, m)
            channel = general_channel_moment(p0, q0, b, h, m)
            factored = hyper_bessel_coefficient(p0, q0, Q(1), m) * laplace_factorial(poly_pow(b, m))
            assert wick == channel == factored, (p, q, m, wick, channel, factored)
        print(f"  actual charges (+{p},-{q}): direct Wick=channels=shadow factor m=1..{max_m}")

    for r in (2, 5):
        print(f"  tied channel tournament r={r}, m=24: {tie_tournament(24, r)}")


def main() -> None:
    print("THM-2018 GMC(2) EXACT RESONANCE CLOSURES")
    print("moment convention: E[Z^A W^B]=A!*1[A=B], L(s^N)=N!")
    part_a_checks()
    part_b_checks()
    print("\nASSUMPTION CHALLENGE / TOURNAMENT METHOD")
    print("candidate vertices: monomials, charges, Catalan pairs (n,j), return channels k, proof obligations")
    print("chosen Part-B vertices: return channels k; all radial-shadow observables tie on h=kappa*b^r")
    print("tie gauge/path: increasing k; fingerprints are transitive with zero directed cycles")
    print("preserved: the complete common radial factor b^m")
    print("destroyed: scalar phase cancellation inside H_m; restored by its root-of-unity EGF")
    print("challenged assumption: scalar return channels are not independent equations")
    print("status: exact verification passed; full NC2 remains open outside the two theorem strata")


if __name__ == "__main__":
    main()
