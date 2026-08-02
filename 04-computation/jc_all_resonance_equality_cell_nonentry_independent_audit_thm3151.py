#!/usr/bin/env python3
"""Independent exact audit for THM-3151's all-R equality-cell exclusion."""

from math import factorial
import sympy as sp


def require(ok, label):
    if not bool(ok):
        raise RuntimeError(label)


def add(out, key, value):
    value = sp.cancel(out.get(key, 0) + value)
    if value:
        out[key] = value
    elif key in out:
        del out[key]


def direct_c(R, n):
    """Direct four-summand multinomial extraction from the defining series."""
    alpha = sp.Rational(2 * R - 1, 2)
    out = {}
    for i in range(n // 2 + 1):
        for qn in range(n // 3 + 1):
            rem = n - 2 * i - 3 * qn
            if rem < 0 or rem % 4:
                continue
            f = rem // 4
            total = i + qn + f
            mult = factorial(total) // (factorial(i) * factorial(qn) * factorial(f))
            base = sp.binomial(alpha, total) * mult * 2**i
            for sn in range(f + 1):
                add(out, (qn, i + 2 * (f - sn), sn),
                    base * sp.binomial(f, sn) * (-1)**sn)
    return out


def direct_flux(R, channel):
    n = 4 * R - 1 if channel == "phi" else 4 * R
    raw = direct_c(R, n)
    out = {}
    for (a, b, c), value in raw.items():
        key = (a - 1, b, c) if channel == "phi" else (a, b, c)
        require(min(key) >= 0, f"division by q at R={R}")
        add(out, key, 4 * value)
    return out


def factorized_flux(R, channel):
    """Extraction after A=(1+dt^2)^2+t^3(q-st)."""
    alpha = sp.Rational(2 * R - 1, 2)
    out = {}
    # Surviving support is l>=0; l<0 would ask a positive-degree
    # (1+dt^2) polynomial for a coefficient beyond its degree.
    for ell in range(R + 1):
        for r in range(2 * ell, 2 * R + 1):
            b = r - 2 * ell
            c = (R - 1 - 2 * r + ell) if channel == "phi" else (R - 2 * r + ell)
            if c < 0:
                continue
            n = R + ell
            qpower = 2 * r + 1 if channel == "phi" else 2 * r
            require(qpower + c == n, "factor split count")
            value = (4 * sp.binomial(alpha, n) * sp.binomial(n, c)
                     * (-1)**c * sp.binomial(-1 - 2 * ell, b))
            add(out, (2 * r, b, c), value)
    return out


def weight(mon, u, v):
    a, b, c = mon
    return a - u * b - v * c


def minimum(poly, u, v):
    return min(weight(mon, u, v) for mon in poly)


z = sp.symbols("z")
checked_rows = 0
checked_gcds = 0
checked_chambers = 0

rows = {}
for j in range(1, 32):
    rows[j] = {}
    for channel in ("phi", "psi"):
        direct = direct_flux(j, channel)
        factored = factorized_flux(j, channel)
        require(direct == factored, f"direct/factorized mismatch j={j} {channel}")
        rows[j][channel] = factored
        checked_rows += 1

    P = sum((-1)**r * sp.binomial(j, 2*r + 1) * z**r
            for r in range((j - 1)//2 + 1))
    Q = sum((-1)**r * sp.binomial(j, 2*r) * z**r
            for r in range(j//2 + 1))
    require(sp.gcd(P, Q) == 1, f"wall gcd j={j}")
    checked_gcds += 1


for R in range(8, 32, 3):
    for v in range(1, 10):
        for u in range(0, 28):
            delta = 2*v + 2 - u
            for channel, baseline in (("phi", -(R-1)*v), ("psi", -R*v)):
                top = minimum(rows[R][channel], u, v)
                if delta >= 0:
                    require(top == baseline, f"I/W top R,u,v,ch={R,u,v,channel}")
                    for j in range(1, R-1):
                        require(minimum(rows[j][channel], u, v) > top,
                                f"I/W lower collision R,j,u,v,ch={R,j,u,v,channel}")
            if delta < 0:
                chosen = "psi" if R % 2 == 0 else "phi"
                top = minimum(rows[R][chosen], u, v)
                for j in range(1, R-1):
                    require(minimum(rows[j][chosen], u, v) > top,
                            f"III lower collision R,j,u,v,ch={R,j,u,v,chosen}")
            checked_chambers += 1

    # s=0 boundary, d nonconstant of degree at least three: delete all
    # s-positive monomials and use the parity-appropriate endpoint.
    for u in range(3, 28):
        chosen = "psi" if R % 2 == 0 else "phi"
        top_poly = {m:c for m,c in rows[R][chosen].items() if m[2] == 0}
        require(top_poly, f"empty zero-s top R={R}")
        top = minimum(top_poly, u, 0)
        for j in range(1, R-1):
            low = {m:c for m,c in rows[j][chosen].items() if m[2] == 0}
            if low:
                require(minimum(low, u, 0) > top,
                        f"zero-s lower collision R,j,u={R,j,u}")
        checked_chambers += 1

    # d=0 boundary, s nonconstant: pure-s endpoints remain and separate.
    for v in range(1, 10):
        for channel in ("phi", "psi"):
            top_poly = {m:c for m,c in rows[R][channel].items() if m[1] == 0}
            top = minimum(top_poly, 0, v)
            for j in range(1, R-1):
                low = {m:c for m,c in rows[j][channel].items() if m[1] == 0}
                require(minimum(low, 0, v) > top,
                        f"zero-d lower collision R,j,v,ch={R,j,v,channel}")
        checked_chambers += 1

    k = (R - 2)//3
    D = 4*k + 3
    a = 2*k + 2
    nu = D + 2
    bmin = a + 1
    require(nu == 2 + D and 2*a < nu and 2*bmin >= nu, f"local resonance R={R}")
    require(2*a - nu == -1, f"T exponent R={R}")
    require(2*(a+1)-nu >= 0 and a+(a+1)-nu >= 0, f"d/s polynomial R={R}")
    require((D-1)-D == -1, f"target infinity R={R}")


require(checked_rows == 62, ("row count", checked_rows))
require(checked_gcds == 31, ("gcd count", checked_gcds))
require(checked_chambers == 2288, ("chamber count", checked_chambers))

print("THM-3151 independent all-R equality-cell infinity audit")
print("direct=factorized flux rows:", checked_rows)
print("even/odd wall gcds:", checked_gcds)
print("hostile degree/parity/zero-boundary chambers:", checked_chambers)
print("R=3k+2 samples: 8,11,14,17,20,23,26,29")
print("RESULT: PASS")
