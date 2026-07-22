#!/usr/bin/env python3
"""Numerical controls for THM-2101's additive orbit-residue identity.

The theorem is exact and does not depend on this script.  These high-precision
checks compare direct Laurent constant terms with the small-root residue sum,
verify the full-root Lagrange sum, and retain a one-sided hostile boundary.
"""
import sympy as sp


u = sp.symbols("u")


def constant_terms(R, M, order):
    return [sp.expand(R**m).coeff(u, M * m) for m in range(order + 1)]


def residue_check(name, R, M, t_value=sp.Rational(1, 1000), order=9):
    phi = sp.Poly(u**M - t_value * R, u)
    derivative = sp.diff(phi.as_expr(), u)
    roots = [complex(z) for z in sp.nroots(phi, n=70, maxsteps=300)]
    roots.sort(key=abs)
    small = roots[:M]

    def weight(alpha):
        return alpha ** (M - 1) / complex(sp.N(derivative.subs(u, alpha), 70))

    small_sum = sum(weight(alpha) for alpha in small)
    full_sum = sum(weight(alpha) for alpha in roots)
    coeffs = constant_terms(R, M, order)
    truncated = sum(complex(c) * float(t_value) ** m for m, c in enumerate(coeffs))
    first = next((m for m, c in enumerate(coeffs[1:], 1) if c != 0), None)
    small_error = abs(small_sum - truncated)
    full_error = abs(full_sum)
    if small_error > 1e-12 or full_error > 1e-12:
        raise RuntimeError((name, small_error, full_error, small_sum, full_sum))
    print(
        f"{name}: M={M} degree={phi.degree()} first_nonzero={first} "
        f"small_minus_series={small_error:.3e} full_sum={full_error:.3e}"
    )


def one_sided_boundary():
    # Lambda=u^-2+u^-1: all positive constant terms vanish, but deg(Phi)=M,
    # so k=M-1 is the top Lagrange degree and the full residue sum equals 1.
    M = 2
    R = 1 + u
    t_value = sp.Rational(1, 1000)
    phi = sp.Poly(u**M - t_value * R, u)
    derivative = sp.diff(phi.as_expr(), u)
    roots = [complex(z) for z in sp.nroots(phi, n=70, maxsteps=300)]
    full_sum = sum(
        alpha ** (M - 1) / complex(sp.N(derivative.subs(u, alpha), 70))
        for alpha in roots
    )
    coeffs = constant_terms(R, M, 8)
    if coeffs != [1] + [0] * 8 or abs(full_sum - 1) > 1e-12:
        raise RuntimeError((coeffs, full_sum))
    print(
        "one_sided_boundary: M=2 degree=2 all_positive_CT_zero=True "
        f"full_sum_minus_one={abs(full_sum - 1):.3e}"
    )


def main():
    residue_check("symmetric", 1 + u**2, 1)
    residue_check("three_charge", 1 + u + u**3, 2)
    residue_check("dihedral_signed", 1 + u - u**3 - u**4, 2)
    residue_check("asymmetric_extent", 2 - u + 3 * u**2 + u**5, 3)
    one_sided_boundary()
    print("PASS")


if __name__ == "__main__":
    main()
