#!/usr/bin/env python3
"""
Cancellation-gate atlas: eta powers and extremal Type II formal enumerators.

codex-2026-06-11-P2

This extends the pentagonal/[72,36,16] packet in two directions:

1. Eta powers:
   eta^1 is Euler pentagonal sparsity; eta^3 has Jacobi triangular sparsity;
   eta^24 is the Ramanujan Delta gate, dense but modularly controlled.

2. Type II formal weight enumerators:
   The Gleason ring gives a scalar extremal enumerator at every length 24m.
   For lengths 24..240 in this stored run the scalar enumerator stays integral
   and nonnegative.  Thus scalar feasibility is a weak gate: support realization
   is the hard problem already at length 72 and remains the correct obstruction
   language higher up.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb


ETA_N = 120
TYPEII_LENGTHS = list(range(24, 241, 24))


def poly_add(a: list[int], b: list[int]) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul_trunc(a: list[int], b: list[int], nmax: int) -> list[int]:
    out = [0] * (min(len(a) + len(b) - 1, nmax + 1))
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j > nmax:
                break
            out[i + j] += ai * bj
    return out


def eta_power_coeffs(power: int, nmax: int = ETA_N) -> list[int]:
    coeffs = [1] + [0] * nmax
    # prod_{m<=nmax} (1-q^m)^power, sufficient through q^nmax.
    for m in range(1, nmax + 1):
        factor = [0] * (power * m + 1)
        for j in range(power + 1):
            factor[j * m] = (-1) ** j * comb(power, j)
        coeffs = poly_mul_trunc(coeffs, factor, nmax)
    return coeffs


def triangular_eta3_coeff(n: int) -> int:
    # Jacobi: prod(1-q^m)^3 = sum_{k>=0} (-1)^k (2k+1) q^{k(k+1)/2}.
    k = 0
    while k * (k + 1) // 2 < n:
        k += 1
    if k * (k + 1) // 2 == n:
        return (-1) ** k * (2 * k + 1)
    return 0


def fadd(a: list[Fraction], b: list[Fraction], scale: Fraction = Fraction(1)) -> list[Fraction]:
    n = max(len(a), len(b))
    out = [Fraction(0)] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += scale * b[i]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def fmul(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def fpow(a: list[Fraction], e: int) -> list[Fraction]:
    out = [Fraction(1)]
    base = a[:]
    while e:
        if e & 1:
            out = fmul(out, base)
        e >>= 1
        if e:
            base = fmul(base, base)
    return out


def typeii_formal(n: int) -> dict[str, object]:
    if n % 24 != 0:
        raise ValueError("stored atlas uses n divisible by 24")
    m = n // 8
    target_d = 4 * (n // 24) + 4
    A = [Fraction(1), Fraction(14), Fraction(1)]
    B = [Fraction(0), Fraction(1), Fraction(-4), Fraction(6), Fraction(-4), Fraction(1)]
    basis = [fmul(fpow(A, m - 3 * j), fpow(B, j)) for j in range(m // 3 + 1)]
    W = basis[0][:]
    coeffs = [Fraction(1)]
    for r in range(1, target_d // 4):
        current = W[r] if r < len(W) else Fraction(0)
        c = -current / basis[r][r]
        coeffs.append(c)
        W = fadd(W, basis[r], c)
    integral = all(c.denominator == 1 for c in W)
    W_int = [int(c) for c in W] if integral else []
    neg = [(4 * i, c) for i, c in enumerate(W) if c < 0]
    first_weight = next((4 * i for i, c in enumerate(W) if i and c), None)
    first_coeff = W[first_weight // 4] if first_weight is not None else None
    return {
        "n": n,
        "d": target_d,
        "coeffs": coeffs,
        "integral": integral,
        "nonnegative": not neg,
        "first_weight": first_weight,
        "first_coeff": first_coeff,
        "sum": sum(W_int) if integral else None,
        "neg": neg,
    }


def fmt_fraction(x: Fraction | int | None) -> str:
    if x is None:
        return "-"
    if isinstance(x, int):
        return str(x)
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def print_eta_section() -> None:
    print("[1] Eta-power cancellation gates")
    for power in [1, 3, 8, 24]:
        coeffs = eta_power_coeffs(power, ETA_N)
        support = sum(1 for c in coeffs if c)
        max_abs = max(abs(c) for c in coeffs)
        print(
            f"  eta^{power:<2d}: support<=q^{ETA_N}: {support:3d}/{ETA_N+1} "
            f"density={support/(ETA_N+1):.3f} max|coeff|={max_abs}"
        )
        nz = [(i, c) for i, c in enumerate(coeffs) if c][:12]
        print("    first nonzero:", nz)
    eta3 = eta_power_coeffs(3, ETA_N)
    ok3 = all(eta3[n] == triangular_eta3_coeff(n) for n in range(ETA_N + 1))
    print(f"  Jacobi eta^3 triangular formula verified through q^{ETA_N}: {ok3}")
    print("  Delta/q = eta^24 is dense: modular control replaces literal sparsity.")
    print()


def print_typeii_section() -> None:
    print("[2] Extremal Type II formal weight-enumerator ladder")
    print("  n    d    A_d                 Gleason coeffs c_j (first few)        nonneg sum-ok")
    for n in TYPEII_LENGTHS:
        row = typeii_formal(n)
        coeffs = row["coeffs"]
        coeff_preview = ",".join(fmt_fraction(c) for c in coeffs[:6])
        if len(coeffs) > 6:
            coeff_preview += ",..."
        sum_ok = row["sum"] == 2 ** (n // 2)
        print(
            f"  {n:3d}  {row['d']:3d}  {fmt_fraction(row['first_coeff']):>18s}  "
            f"{coeff_preview:<38s} {str(row['nonnegative']):>6s} {str(sum_ok):>6s}"
        )
    print()
    print("  Landmarks:")
    for n in [24, 48, 72, 96, 120]:
        row = typeii_formal(n)
        print(
            f"    n={n}: formal extremal d={row['d']}, "
            f"A_d={fmt_fraction(row['first_coeff'])}"
        )
    print("  n=24 gives the Golay octad number 759; n=72 gives 249849.")
    print("  The scalar gate stays open through this ladder; support realization is separate.")
    print()


def print_cross_domain_section() -> None:
    print("[3] Cross-domain reading")
    print("  eta^1: sparse infinite product gate -> partition reciprocal.")
    print("  eta^3: Jacobi triangular gate -> still sparse, now weighted by 2k+1.")
    print("  eta^24: Delta gate -> dense tau coefficients but modularly controlled.")
    print("  Type II W_n: finite invariant-ring gate -> low-weight zeros by Gleason.")
    print("  Tutte/Greene direction: code weight enumerators are matroid partition")
    print("  functions, so support realization should be audited as a matroid gate,")
    print("  not just as a polynomial gate.")


def main() -> None:
    print("Cancellation-gate atlas: eta powers and Type II formal enumerators")
    print(f"ETA_N={ETA_N}; Type II lengths={TYPEII_LENGTHS[0]}..{TYPEII_LENGTHS[-1]}")
    print()
    print_eta_section()
    print_typeii_section()
    print_cross_domain_section()


if __name__ == "__main__":
    main()
