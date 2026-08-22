#!/usr/bin/env python3
"""Numerical control only for THM-3213's analytically proved diagonal asymptotic."""

from math import comb

import mpmath as mp


mp.mp.dps = 80


def ordered_stirling_profile_scaled(n):
    # b[n,k] = k! S(n,k) / n! = [x^n](exp(x)-1)^k.
    row = [mp.mpf(1)]
    for m in range(1, n + 1):
        nxt = [mp.mpf(0)] * (m + 1)
        for k in range(1, m + 1):
            left = row[k] if k < len(row) else 0
            right = row[k - 1] if k - 1 < len(row) else 0
            nxt[k] = mp.mpf(k) * (left + right) / m
        row = nxt
    return row


def numerator_power(d):
    base = {
        (1, 0, 0): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (0, 1, 1): 1,
        (1, 1, 1): 3,
    }
    poly = {(0, 0, 0): 1}
    for _ in range(d):
        out = {}
        for a, va in poly.items():
            for b, vb in base.items():
                e = tuple(a[i] + b[i] for i in range(3))
                out[e] = out.get(e, 0) + va * vb
        poly = out
    return poly


def normalized_coordinate(row, d):
    n = len(row) - 1
    ans = mp.mpf(0)
    for (a, b, c), coeff in numerator_power(d).items():
        for q in range(0, n + 1):
            if q + max(a, b, c) > n:
                break
            ans += (
                coeff
                * comb(d + q - 1, d - 1)
                * row[q + a]
                * row[q + b]
                * row[q + c]
            )
    return ans


def predicted_constant(d):
    ell = mp.log(2)
    return 9**d / (
        (2 * ell) ** (d - 1) * 4 * mp.pi * mp.sqrt(3) * ell * (1 - ell)
    )


ratios_by_n = {}
for n in (40, 80, 160, 320, 640):
    row = ordered_stirling_profile_scaled(n)
    ell = mp.log(2)
    print(f"n={n}")
    ratios_by_n[n] = []
    for d in (1, 2, 3, 4):
        value = normalized_coordinate(row, d)
        scaled = value * ell ** (3 * n) * n / comb(n - 1, d - 1)
        constant = predicted_constant(d)
        ratio = scaled / constant
        ratios_by_n[n].append(ratio)
        print(
            f" d={d} scaled={mp.nstr(scaled, 20)} "
            f"constant={mp.nstr(constant, 20)} ratio={mp.nstr(ratio, 16)}"
        )


if any(abs(ratio - 1) >= mp.mpf(1) / 200 for ratio in ratios_by_n[640]):
    raise RuntimeError(('final asymptotic control outside 1/200', ratios_by_n[640]))
print('status=NUMERICAL_CONTROL_ONLY;proof=analytic_smooth_diagonal_in_THM3213')
print('scope=K1_cyclic_diagonal_fixed_d_1_through_4_n_40_through_640')
print('FAILED_CHECKS=NONE')
