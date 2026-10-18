#!/usr/bin/env python3
"""Unflavored Schur index of 4d N=4 SU(2) SYM (N=2 vector + adjoint hyper), in s=q^{1/2}:
I = oint dz/(2 pi i z) (1/2)(1-z^2)(1-z^-2) prod_{n>=0,b in {2,0,-2}} (1-s^{2n+1}z^b)^{-2}
                                           prod_{n>=1,b} (1-s^{2n}z^b)^{2}."""
import sys
from collections import defaultdict
Smax = int(sys.argv[1]) if len(sys.argv) > 1 else 16   # max power of s
ser = {(0, 0): 1}
def mul_factor(ser, a, b, power):
    # multiply by (1 - s^a z^b)^power, power = +-1
    out = defaultdict(int)
    if power == 1:
        for (sp, zp), c in ser.items():
            out[(sp, zp)] += c
            if sp + a <= Smax:
                out[(sp + a, zp + b)] -= c
    else:
        # 1/(1 - x) = sum x^k
        for (sp, zp), c in ser.items():
            k = 0
            while sp + k * a <= Smax:
                out[(sp + k * a, zp + k * b)] += c
                k += 1
    return {kk: v for kk, v in out.items() if v}
for n in range(0, Smax):
    for b in (2, 0, -2):
        if 2 * n + 1 <= Smax:
            for _ in range(2):
                ser = mul_factor(ser, 2 * n + 1, b, -1)
        if n >= 1 and 2 * n <= Smax:
            for _ in range(2):
                ser = mul_factor(ser, 2 * n, b, 1)
# Haar measure (1/2)(1-z^2)(1-z^-2) = (1/2)(2 - z^2 - z^-2): take z^0 coefficient
res = defaultdict(int)
for (sp, zp), c in ser.items():
    if zp == 0: res[sp] += 2 * c
    if zp == -2: res[sp] -= c   # times z^2
    if zp == 2: res[sp] -= c    # times z^-2
coef = [res.get(i, 0) // 2 for i in range(Smax + 1)]
print("Schur index coefficients in q^{1/2}:", coef)
print("|coeffs| by weight (D = i/2):", [abs(x) for x in coef])
