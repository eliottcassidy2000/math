#!/usr/bin/env python3
"""Independent check: Kac-Wakimoto (integral-Weyl-group) vacuum character of L_k(sl2),
k = -2 + p/u admissible.  ch = N/D,
  N = sum_j q^{p u j^2 + u j} (x^{p j} - x^{-p j - 1}),
  D = (1 - x^{-1}) prod_{n>=1} (1-q^n)(1-x q^n)(1-x^{-1} q^n),
x = e^{alpha_1} (h_0-charge 2).  Prints dim L_D (summed over charge) and the charge-0 part."""
import sys
from fractions import Fraction
from collections import defaultdict

def kw(p, u, Dmax):
    # numerator as dict (qpow, xpow) -> coeff
    N = defaultdict(int)
    J = 0
    while True:
        done = True
        for j in ({J, -J} if J else {0}):
            e = p * u * j * j + u * j
            if e <= Dmax:
                done = False
                N[(e, p * j)] += 1
                N[(e, -p * j - 1)] -= 1
        if done and J > 2:
            break
        J += 1
    # divide by (1 - x^{-1}) per q-power:  A(x) = B(x)(1 - x^{-1}) -> B = sum_{m} ... (exact division)
    B = defaultdict(int)
    for D in range(Dmax + 1):
        poly = {xp: c for (qp, xp), c in N.items() if qp == D and c}
        if not poly:
            continue
        # B(x) = A(x)/(1 - x^{-1}); B_m - B_{m+1} = A_m  ... solve from top degree down: B_m = A_m + B_{m+1}
        lo, hi = min(poly), max(poly)
        acc = 0
        res = {}
        for m in range(hi, lo - 1, -1):
            acc = poly.get(m, 0) + acc
            res[m] = acc
        assert acc == 0 or True
        # verify exactness: B_{lo-1} contribution should vanish -> B_lo - 0 = A_{lo-1}=0 means acc at lo-1 = 0
        # (since B is a Laurent polynomial, the running sum must return to 0 below lo)
        assert sum(poly.values()) == 0, "numerator not divisible by (1 - x^-1)"
        for m, c in res.items():
            if c:
                B[(D, m)] += c
    # multiply by prod 1/((1-q^n)(1-xq^n)(1-x^{-1}q^n))
    ser = dict(B)
    for n in range(1, Dmax + 1):
        for (dx) in ((0,), (1,), (-1,)):
            d = dx[0]
            new = defaultdict(int)
            for (qp, xp), c in ser.items():
                k = 0
                while qp + k * n <= Dmax:
                    new[(qp + k * n, xp + k * d)] += c
                    k += 1
            ser = {kk: v for kk, v in new.items() if v}
    tot = [0] * (Dmax + 1)
    ch0 = [0] * (Dmax + 1)
    for (qp, xp), c in ser.items():
        tot[qp] += c
        if xp == 0:
            ch0[qp] += c
    return tot, ch0, ser

if __name__ == '__main__':
    p, u, Dmax = map(int, sys.argv[1:4])
    tot, ch0, ser = kw(p, u, Dmax)
    print("k =", Fraction(p, u) - 2)
    print("dim L_D:", tot)
    print("charge-0:", ch0)
    print("min coeff:", min(ser.values()))
