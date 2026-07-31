#!/usr/bin/env python3
"""Verify the reductions derived by hand for S(4)."""
import mpmath as mm
from mpmath import mp

mp.dps = 60


def F(z):
    w = mm.sqrt(1 - z)
    m = (1 - w) / (1 + w)
    return mm.sqrt(2 / (1 + w)) * (2 / mm.pi) * mm.ellipk(m)


S4 = mm.quad(lambda t: F(t ** 4), [0, 1], maxdegree=14)
print("S(4) reference          =", mm.nstr(S4, 45))

# form 1:  S(4) = (sqrt2/pi) int_0^1 [atanh b - atan b] (1-b^4)^{-3/4} db
f1 = (mm.sqrt(2) / mm.pi) * mm.quad(
    lambda b: (mm.atanh(b) - mm.atan(b)) / (1 - b ** 4) ** mm.mpf(0.75), [0, 1], maxdegree=14)
print("form1 atanh-atan quartic=", mm.nstr(f1, 45), "  diff", mm.nstr(abs(f1 - S4), 5))

# form 2:  S(4) = (2 sqrt2/pi) 3F2(1,1,3/4; 5/4,7/4; 1)
f2 = (2 * mm.sqrt(2) / mm.pi) * mm.hyper([1, 1, mm.mpf(3) / 4],
                                         [mm.mpf(5) / 4, mm.mpf(7) / 4], 1)
print("form2 3F2(1,1,3/4;..)   =", mm.nstr(f2, 45), "  diff", mm.nstr(abs(f2 - S4), 5))

# form 3:  S(4) = (4/(pi sqrt2)) int_0^inf x^2 R(x)/(1+x^4) dx,  R=int_0^x (1+y^4)^{-3/4}
R = lambda x: mm.quad(lambda y: (1 + y ** 4) ** mm.mpf(-0.75), [0, x])
f3 = (4 / (mm.pi * mm.sqrt(2))) * mm.quad(
    lambda x: x ** 2 * R(x) / (1 + x ** 4), [0, 1, mm.inf], maxdegree=8)
print("form3 double integral   =", mm.nstr(f3, 30), "  diff", mm.nstr(abs(f3 - S4), 5))

# form 4:  S(4) = (2 sqrt2/pi) int_0^inf mu (1+mu^4)^{-3/4} P(mu) dmu, P=int_0^mu (1+s^4)^{-1}
P = lambda x: mm.quad(lambda s: 1 / (1 + s ** 4), [0, x])
f4 = (2 * mm.sqrt(2) / mm.pi) * mm.quad(
    lambda u: u * (1 + u ** 4) ** mm.mpf(-0.75) * P(u), [0, 1, mm.inf], maxdegree=8)
print("form4                   =", mm.nstr(f4, 30), "  diff", mm.nstr(abs(f4 - S4), 5))

# the analogous 3F2 for general k, guessed from the k=4 Thomae step:
# S(k) =? c_k * 3F2(...)  -- test the direct 3F2 (iii) for sanity
for k in (1, 2, 3, 4, 5):
    v = mm.hyper([mm.mpf(1) / 4, mm.mpf(3) / 4, mm.mpf(1) / k],
                 [1, 1 + mm.mpf(1) / k], 1)
    s = mm.quad(lambda t: F(t ** k), [0, 1], maxdegree=14)
    print(f"  k={k}: 3F2(iii) - S(k) = {mm.nstr(abs(v-s),5)}")
