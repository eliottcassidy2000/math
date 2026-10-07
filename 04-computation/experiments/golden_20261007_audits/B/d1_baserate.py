# Independent re-run of the Platonic base-rate statistic (own CF code with exact Fractions of a
# 60-bit random alpha), seed differs from the original.
import random
from fractions import Fraction
from math import log2
NMAX = 120
def semis(alpha):
    S = set(); x = alpha
    a = []
    for _ in range(30):
        q = x.numerator // x.denominator; a.append(q); f = x - q
        if f == 0: break
        x = 1 / f
    p2, q2, p1, q1 = 0, 1, 1, 0
    for i, ai in enumerate(a):
        for j in range(1 if i > 0 else ai, ai + 1):
            p, q = j * p1 + p2, j * q1 + q2
            if p <= NMAX: S.add(p)
            if q <= NMAX: S.add(q)
        p2, q2, p1, q1 = p1, q1, ai * p1 + p2, ai * q1 + q2
        if p1 > NMAX and q1 > NMAX: break
    return S
alpha0 = Fraction(log2(3)).limit_denominator(10**15)
C0 = semis(alpha0)
P = {4, 6, 8, 12, 20, 30, 24, 60, 48, 120}
obs = len(P & C0)
print("C(log2 3) <= 120:", sorted(C0), " obs =", obs, sorted(P & C0))
random.seed(123457)
N = 100000; ge = 0; tot = 0
for _ in range(N):
    al = Fraction(3, 2) + Fraction(random.getrandbits(60), 2**60) * Fraction(1, 6)
    v = len(P & semis(al)); tot += v; ge += (v >= obs)
print(f"null U[1.5,5/3]: mean {tot/N:.3f}  P(>=obs) = {ge/N:.4f}")
