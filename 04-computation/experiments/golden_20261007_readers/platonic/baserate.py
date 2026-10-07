"""Base-rate test: Platonic / polyhedral numbers vs. the convergents and semiconvergents of log_2 3.
Statistic: |P ∩ C(alpha)|, where C(alpha) = numerators and denominators (<= NMAX) of all
convergents and semiconvergents (intermediate fractions) of alpha.  Null: alpha uniform in
[1.5, 5/3] (same first partial quotients [1;1,1,...] as log_2 3) and, separately, in [1.2, 1.9]."""
import random, math
from fractions import Fraction
NMAX = 120
def cf(alpha, n=40):
    a = []; x = alpha
    for _ in range(n):
        q = math.floor(x); a.append(q); f = x - q
        if f < 1e-12: break
        x = 1 / f
    return a
def semis(alpha):
    a = cf(alpha)
    p_2, q_2, p_1, q_1 = 0, 1, 1, 0
    S = set()
    for i, ai in enumerate(a):
        for j in range(1 if i > 0 else ai, ai + 1):   # intermediate fractions j = 1..a_i (j = a_i is the convergent)
            p = j * p_1 + p_2; q = j * q_1 + q_2
            if p <= NMAX: S.add(p)
            if q <= NMAX: S.add(q)
        p_2, q_2, p_1, q_1 = p_1, q_1, ai * p_1 + p_2, ai * q_1 + q_2
        if q_1 > NMAX and p_1 > NMAX: break
    return S
alpha0 = math.log2(3)
C0 = semis(alpha0)
print("C(log2 3) up to", NMAX, ":", sorted(C0))
sets = {
 "Platonic V/E/F + rotation orders + binary orders": {4, 6, 8, 12, 20, 30, 24, 60, 48, 120},
 "  ... plus Galois primes 5,7,11": {4, 5, 6, 7, 8, 11, 12, 20, 24, 30, 48, 60, 120},
 "rotation group orders only {12,24,60}": {12, 24, 60},
 "Ellison x,y values {13,14,16,19,27,8,9,10,12,17}": {13, 14, 16, 19, 27, 8, 9, 10, 12, 17},
}
random.seed(20261007)
N = 200000
for name, P in sets.items():
    obs = len(P & C0)
    for lo, hi in [(1.5, 5/3), (1.2, 1.9)]:
        cnt_ge = 0; tot = 0
        for _ in range(N):
            al = random.uniform(lo, hi)
            v = len(P & semis(al)); tot += v
            if v >= obs: cnt_ge += 1
        print("%-52s obs=%d (%s)  null[%g,%g]: mean=%.2f  P(>=obs)=%.3f" % (name, obs, sorted(P & C0), lo, hi, tot / N, cnt_ge / N))
