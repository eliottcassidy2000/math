# m_a = Lebesgue measure of {x in [0,1): the first a significant digits of 10^x are all nonzero}
#     = sum_{s zeroless a-digit} log10(1+1/s)   (= density of n with first a digits of 2^n zeroless, by Weyl since log10 2 irrational)
# computed with the Kempner-block power-sum recursion (Baillie / Burnol style):
#   u_{a,j} = sum_{s zeroless a-digit} s^-j ;  u_{a+1,j} = sum_{m>=0} binom(-j,m) P_m 10^{-j-m} u_{a,j+m},  P_m = sum_{d=1}^9 d^m
from mpmath import mp, mpf, binomial, log, nstr, power
import itertools, json
mp.dps = 40
JMAX = 60
A0 = 4
def direct_u(a, J):
    u = [mpf(0)]*(J+1)
    for digs in itertools.product(range(1,10), repeat=a):
        s = 0
        for d in digs: s = 10*s + d
        x = mpf(1)/s; p = x
        for j in range(1, J+1):
            u[j] += p; p *= x
    return u
P = [sum(mpf(d)**m for d in range(1,10)) for m in range(JMAX+80)]
u = direct_u(A0, JMAX+70)
m_a = {}
def m_from_u(u, J):
    tot = mpf(0)
    for j in range(1, J+1):
        tot += (-1)**(j+1) * u[j] / j
    return tot / log(10)
# exact small a by direct summation
for a in range(1, A0+1):
    ua = direct_u(a, 50) if a < A0 else u
    m_a[a] = m_from_u(ua, 50) if a >= 2 else mpf(1)
# check a=1: sum_{s=1}^9 log10(1+1/s) = 1 exactly
# (the series for s=1 converges slowly; set exactly)
cur = u
Jcur = JMAX + 70
for a in range(A0+1, 61):
    Jn = max(JMAX, Jcur - 2)
    new = [mpf(0)]*(Jn+1)
    for j in range(1, Jn+1):
        tot = mpf(0)
        for m in range(0, Jcur - j + 1):
            term = binomial(-j, m) * P[m] * power(10, -j-m) * cur[j+m]
            tot += term
            if m > 8 and abs(term) < mpf(10)**(-45) * abs(tot): break
        new[j] = tot
    cur = new; Jcur = Jn
    m_a[a] = m_from_u(cur, min(Jcur, 30))
# direct check a=5 by brute force (9^5 = 59049 terms)
import math
b5 = sum(math.log10(1 + 1/s) for s in range(10000, 100000) if '0' not in str(s))
print('check a=5: recursion %s  brute %.15f' % (nstr(m_a[5], 18), b5))
lam = {a: m_a[a] / mpf('0.9')**a for a in m_a}
for a in [1,2,3,4,5,6,8,10,15,20,25,30,40,50,60]:
    print('a=%2d  m_a=%s   m_a/0.9^a=%s' % (a, nstr(m_a[a], 20), nstr(lam[a], 20)))
json.dump({str(a): str(m_a[a]) for a in m_a}, open('leading_measure.json','w'), indent=0)
