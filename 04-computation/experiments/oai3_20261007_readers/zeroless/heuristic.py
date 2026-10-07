# Refined heuristic for the number of zeroless 2^n.
from mpmath import mp, mpf, log10, floor, nstr, exp, log
import json
mp.dps = 40
Z = {int(k): int(v) for k, v in json.load(open('zk_values.json'))['Z'].items()}
mlead = {int(a): mpf(v) for a, v in json.load(open('leading_measure.json')).items()}
c = mpf('0.887694043115148264')          # lim Z_k (2/9)^k  (k=40 value, ~18 digits)
tau = mpf(5)/4 * c                          # lim (Z_k/T_k)/0.9^k
lam = mpf('1.0845342222539724319')          # lim m_a/0.9^a
print('tau = (5/4)c =', nstr(tau, 18), '  lambda =', nstr(lam, 18), '  tau*lambda =', nstr(tau*lam, 18))
print('tau*0.9 =', nstr(tau*mpf('0.9'), 12), ' lambda*0.9 =', nstr(lam*mpf('0.9'), 12))
L2 = log10(2)
def D(n): return int(floor(n * L2)) + 1
def trail(k):   # P(last k digits zeroless) for a uniformly random residue
    if k == 0: return mpf(1)
    if k in Z: return mpf(Z[k]) / (4 * mpf(5)**(k-1))
    return tau * mpf('0.9')**k
def lead(a):
    if a == 0: return mpf(1)
    if a in mlead: return mlead[a]
    return lam * mpf('0.9')**a
def P(n):
    d = D(n)
    if d == 1: return mpf(1)  # single digit 1,2,4,8
    k = d // 2
    return lead(d - k) * trail(k)
# calibration: n <= 86
E86 = sum(P(n) for n in range(0, 87))
print('model expected count for 0<=n<=86: %s   (actual: 36)' % nstr(E86, 8))
# beyond 86
def tail_sum(N0, N1=None):
    s = mpf(0); n = N0
    while True:
        t = P(n); s += t; n += 1
        if (N1 is not None and n > N1) or t < mpf(10)**-30: break
    return s
E87 = tail_sum(87)
print('model expected count for n>=87: %s ; Poisson P(none) = %s' % (nstr(E87, 8), nstr(exp(-E87), 6)))
# simple version tau*lam*0.9^D
E87s = sum(tau*lam*mpf('0.9')**D(n) for n in range(87, 4000))
print('asymptotic-form tau*lambda*0.9^D(n), n>=87: %s' % nstr(E87s, 8))
naive = sum(mpf('0.9')**D(n) for n in range(87, 4000))
naive2 = sum(mpf('0.9')**(D(n)-2) for n in range(87, 4000))
print('naive 0.9^D: %s ; naive 0.9^(D-2) (first/last digit never 0): %s' % (nstr(naive, 8), nstr(naive2, 8)))
for N in [10**3, 10**4, 46*10**6, 10**10, 11*10**10, 7879942137257]:
    # tail beyond N: geometric in n with ratio 0.9^log10(2)
    r = mpf('0.9')**L2
    t = tau*lam*mpf('0.9')**(N*L2 + 1) / (1 - r)
    print('expected number with n > %d : 10^(%s)' % (N, nstr(log10(t), 10)))
# expected count in windows, and "first-zero" positional model for the verification statistics
# --- cluster adjustment (HEURISTIC): 2x is zeroless iff zeroless x has no '5' immediately left of a digit 1-4.
# Random-digit transfer matrix for "all digits nonzero and no 5[1-4]": largest eigenvalue rho2 solves r^2 - 0.9 r + 0.04 = 0.
from mpmath import sqrt
rho2 = (mpf('0.9') + sqrt(mpf('0.65'))) / 2
print('rho2 =', nstr(rho2, 15), ' conditional P(2^(n+1) zl | 2^n zl) ~ (rho2/0.9)^D, at D=27:', nstr((rho2/mpf('0.9'))**27, 6))
Ecl = sum(tau*lam*(mpf('0.9')**D(n) - rho2**D(n-1)) for n in range(87, 4000))
print('expected number of zeroless CLUSTERS (runs) starting at n>=87: %s ; Poisson P(none) = %s' % (nstr(Ecl, 6), nstr(exp(-Ecl), 4)))
zl86 = [0,1,2,3,4,5,6,7,8,9,13,14,15,16,18,19,24,25,27,28,31,32,33,34,35,36,37,39,49,51,67,72,76,77,81,86]
runs = sum(1 for i, n in enumerate(zl86) if i == 0 or zl86[i-1] != n-1)
Ecl86 = sum(P(n) for n in range(0, 87)) - sum(tau*lam*rho2**D(n-1) for n in range(1, 87))
print('actual runs for n<=86: %d ; model runs (rough, small-D regime): %s' % (runs, nstr(Ecl86, 6)))
