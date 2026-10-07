# Independent audit of THM-4568 (P_min least element of #107 Lemma 4.1 hypotheses)
from fractions import Fraction as Fr
import math, sys
import mpmath as mp

N = 70
# sigma_m = prod_{j=1}^{m-1} (1 + 1/(3j))
sig = {1: Fr(1)}
for m in range(2, 12*N):
    sig[m] = sig[m-1] * Fr(3*(m-1)+1, 3*(m-1))
def P(a, b):
    m, M = min(a, b), max(a, b)
    return sig[m] * Fr(2*M + m - 1, 2)

bad = []
# (B)
for b in range(1, 8*N):
    if P(1, b) != b: bad.append(('B', b))
# (C) in second arg for all a (first arg by symmetry)
for a in range(1, 2*N):
    for b in range(2, 8*N):
        if 2*P(a, b) < P(a, b+1) + P(a, b-1): bad.append(('C', a, b))
# (T)
for a in range(1, 2*N):
    for h in range(1, 2*N):
        if P(a, 3*h + a - 1) < 3*P(a, h): bad.append(('T', a, h))
print("admissibility violations:", bad[:5], len(bad))

# rank bound P^3 <= (a+b-1)^4 exact; equality cases
eq = []; viol = []
for a in range(1, 2*N):
    for b in range(a, 8*N):
        L = P(a, b)**3; R = Fr(a+b-1)**4
        if L > R: viol.append((a, b))
        if L == R: eq.append((a, b))
print("rank-bound violations:", viol[:5], "equalities:", eq[:5])

# which steps of the paper's chain are equalities on P_min
strict_cube = [m for m in range(1, 20) if (1 + Fr(1, 3*m))**3 == 1 + Fr(1, m)]
print("(1+1/(3m))^3 == 1+1/m for m in 1..19:", strict_cube, " -> the cube step is STRICT")
print("D_a^min / a^(4/3) at a=1,2,10,100:", [float(P(a, a)) / a**(4/3) for a in (1, 2, 10, 100)])

# gadget barrier: P_min(a, m h + c) >= m P_min(a,h), c = ceil((m-1)(a-1)/2), m=2..10
gv = []
for m in range(2, 11):
    for a in range(1, 120):
        c = -(-(m-1)*(a-1)//2)
        for h in range(1, 120):
            B = m*h + c
            if B >= 12*N: continue
            if P(a, B) < m*P(a, h): gv.append((m, a, h))
print("min-shift gadget violations (m<=10, a,h<120):", gv[:5], len(gv))
# and is c = ceil((m-1)(a-1)/2) - 1 refuted by P_min itself? (it should be for large h when (m-1)(a-1) odd/even...)
gv2 = []
for m in (2, 3):
    for a in range(2, 30):
        c = -(-(m-1)*(a-1)//2) - 1
        h = 200
        B = m*h + c
        if P(a, B) < m*P(a, h): gv2.append((m, a))
print("shift one less than minimal violated by P_min at (m,a) [h=200]:", gv2[:12], len(gv2))

# OEIS A004991 check: 9^(m-1) sigma_m
print("9^(m-1) sigma_m:", [9**(m-1)*sig[m] for m in range(1, 11)])
# generating function (1+x)(1-x)^(-7/3)
c73 = [Fr(1)]
for n in range(1, 60): c73.append(c73[-1]*(Fr(n-1) + Fr(7, 3))/n)
print("GF check:", all(P(a, a) == c73[a-1] + (c73[a-2] if a >= 2 else 0) for a in range(1, 60)))

# omega_a with mpmath
mp.mp.dps = 60
def lD(a):
    a = mp.mpf(a)
    return mp.log((3*a-1)/2) + mp.loggamma(a+mp.mpf(1)/3) - mp.loggamma(mp.mpf(4)/3) - mp.loggamma(a)
def om(a): return 3*mp.log(2*mp.mpf(a)-1)/lD(a)
# cross-check lD vs exact for small a
print("lD check a=50:", mp.nstr(lD(50) - mp.log(mp.mpf(P(50,50).numerator)/P(50,50).denominator), 5))
for a in (2, 186, 187, 188, 189, 190, 191, 595974, 595975, 595976):
    print(f"omega_{a} = {mp.nstr(om(a), 12)}")
print("omega_2 exact = 3 ln3/ln(10/3) =", mp.nstr(3*mp.log(3)/mp.log(mp.mpf(10)/3), 10), " log2 7 =", mp.nstr(mp.log(7, 2), 10))
# monotonicity check over a range (sampled + full small range)
prev = om(2); mono = True
for a in range(3, 5000):
    v = om(a)
    if v >= prev: mono = False; print("non-monotone at", a); break
    prev = v
print("omega_a strictly decreasing for 2<=a<5000:", mono)
# K constant
K = mp.mpf(9)/4*(mp.log(2) - mp.mpf(3)/4*mp.log(3/(2*mp.gamma(mp.mpf(4)/3))))
print("K =", mp.nstr(K, 12))
for a in (10**3, 10**6, 10**12, 10**30):
    print(f"  a=1e{int(mp.log10(a))}: (omega_a - 9/4) ln a = {mp.nstr((om(a) - mp.mpf(9)/4)*mp.log(a), 10)}")
