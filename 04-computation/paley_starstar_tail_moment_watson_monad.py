#!/usr/bin/env python3
"""
THM-438 ADD-16 (monad-explorer, deep-research 18th session).

GOAL: make ADD-15's qualitative "the rho(x)e^x hump and the A088368(k)/k! hump are
ONE phenomenon" PRECISE and QUANTITATIVE -- as a WATSON-LEMMA bridge.

THE FREE FACTORIAL LAW mu_free (free CP of nu=e^{-x}dx, free cumulants kappa_n=n!).

SPECTRAL TAIL (ADD-15, Stokes-suppressed to all orders in 1/x):
  K(w)=1/w+R(w)-i*pi*w^{-2}e^{-1/w},  R(w)=sum_{n>=1} n! w^{n-1}=1+2w+6w^2+24w^3+...
  with s=1/w_r (real ->+oo), sigma=w_r (small):
        x        = 1/sigma + R(sigma)
        rho(x)   = e^{-1/w_r} = e^{-s}            (Stokes term, rho=-w_i/pi)
   =>   rho(x)e^x= e^{x-s} = e^{R(sigma)}.
  Expanding sigma as a series in t=1/x gives
        rho(x) ~ e^{-x} * sum_{j>=0} a_j x^{-j},  a_0 = e,   a_j = e*b_j.

MOMENTS  m_k = A088368(k) = INT_0^oo x^k rho(x) dx.
Watson's lemma applied to rho(x)=e^{-x} sum a_j x^{-j}:
        m_k ~ sum_j a_j * Gamma(k+1-j) = sum_j a_j (k-j)!
   =>   m_k/k! ~ sum_j a_j (k-j)!/k! = sum_j a_j/(k)_j = e*[ b_0 + b_1/k + b_2/(k(k-1)) + ...]
  where (k)_j = k(k-1)...(k-j+1) is the falling factorial.

  >>> THE BRIDGE: the SAME coefficients {a_j} = e*{b_j} govern BOTH the spectral tail
      of the density AND the large-k moment ratio.  The resurgent factorial series R
      is read once on the spectral side (tail of rho) and once on the cumulant side
      (large-k moments) -- they are Watson-lemma images of each other. <<<

Everything below is exact (Fraction) power-series arithmetic + a high-precision
numerical cross-check against the genuine ADD-15 parametric density.
"""
from fractions import Fraction as Fr
from math import factorial
import mpmath as mp

# --------------------------------------------------------------------------- #
#  minimal truncated power-series toolkit over Fraction (index = power)        #
# --------------------------------------------------------------------------- #
def trunc(a, N):
    a = a[:N + 1]
    return a + [Fr(0)] * (N + 1 - len(a))

def mul(a, b, N):
    c = [Fr(0)] * (N + 1)
    for i, ai in enumerate(a[:N + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[:N + 1 - i]):
            if bj:
                c[i + j] += ai * bj
    return c

def compose(a, b, N):
    """a(b(t)) with b[0]==0, truncated to order N."""
    assert b[0] == 0
    res = [Fr(0)] * (N + 1)
    res[0] = a[0]
    bp = [Fr(1)] + [Fr(0)] * N   # b^0
    for k in range(1, len(a)):
        bp = mul(bp, b, N)
        ak = a[k]
        if ak:
            for i in range(N + 1):
                res[i] += ak * bp[i]
    return res

def series_exp(a, N):
    """exp(a) with a[0]==0."""
    assert a[0] == 0
    res = [Fr(0)] * (N + 1)
    res[0] = Fr(1)
    term = [Fr(1)] + [Fr(0)] * N
    for k in range(1, N + 1):
        term = mul(term, a, N)
        term = [x / k for x in term]
        for i in range(N + 1):
            res[i] += term[i]
    return res

def revert(t_of_s, N):
    """given t = s + c2 s^2 + ... (t_of_s[0]=0, t_of_s[1]=1), return s(t) to order N."""
    assert t_of_s[0] == 0 and t_of_s[1] == 1
    s = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)   # initial s = t
    for _ in range(N + 1):
        comp = compose(t_of_s, s, N)         # t_of_s(s(t)) should equal t
        corr = [comp[i] - (Fr(1) if i == 1 else Fr(0)) for i in range(N + 1)]
        s = [s[i] - corr[i] for i in range(N + 1)]
    return s

# --------------------------------------------------------------------------- #
#  (A) A088368(k) = free moments, M(z)=F(zM), F(w)=sum n! w^n  (pure int)      #
# --------------------------------------------------------------------------- #
def free_moments(kmax):
    fact = [1] * (kmax + 1)
    for n in range(1, kmax + 1):
        fact[n] = fact[n - 1] * n
    def pmul(A, B, K):
        C = [0] * (K + 1)
        for i, ai in enumerate(A):
            if ai == 0 or i > K:
                continue
            for j, bj in enumerate(B):
                if i + j > K:
                    break
                if bj:
                    C[i + j] += ai * bj
        return C
    M = [0] * (kmax + 1); M[0] = 1
    for _ in range(kmax + 1):
        zM = [0] + M[:kmax]
        Fval = [0] * (kmax + 1); Fval[0] += fact[0]
        zMpow = [1] + [0] * kmax
        for n in range(1, kmax + 1):
            zMpow = pmul(zMpow, zM, kmax)
            if not any(zMpow):
                break
            fn = fact[n]
            for p in range(kmax + 1):
                if zMpow[p]:
                    Fval[p] += fn * zMpow[p]
        if Fval == M:
            break
        M = Fval
    return M

KMAX = 60
print("Computing free moments m_k = A088368(k) (k<=%d) via M=F(zM)..." % KMAX, flush=True)
m = free_moments(KMAX)
print("m[0..8] =", m[:9])
# ground truth (brute NC enumeration, this session): Sum_{NC(k)} prod|B|!
nc_truth = [1, 1, 3, 13, 69, 421, 2867, 21477, 175769]
print("matches brute NC-partition sum 1,1,3,13,69,421,2867,21477,175769:", m[:9] == nc_truth)

# --------------------------------------------------------------------------- #
#  (B) tail coefficients  rho*e^x = e * sum b_j t^j,  t=1/x                    #
#      x = (1 + sum_{n>=1} n! sigma^n)/sigma  =>  t = sigma/(1+sum n! sigma^n) #
#      rho*e^x = e * exp(P(sigma)),  P=sum_{n>=2} n! sigma^{n-1}               #
# --------------------------------------------------------------------------- #
N = 12
# Q(sigma) = 1 + sum_{n>=1} n! sigma^n  (so x = Q/sigma, t = sigma/Q)
Q = [Fr(0)] * (N + 3)
Q[0] = Fr(1)
for n in range(1, N + 3):
    Q[n] = Fr(factorial(n))
# t(sigma) = sigma / Q(sigma) = sigma * (1/Q)
# 1/Q :
invQ = [Fr(0)] * (N + 3); invQ[0] = Fr(1)
for k in range(1, N + 3):
    s = Fr(0)
    for j in range(1, k + 1):
        s += Q[j] * invQ[k - j]
    invQ[k] = -s  # since Q[0]=1
t_of_s = [Fr(0)] * (N + 3)
for k in range(N + 2):
    t_of_s[k + 1] = invQ[k]          # multiply by sigma
t_of_s = trunc(t_of_s, N)
s_of_t = revert(t_of_s, N)           # sigma(t)

# P(sigma) = sum_{n>=2} n! sigma^{n-1}  (note: constant term 0, linear coeff = 2)
P = [Fr(0)] * (N + 1)
for n in range(2, N + 2):
    if n - 1 <= N:
        P[n - 1] = Fr(factorial(n))
P_of_t = compose(P, s_of_t, N)       # P(sigma(t)) ; P_of_t[0]=0
b = series_exp(P_of_t, N)            # rho*e^x / e = exp(P) = sum b_j t^j

print("\n(B) tail coefficients  b_j  (rho*e^x = e * sum_j b_j x^{-j}):")
for j in range(N + 1):
    print(f"   b_{j} = {b[j]}   (~ {float(b[j]):.6f})")

# OEIS check on the integer/rational b_j numerators
print("\n   b_j numerators:", [b[j].numerator for j in range(N + 1)])
print("   b_j as floats :", [round(float(b[j]), 4) for j in range(N + 1)])

# --------------------------------------------------------------------------- #
#  (C1) WATSON BRIDGE: m_k/k!  vs  e*sum_{j<=J} b_j/(k)_j                      #
# --------------------------------------------------------------------------- #
mp.mp.dps = 40
bmp = [mp.mpf(b[j].numerator) / mp.mpf(b[j].denominator) for j in range(N + 1)]
print("\n(C1) WATSON BRIDGE:  m_k/k!  vs  e * sum_{j<=J} b_j / (k)_j")
print(f"{'k':>4} {'m_k/k! exact':>20} {'Watson J=4':>18} {'Watson J=8':>18} "
      f"{'Watson J=12':>18} {'relerr J=12':>13}")
for k in [10, 15, 20, 25, 30, 40, 50, 60]:
    ratio = mp.mpf(m[k]) / mp.factorial(k)
    def watson(J):
        s = mp.mpf(0)
        for j in range(min(J, k) + 1):     # (k)_j = 0 for j>k; asymptotic needs k>=J
            fall = mp.mpf(1)
            for i in range(j):
                fall *= (k - i)
            s += bmp[j] / fall
        return mp.e * s
    w4, w8, w12 = watson(4), watson(8), watson(12)
    rel = abs(ratio - w12) / ratio
    print(f"{k:>4} {mp.nstr(ratio,14):>20} {mp.nstr(w4,12):>18} {mp.nstr(w8,12):>18} "
          f"{mp.nstr(w12,12):>18} {mp.nstr(rel,3):>13}")

# --------------------------------------------------------------------------- #
#  (C2) cross-check: a_j against the GENUINE ADD-15 parametric density tail    #
#       rho*e^x = e^{R(w_r)} (R Borel/optimal-truncation summed), x=1/w_r+R    #
# --------------------------------------------------------------------------- #
def R_opt(wr):
    """R(w)=sum n! w^{n-1}, optimal (smallest-term) truncation for wr>0 small."""
    s = mp.mpf(0); prev = None
    for n in range(1, 400):
        term = mp.factorial(n) * wr ** (n - 1)
        if prev is not None and abs(term) > abs(prev):
            break
        s += term; prev = term
    return s
print("\n(C2) parametric density tail rho*e^x = e^{R(w_r)}  vs  series e*sum b_j x^{-j}:")
print(f"{'w_r':>8} {'x':>13} {'rho*e^x param':>18} {'series':>18} {'relerr':>11}")
for wr in [mp.mpf('0.05'), mp.mpf('0.03'), mp.mpf('0.02'), mp.mpf('0.012'), mp.mpf('0.008')]:
    Rv = R_opt(wr)
    x = 1 / wr + Rv
    param = mp.e ** Rv                       # rho*e^x
    ser = mp.e * sum(bmp[j] / x ** j for j in range(N + 1))
    rel = abs(param - ser) / param
    print(f"{mp.nstr(wr,3):>8} {mp.nstr(x,9):>13} {mp.nstr(param,11):>18} "
          f"{mp.nstr(ser,11):>18} {mp.nstr(rel,3):>11}")

# --------------------------------------------------------------------------- #
#  (D) RESURGENCE diagnostic: the moment-ratio asymptotic is ITSELF a          #
#      divergent (Gevrey-1) series -- the resurgence of R propagates through    #
#      Watson's lemma. Show b_j growth + optimal truncation of the moment       #
#      series at large fixed k.                                                 #
# --------------------------------------------------------------------------- #
print("\n(D) RESURGENCE of the tail/moment coefficients b_j:")
print(f"{'j':>4} {'b_j':>18} {'b_j/b_{j-1}':>14} {'b_j/j!':>16}")
for j in range(1, N + 1):
    rj = float(b[j]) / float(b[j - 1])
    bf = float(b[j]) / factorial(j)
    print(f"{j:>4} {float(b[j]):>18.4g} {rj:>14.4f} {bf:>16.4f}")

print("\n   optimal truncation of  m_k/k! = e*sum b_j/(k)_j  at fixed k")
print("   (error of best partial sum vs exact -> the resurgent floor):")
for k in [40, 50, 60]:
    ratio = mp.mpf(m[k]) / mp.factorial(k)
    best = None; bestJ = None
    partial = mp.mpf(0)
    for j in range(0, min(N, k) + 1):
        fall = mp.mpf(1)
        for i in range(j):
            fall *= (k - i)
        partial += mp.e * bmp[j] / fall
        err = abs(partial - ratio)
        if best is None or err < best:
            best, bestJ = err, j
    print(f"   k={k}: best partial sum at J={bestJ}, |err|={mp.nstr(best,3)} "
          f"(rel {mp.nstr(best/ratio,3)})")

# --------------------------------------------------------------------------- #
#  (E) INDEPENDENT cross-check of b_j by Lagrange-Burmann inversion.           #
#      t = sigma/Q(sigma) = sigma*phi(sigma), phi=1/Q.                          #
#      H(sigma)=exp(R(sigma)-1).  b_j=[t^j]H(sigma(t)) = (1/j)[s^{j-1}] H'(s) Q(s)^j.
# --------------------------------------------------------------------------- #
print("\n(E) Lagrange-Burmann cross-check of b_j (independent of the reversion):")
M2 = N + 2
# R(s)-1 = P(s) ; H=exp(P); H' = P' * H
Pser = [Fr(0)] * (M2 + 1)
for n in range(2, M2 + 2):
    if n - 1 <= M2:
        Pser[n - 1] = Fr(factorial(n))
Hser = series_exp(Pser, M2)               # exp(P)
# H'(s):
Hp = [Fr(0)] * (M2 + 1)
for k in range(1, M2 + 1):
    Hp[k - 1] = Hser[k] * k
# Q(s) coeff list to order M2
Qs = [Fr(0)] * (M2 + 1)
Qs[0] = Fr(1)
for n in range(1, M2 + 1):
    Qs[n] = Fr(factorial(n))
ok = True
bL = [Fr(1)]                              # b_0
for j in range(1, N + 1):
    # Q^j to order j-1
    Qpow = [Fr(1)] + [Fr(0)] * (j - 1 + 1)
    for _ in range(j):
        Qpow = mul(Qpow, Qs, j - 1)
    val = Fr(0)
    for i in range(j):                    # [s^{j-1}] Hp * Qpow
        val += Hp[i] * Qpow[j - 1 - i]
    bj = val / j
    bL.append(bj)
    if bj != b[j]:
        ok = False
        print(f"   MISMATCH at j={j}: lagrange {bj} vs reversion {b[j]}")
print("   Lagrange-Burmann b_j == reversion b_j for all j<=%d : %s" % (N, ok))

print("\nDONE.")
