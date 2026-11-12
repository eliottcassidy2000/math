"""Auditor's own checks of Part B (independent code paths; does not import the packet's code).
 1. exact f_n by an integer Burnside sum over odd-part partitions (own enumeration); A000568 small terms.
 2. per cycle type: bit loss D vs the claimed level bounds; count of permutations of defect d vs C(n,2)^d/d!.
 3. T_K(n) vs B_K(n) for (12,0), (20,1), (30,2) and a sweep.
 4. N_odd(n), odd n <= 13, against A049313; clock-shift identity; first correction.
 5. leading digits of a(n) by mpmath (log-gamma at high precision), n = 300, 1000, 1e4, 1e5, 1e6, 1e9, 1e12, 1e100, 1e1000.
 6. Prop B2 moments; strong tournaments ratio; section-4 heuristic densities.
"""
from fractions import Fraction
from math import comb, factorial, gcd, ceil, floor
import sys

def oddparts(n, mx=None):
    if mx is None: mx = n if n % 2 else n - 1
    if n == 0: yield (); return
    for p in range(min(n, mx), 0, -1):
        if p % 2:
            for r in oddparts(n - p, p): yield (p,) + r

def zee(lam):
    z = 1
    for a in set(lam):
        c = lam.count(a); z *= a**c * factorial(c)
    return z

def orbits(lam):
    return sum((a - 1)//2 for a in lam) + sum(gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i+1, len(lam)))

def a_exact(n):
    s = sum(Fraction(2**orbits(l), zee(l)) for l in oddparts(n))
    assert s.denominator == 1
    return s.numerator

out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

# 1
A = [a_exact(n) for n in range(1, 13)]
P("1. A000568(1..12) =", A)
P("   matches known 1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168:",
  A == [1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168])

# per-permutation data at order n, grouped by level
def level_data(n):
    """returns dict k -> list of (count, D) over odd cycle types (count = number of permutations)"""
    d = {}
    for lam in oddparts(n):
        m = sum(a for a in lam if a > 1); r = sum(1 for a in lam if a > 1)
        k = (m - r)//2
        cnt = Fraction(factorial(n), zee(lam)); assert cnt.denominator == 1
        D = comb(n, 2) - orbits(lam)
        d.setdefault(k, []).append((cnt.numerator, D, m))
    return d

# 2
bad = 0; badc = 0; ncase = 0
for n in range(3, 41):
    d = level_data(n)
    for k, lst in d.items():
        if k == 0: continue
        bound = 2*k*n - 3*k*k - k if 3*k <= n else ceil(n*(n-1)/3)
        for cnt, D, m in lst:
            ncase += 1
            if D < bound: bad += 1
            if D < 2*k*(n - m) + Fraction(2, 3)*comb(m, 2): bad += 1
        tot = sum(c for c, D, m in lst)
        # all permutations of defect 2k (not only odd ones) <= C(n,2)^(2k)/(2k)!
        # unsigned Stirling c(n, n-2k)
    # unsigned Stirling numbers via e_d(1..n-1)
    e = [1] + [0]*n
    for x in range(1, n):
        for dd in range(n, 0, -1): e[dd] += e[dd-1]*x
    for dd in range(1, n):
        if Fraction(e[dd]) > Fraction(comb(n, 2)**dd, factorial(dd)): badc += 1
P("2. bit-loss bound by level (n=3..40, every odd cycle type): cases", ncase, "violations", bad,
  "; e_d(1..n-1) <= C(n,2)^d/d! violations:", badc)

def f_exact(n):
    return Fraction(a_exact(n)*factorial(n), 2**comb(n, 2))

def L(k, n):
    d = level_data(n)
    return sum((Fraction(c, 2**D) for c, D, m in d.get(k, [])), Fraction(0))

def t(k, n):
    return Fraction(comb(n, 2)**(2*k) * 2**(3*k*k + k), factorial(2*k) * 2**(2*k*n))

def B(K, n):
    return sum((t(k, n) for k in range(K+1, n//3 + 1)), Fraction(0)) + Fraction(factorial(n), 2**ceil(Fraction(n*(n-1), 3)))

def sci(x):
    from mpmath import mp, mpf, log10
    mp.dps = 30
    return mp.nstr(mpf(x.numerator)/mpf(x.denominator), 6)

# 3
for (n, K) in [(12, 0), (20, 1), (30, 2)]:
    f = f_exact(n); TK = f - sum(L(k, n) for k in range(0, K+1)); BK = B(K, n)
    P("3. n=%d K=%d  T_K=%s  B_K=%s  T<=B: %s  T>=0: %s" % (n, K, sci(TK), sci(BK), TK <= BK, TK >= 0))
viol = 0; cnt = 0
for n in range(3, 41):
    f = f_exact(n); part = Fraction(0)
    for K in range(0, n//2 + 1):
        part += L(K, n); TK = f - part; cnt += 1
        if TK < 0 or TK > B(K, n): viol += 1
    # log-convexity of t_k
for n in range(6, 200):
    for k in range(2, n//3):
        if t(k+1, n)*t(k-1, n) < t(k, n)**2: viol += 1
P("   sweep 3<=n<=40, all K: cases", cnt, "violations", viol, "(incl. log-convexity of t_k, n<200)")
P("   B_0(n) < 1 for n = 8..60:", all(B(0, n) < 1 for n in range(8, 61)), "; B_0(7) =", sci(B(0, 7)))
P("   n=100: T_0 =", sci(f_exact(100) - 1), " B_0 =", sci(B(0, 100)))

# 4
A049313 = {3: 1, 5: 2, 7: 12, 9: 792, 11: 886288, 13: 11856006240}
for n in sorted(A049313):
    No = sum(Fraction(2**(orbits(l) - len(l) + 1), zee(l)) for l in oddparts(n))
    # clock shift: N_odd n!/2^C(n-1,2) = sum_k 4^k L_k(n)
    lhs = No * factorial(n) / 2**comb(n-1, 2)
    rhs = sum(4**k * L(k, n) for k in range(0, n//2 + 1))
    P("4. n=%2d N_odd=%s A049313=%d equal=%s clock-shift identity=%s" % (n, No, A049313[n], No == A049313[n], lhs == rhs))

# 5 digits via mpmath
from mpmath import mp, mpf, log10, loggamma, log, floor as mfloor, power
def lead_digits(n, extra_dps=60):
    digits_n = len(str(n))
    mp.dps = 2*digits_n + extra_dps
    l10 = (mpf(n)*(n-1)/2)*log10(2) - loggamma(mpf(n)+1)/log(10)
    # correction: level 1 (one 3-cycle) relative size (16/3)(n)_3 4^-n; add when not negligible
    if n <= 5000:
        L1 = mpf(16)/3*mpf(n)*(n-1)*(n-2)*power(4, -n)
        l10 += log10(1 + L1)
    E = int(mfloor(l10)); mant = power(10, l10 - E)
    return E, mp.nstr(mant, 32, strip_zeros=False)
for n in [300, 1000, 10**4, 10**5, 10**6, 10**9, 10**12, 10**100, 10**1000]:
    E, m = lead_digits(n)
    Es = str(E)
    P("5. n=%s a(n) = %s x 10^%s" % (n if n < 10**13 else "10^%d" % (len(str(n))-1), m, Es if len(Es) < 40 else Es[:20] + "..." + Es[-8:] + " (%d digits)" % len(Es)))

# 6 Prop B2 moments (exact sums truncated at |j| <= 200)
mp.dps = 50
def mom_tri(M):
    Z = sum(power(2, -mpf(j)*(j+1)/2) for j in range(-200, 201))
    return sum(power(2, M*j - mpf(j)*(j+1)/2) for j in range(-200, 201))/Z
worst = max(abs(mom_tri(M)/power(2, comb(M, 2)) - 1) for M in range(0, 13))
P("6. E[2^(M J_1)] / 2^C(M,2) - 1, M=0..12: max", mp.nstr(worst, 3))
def mom_phi(M, ph):
    Z = sum(power(2, -mpf(j)**2/(2*ph)) for j in range(-400, 401))
    return sum(power(2, M*j - mpf(j)**2/(2*ph)) for j in range(-400, 401))/Z
worst = max(abs(mom_phi(M, ph)/power(2, mpf(ph)*M*M/2) - 1) for ph in (2, 4, 6) for M in range(-3, 7))
P("   E[2^(M J_d)] / 2^(phi M^2/2) - 1, phi in {2,4,6}, M=-3..6: max", mp.nstr(worst, 3))
# strong tournaments
NN = 41
Aser = [1] + [a_exact(n) for n in range(1, NN)]
# A(x) = 1/(1 - S(x)) with A(0)=1  =>  S = 1 - 1/A
inv = [Fraction(0)]*NN; inv[0] = Fraction(1)
for i in range(1, NN): inv[i] = -sum(Aser[j]*inv[i-j] for j in range(1, i+1))
S = [-x for x in inv]; S[0] = 0
P("   strong s(n), n=1..8:", [int(x) for x in S[1:9]])
for n in (10, 20, 40):
    r = 1 - Fraction(S[n]) / Aser[n]
    P("   n=%d 1-s/a=%s  n 2^(2-n)=%s ratio=%s" % (n, sci(r), sci(Fraction(n, 2**(n-2))), sci(r / Fraction(n, 2**(n-2)))))
# section 4 densities
maxH = {5: 15, 6: 45, 7: 189, 8: 661}
for n in (5, 6, 7, 8):
    an = a_exact(n); nodd = (maxH[n] + 1)//2
    P("   order %d: a(n)=%d, odd values up to max H=%d: %d, density %.2f" % (n, an, maxH[n], nodd, an/nodd))
a18 = a_exact(18); mean18 = factorial(18)/2**17
P("   order 18: a(18)=%.4e, mean H=%.4e (log2 %.2f), a/(mean/2)=%.3e, a/mean=%.3e" % (a18, mean18, __import__('math').log2(mean18), a18/(mean18/2), a18/mean18))
open(sys.argv[1] if len(sys.argv) > 1 else "audit_census.out", "w").write("\n".join(out) + "\n")
