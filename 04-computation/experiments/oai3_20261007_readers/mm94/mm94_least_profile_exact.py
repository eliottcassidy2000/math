#!/usr/bin/env python3
"""
mm94 lane: the least element of the hypotheses of Lemma 4.1 (Diagonal growth) of
OpenAI, "An Upper Bound of 9/4 for the Matrix Multiplication Exponent" (openai/math #107).

Hypotheses (the admissible set A):  P : N_{>0}^2 -> R_{>0}
  (S)  P(a,b) = P(b,a)
  (B)  P(1,b) = b
  (C)  2P(a,b) >= P(a,b+1) + P(a,b-1)        (a >= 1, b >= 2)   [and by (S) in a]
  (T)  P(a,3h+a-1) >= 3 P(a,h)                (a, h >= 1)
Rank bound used by the paper:  P(a,b) <= (a+b-1)^(1/t), t = mean dot-product exponent.

Claim (mm94, PROVED in the report): A has a least element
  P_min(a,b) = sigma_m * (2M + m - 1)/2,  m = min(a,b), M = max(a,b),
  sigma_m = prod_{j=1}^{m-1} (1 + 1/(3j)) = [x^(m-1)] (1-x)^(-4/3) = Gamma(m+1/3)/(Gamma(4/3)Gamma(m)),
and P_min(a,b) <= (a+b-1)^(4/3) for all a,b.  So the paper's hypotheses + rank bound are
satisfiable exactly up to t = 3/4: the argument's ceiling is omega <= 9/4.

Everything below is exact rational arithmetic unless marked FLOAT.
"""
from fractions import Fraction as Fr
import math, random, sys

N_EXACT = int(sys.argv[1]) if len(sys.argv) > 1 else 90     # exact grid
N_FLOAT = int(sys.argv[2]) if len(sys.argv) > 2 else 1500   # float grid

# ---------------------------------------------------------------- sigma, P_min
def make_sigma(n):
    s = [None, Fr(1)]
    for a in range(2, n + 1):
        s.append(s[-1] * Fr(3 * a - 2, 3 * a - 3))
    return s

SIG = make_sigma(8 * N_EXACT + 10)

def Pmin(a, b):
    m, M = min(a, b), max(a, b)
    return SIG[m] * Fr(2 * M + m - 1, 2)

def Pflat(a, b):  # FLOAT: symmetrized profile of the three flattening-rank characters
    return math.sqrt(a * b * (a + b - 1))

fails = []
def check(cond, msg):
    if not cond:
        fails.append(msg)

# ---------------------------------------------------------------- (1) P_min in A (exact)
N = N_EXACT
for b in range(1, 4 * N):
    check(Pmin(1, b) == b, f"B fails b={b}")
for a in range(1, N + 1):
    for b in range(2, 4 * N):
        check(2 * Pmin(a, b) >= Pmin(a, b + 1) + Pmin(a, b - 1), f"C fails {a},{b}")
for a in range(1, N + 1):
    for h in range(1, N + 1):
        lhs, rhs = Pmin(a, 3 * h + a - 1), 3 * Pmin(a, h)
        check(lhs >= rhs, f"T fails {a},{h}")
        if h >= a:
            check(lhs == rhs, f"T not tight at h>=a {a},{h}")
print(f"[1] P_min satisfies (S),(B),(C),(T) on a<={N}, b<4N (exact): {'OK' if not fails else fails[:5]}")

# ---------------------------------------------------------------- (2) the paper's chain is tight on P_min
ok = True
for a in range(1, N + 1):
    D = Pmin(a, a)
    H = 2 * D / (3 * a - 1)
    ok &= (H == SIG[a])                                   # H_a = sigma_a
    ok &= (Pmin(a, a + 1) - Pmin(a, a) == SIG[a])          # Delta_{a,a} = H_a
    if a >= 2:
        ok &= (D - Pmin(a - 1, a - 1) == SIG[a] + SIG[a - 1])   # (4.1) with equality
print(f"[2] every inequality in the proof of Lemma 4.1 is an equality on P_min: {ok}")

# ---------------------------------------------------------------- (3) generating function (1+x)(1-x)^(-7/3)
c = [Fr(1)]
for n in range(1, N + 2):
    c.append(c[-1] * (Fr(n - 1) + Fr(7, 3)) / n)       # [x^n](1-x)^(-7/3)
gf_ok = all(Pmin(a, a) == c[a - 1] + (c[a - 2] if a >= 2 else 0) for a in range(1, N + 1))
s43 = [Fr(1)]
for n in range(1, N + 2):
    s43.append(s43[-1] * (Fr(n - 1) + Fr(4, 3)) / n)   # [x^n](1-x)^(-4/3)
gf2_ok = all(SIG[a] == s43[a - 1] for a in range(1, N + 1))
print(f"[3] sum_a D_a^min x^(a-1) = (1+x)(1-x)^(-7/3): {gf_ok};  sigma_a = [x^(a-1)](1-x)^(-4/3): {gf2_ok}")
print("    D_1..D_8 =", [str(Pmin(a, a)) for a in range(1, 9)])
print("    9^(a-1) sigma_a (integers?) =", [str(9 ** (a - 1) * SIG[a]) for a in range(1, 10)])

# ---------------------------------------------------------------- (4) rank bound at t = 3/4, exact: P^3 <= (a+b-1)^4
rb_ok, tight = True, []
for a in range(1, N + 1):
    for b in range(a, 4 * N):
        P = Pmin(a, b)
        R = Fr((a + b - 1) ** 4)
        if P ** 3 > R:
            rb_ok = False
        if P ** 3 == R:
            tight.append((a, b))
print(f"[4] P_min(a,b)^3 <= (a+b-1)^4 on a<={N}, a<=b<4N (exact): {rb_ok}; equality only at {tight[:5]}")

# ---------------------------------------------------------------- (5) extended gadgets (our extension of the paper's Lemma 3.2 construction)
# odd k-block alternating F,T,F,...,F:  P(a, k h + (k-1)(a-1)/2) >= k P(a,h)
# arbitrary alternating sizes:          P(a, sum_F h_i + sum_T (h'_j + a - 1)) >= sum P(a, .)
random.seed(94)
ext_ok = True
for k in (3, 5, 7, 9):
    for a in range(1, 40):
        for h in range(1, 40):
            B = k * h + (k - 1) * (a - 1) // 2
            ext_ok &= Pmin(a, B) >= k * Pmin(a, h)
for _ in range(4000):
    a = random.randint(1, 30)
    mF = random.randint(2, 4)
    hs = [random.randint(1, 30) for _ in range(2 * mF - 1)]
    B = sum(hs[0::2]) + sum(h + a - 1 for h in hs[1::2])
    ext_ok &= Pmin(a, B) >= sum(Pmin(a, h) for h in hs)
print(f"[5] P_min satisfies all odd k-block gadgets and random alternating F/T decompositions: {ext_ok}")

# m = 2 'gadget' with the flattening-minimal shift ceil((a-1)/2) (no such tensor gadget is claimed; informative)
m2 = [(a, h) for a in range(1, 60) for h in range(1, 60)
      if Pmin(a, 2 * h + (a) // 2) < 2 * Pmin(a, h)]
print(f"    hypothetical 2-fold gadget P(a,2h+floor(a/2)) >= 2P(a,h) violated by P_min at {len(m2)} pairs, e.g. {m2[:6]}")

# ---------------------------------------------------------------- (6) flattening-rank profile (a real character) : sanity of the extension
fl_ok, fl_dom = True, True
eps = 1e-9
for a in range(1, 80):
    for b in range(2, 300):
        fl_ok &= 2 * Pflat(a, b) >= Pflat(a, b + 1) + Pflat(a, b - 1) - eps
    for h in range(1, 120):
        for k in (3, 5, 7):
            fl_ok &= Pflat(a, k * h + (k - 1) * (a - 1) // 2) >= k * Pflat(a, h) - eps
    for b in range(1, 300):
        fl_dom &= Pflat(a, b) >= float(Pmin(a, b)) - eps
for _ in range(4000):
    a = random.randint(1, 30)
    mF = random.randint(2, 4)
    hs = [random.randint(1, 30) for _ in range(2 * mF - 1)]
    B = sum(hs[0::2]) + sum(h + a - 1 for h in hs[1::2])
    fl_ok &= Pflat(a, B) >= sum(Pflat(a, h) for h in hs) - eps
print(f"[6] FLOAT flattening profile sqrt(ab(a+b-1)) satisfies (C), all k-block and random F/T gadgets: {fl_ok}; dominates P_min: {fl_dom}")
# flattening forbids a smaller shift: P(a, 3h + a - 2) >= 3P(a,h) must fail for large h
bad = [(a, h) for a in (2, 5, 20) for h in (10**3, 10**4, 10**5) if Pflat(a, 3*h + a - 2) < 3*Pflat(a, h)]
print(f"    shift a-2 instead of a-1 is refuted by flattening ranks at (a,h) = {bad}")

# ---------------------------------------------------------------- (7) FLOAT check on a large grid
F = N_FLOAT
lg = [0.0, 0.0]
for a in range(2, 8 * F + 10):
    lg.append(lg[-1] + math.log1p(1.0 / (3 * (a - 1))))
def Pf(a, b):
    m, M = min(a, b), max(a, b)
    return math.exp(lg[m]) * (2 * M + m - 1) / 2
big_ok = True
for a in range(1, F + 1, 7):
    for b in range(a, 4 * F, 11):
        big_ok &= Pf(a, b) <= (a + b - 1) ** (4.0 / 3.0) * (1 + 1e-12)
        if b >= 2:
            big_ok &= 2 * Pf(a, b) >= Pf(a, b + 1) + Pf(a, b - 1) - 1e-9 * Pf(a, b)
    for h in range(1, F, 13):
        big_ok &= Pf(a, 3 * h + a - 1) >= 3 * Pf(a, h) * (1 - 1e-12)
print(f"[7] FLOAT grid a<={F}: rank bound at t=3/4, concavity, tripling: {big_ok}")

# ---------------------------------------------------------------- (8) asymptotics and finite-size exponents
G43 = math.gamma(4.0 / 3.0)
print(f"[8] D_a^min ~ (3/(2 Gamma(4/3))) a^(4/3) = {3/(2*G43):.6f} a^(4/3);  rank bound at t=3/4 on diagonal: 2^(4/3) a^(4/3) = {2**(4/3):.6f} a^(4/3)")
def logD(a):  # log D_a^min via lgamma
    return math.log((3 * a - 1) / 2) + math.lgamma(a + 1/3) - math.lgamma(4/3) - math.lgamma(a)
def omega_a(a):
    return 3 * math.log(2 * a - 1) / logD(a)
print("    omega_a = 3 log(2a-1)/log D_a^min (each is a valid bound omega <= omega_a, using only C(a',b'), a'<=a, b'<=4a-1):")
for a in (2, 3, 4, 5, 10, 20, 50, 100, 1000, 10**4, 10**6, 10**9, 10**12):
    print(f"      a = {a:>14}: omega_a = {omega_a(a):.6f}   [paper's crude a^(4/3): {3*math.log(2*a-1)/((4/3)*math.log(a)) if a>1 else float('nan'):.6f}]")
print(f"    exact: omega_2 = 3 log 3 / log(10/3) = {3*math.log(3)/math.log(10/3):.6f} < log2 7 = {math.log2(7):.6f}")
import mpmath as mp
mp.mp.dps = 300
def om_mp(a):
    a = mp.mpf(a)
    lD = mp.log((3*a-1)/2) + mp.loggamma(a+mp.mpf(1)/3) - mp.loggamma(mp.mpf(4)/3) - mp.loggamma(a)
    return 3*mp.log(2*a-1)/lD
def first_a(thr):
    thr = mp.mpf(thr)
    lo, hi = 2, 2
    while om_mp(hi) >= thr:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if om_mp(mid) < thr:
            hi = mid
        else:
            lo = mid + 1
    return lo
print("    (thresholds below use mpmath at 300 digits; float lgamma differences are unreliable beyond a ~ 1e12)")
for name, thr in [("Strassen log2 7", mp.log(7, 2)), ("2.5", "2.5"), ("Coppersmith-Winograd 2.375477", "2.375477"),
                  ("Alman et al. 2025: 2.371339", "2.371339"), ("Dupont et al. 2026: 2.371177", "2.371177"),
                  ("companion every-field 2.371054886006746", "2.371054886006746"),
                  ("2.30", "2.30"), ("2.28", "2.28"), ("2.26", "2.26"), ("companion C: 2.258", "2.258")]:
    a = first_a(thr)
    print(f"    first a with omega_a < {name}: a = {a}  (log10 a = {mp.nstr(mp.log10(a), 6)})")
K = mp.mpf(9)/4*(mp.log(2) - mp.mpf(3)/4*mp.log(3/(2*mp.gamma(mp.mpf(4)/3))))
print(f"    omega_a = 9/4 + K/ln a + O(1/ln^2 a), K = (9/4)(ln 2 - (3/4) ln(3/(2 Gamma(4/3)))) = {mp.nstr(K, 10)}")

print("\nALL CHECKS PASSED" if not fails and ok and gf_ok and gf2_ok and rb_ok and ext_ok and fl_ok and fl_dom and big_ok else "\nSOME CHECK FAILED")
