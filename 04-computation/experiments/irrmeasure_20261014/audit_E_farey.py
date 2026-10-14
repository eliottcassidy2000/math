"""Audit E (2026-10-14) of THM-4633 statement 5 / Hercher's Corollary 29.
Rigorous enumeration of ALL reduced fractions K/k with k < q_23 = 137528045312 lying in the admissible window
(theta, theta + 1/(3 M0 ln 2)) for M0 in {2^68, 3*2^69, 2^71}, by the Stern-Brocot 'simplest fraction in an interval'
recursion (every fraction with denominator < D in an open interval is found by repeatedly taking the smallest-denominator
fraction and splitting).  Expected: at 2^71 the window contains exactly one fraction with k < q_23, the semiconvergent
114208327604/72057431991, which the plain inequality therefore does NOT exclude (its budget is 3781*2^60 = 4.36e21 > 2^71);
this is the fraction Hercher excludes in Corollary 29 with his finer (computer-assisted) argument.
Also recomputes the numbers of Hercher's Remark 28: the window width 1.1032e-22 and the thresholds 3781*2^60, 2836*2^60."""
from fractions import Fraction as Fr
from mpmath import mp, mpf, log, floor, nstr
mp.dps = 120
theta = log(3) / log(2); ln2 = log(2); ln3 = log(3)
# rational sandwich theta_lo < theta < theta_hi with 1e-80 precision
D = 10 ** 80
theta_lo = Fr(int(floor(theta * D)), D); theta_hi = theta_lo + Fr(1, D)
assert mpf(theta_lo.numerator) / theta_lo.denominator < theta < mpf(theta_hi.numerator) / theta_hi.denominator

def simplest_between(lo, hi):
    """smallest-denominator fraction strictly between lo < hi (Stern-Brocot), lo, hi Fractions"""
    # standard continued-fraction algorithm
    fl = lo.numerator // lo.denominator
    if fl + 1 < hi:  # an integer strictly inside
        return Fr(fl + 1)
    if lo.denominator == 1 and fl + 1 == hi:  # hmm: open interval, no integer inside, fall through
        pass
    # shift by fl, invert
    lo2, hi2 = lo - fl, hi - fl   # 0 <= lo2 < hi2 <= 1 (hi2 may be > 1 only if integer inside, handled)
    if lo2 == 0:
        # fractions in (0, hi2): smallest denominator is 1/n with 1/n < hi2
        n = hi2.denominator // hi2.numerator + 1
        return fl + Fr(1, n)
    inner = simplest_between(1 / hi2, 1 / lo2)
    return fl + 1 / inner

def all_fractions_below(lo, hi, Dmax):
    """all fractions with denominator < Dmax strictly inside (lo, hi)"""
    out = []
    stack = [(lo, hi)]
    while stack:
        a, b = stack.pop()
        f = simplest_between(a, b)
        if f.denominator >= Dmax: continue
        out.append(f); stack.append((a, f)); stack.append((f, b))
    return sorted(out)

q23 = 137528045312
for name, M0 in (("2^68", mpf(2) ** 68), ("3*2^69", 3 * mpf(2) ** 69), ("2^71", mpf(2) ** 71), ("3781*2^60", 3781 * mpf(2) ** 60), ("2^72", mpf(2) ** 72)):
    eps = 1 / (3 * M0 * ln2)
    eps_r = Fr(int(floor(eps * D)) + 1, D)  # slightly enlarged rational window
    found = all_fractions_below(theta_lo, theta_hi + eps_r, q23)
    print("M0 = %s: window width eps = %s; fractions with k < q_23 in the (slightly enlarged) window:" % (name, nstr(eps, 6)))
    for f in found:
        d = f.numerator * ln2 - f.denominator * ln3
        inside = (d > 0) and (mpf(f.numerator) / f.denominator - theta < eps)
        print("   K/k = %d/%d  K/k - theta = %s  budget m_max = %s  strictly inside window: %s" % (
            f.numerator, f.denominator, nstr(mpf(f.numerator) / f.denominator - theta, 6), nstr(f.denominator / (3 * d), 6) if d > 0 else "-", inside))
    if not found: print("   (none)")

# Hercher's Remark 28 numbers for the semiconvergent s = 114208327604/72057431991
P, Q = 114208327604, 72057431991
dist = mpf(P) / Q - theta
print("\nsemiconvergent s: (K+L)/K - delta in Hercher's notation = K/k - theta =", nstr(dist, 6), " (Hercher Remark 28: 1.1032e-22)")
bud = Q / (3 * (P * ln2 - Q * ln3))
print("plain budget m_max(s) =", nstr(bud, 8), "= %s * 2^60  (Hercher Remark 28: 'previous methods' need X_0 >= 3781 * 2^60)" % nstr(bud / mpf(2) ** 60, 6))
print("with Hercher's Theorem 27 factor 3/4: (3/4) m_max(s) =", nstr(0.75 * bud / mpf(2) ** 60, 6), "* 2^60  (Remark 28: 2836 * 2^60)")
print("Hercher Corollary 29 hypothesis 3*2^69 = 1536*2^60 =", nstr(3 * mpf(2) ** 69, 6), "; Barina 2025 limit 2^71 = 2048*2^60 =", nstr(mpf(2) ** 71, 6), "; 2^71 >= 3*2^69:", 2 ** 71 >= 3 * 2 ** 69)
print("Hercher's unconditional input in the paper: X_0 = 704*2^60 =", nstr(704 * mpf(2) ** 60, 6), "(= 2^%s)" % nstr(log(704 * mpf(2) ** 60) / ln2, 5))
