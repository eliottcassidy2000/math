#!/usr/bin/env python3
"""
Independent audit of THM-4476 (thin divergence) and of the note
05-knowledge/results/collatz_thin_20260925_thin_divergent_orbits.md.

Written WITHOUT reading collatz_thin_20260925_counts.py.

Sections
  A  constants: h*, 1/h*, theta_0, thresholds, un-bootstrapped exponent
  B  the entropy step of the counting lemma (exact binomial sums vs bound)
  C  the no-dip sets F_b(X, theta) on both sheets, counts vs the lemma bound
  D  Proposition 6 identities with exact fractions (minus and plus sheet),
     the factor 1 - 1/(3 m_l), rational examples (vanishing / negative factor)
  E  the dichotomy / pigeonhole of section 1.5 on real finite segments,
     sign-change count of section 1.1, Terras bijection of section 1.2
  F  the withdrawn Corollary 4 (sign strategies): level-2 strategies where
     the parity map is NOT a bijection and where the counting lemma FAILS
  G  the replacement corollaries: injectivity of T_b on periodic points
     (new Cor. 4), the log(1-x) inequality of Cor. 6, the exponent ladder
     a(d) = max(1,a), and "bijection only for a constant shift" at level 3
"""
import math
import sys
from fractions import Fraction
from math import comb, log2

try:
    sys.stdout.reconfigure(newline="\n")   # raw LF output (hash_basis: raw LF bytes)
except Exception:
    pass

ALPHA = log2(3.0)                      # log_2 3
RHO0 = 1.0 / ALPHA                     # log_3 2


def h(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p)


def rho_of(theta):
    return (1.0 - theta) / ALPHA


def k0_of(theta):
    """least k with rho*k - 2 > k/2, i.e. k > 2/(rho - 1/2)."""
    r = rho_of(theta)
    return math.floor(2.0 / (r - 0.5)) + 1


def lemma_bound(X, theta, absb=1):
    r = rho_of(theta)
    A = 2.0 * (r / (1.0 - r)) ** 2
    return (2.0 * absb * X ** (log2(1.5) + theta)
            + A * (log2(X) + 1.0) * X ** h(r))


def T(x, b):
    return x // 2 if x % 2 == 0 else (3 * x + b) // 2


# ----------------------------------------------------------------- A
print("=" * 72)
print("A. constants")
print("=" * 72)
hstar = h(RHO0)
print("log_2 3      = %.9f" % ALPHA)
print("log_3 2      = %.9f" % RHO0)
print("h*           = h(log_3 2) = %.9f" % hstar)
print("1/h*         = %.9f" % (1.0 / hstar))
print("theta_0      = 1 - log_2(3)/2 = %.9f" % (1.0 - ALPHA / 2.0))
print("(1/h*-1)/2   = %.9f   (plus-sheet log-band threshold, Cor. 3)" % ((1.0 / hstar - 1.0) / 2.0))
print("1/h*         = %.9f   (minus-sheet log-band threshold, Cor. 3)" % (1.0 / hstar))
print("log_2(3/2)   = %.9f" % log2(1.5))
for th in (0.03, 0.1):
    r = rho_of(th)
    print("theta=%.2f: rho=%.6f  h(rho)=%.6f  k_0(theta)=%d  A(theta)=2(rho/(1-rho))^2=%.4f  log2(3/2)+theta=%.4f"
          % (th, r, h(r), k0_of(th), 2.0 * (r / (1.0 - r)) ** 2, log2(1.5) + th))
# un-bootstrapped exponent: 1 - theta = h(rho(theta))
lo, hi = 1e-6, 0.2
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if (1.0 - mid) - h(rho_of(mid)) > 0:
        lo = mid
    else:
        hi = mid
print("root of 1-theta = h(rho(theta)): theta=%.6f, exponent 1-theta=%.6f" % (lo, 1.0 - lo))
print("monotonicity check: h decreasing on [1/2,1]:",
      all(h(0.5 + i / 2000.0) > h(0.5 + (i + 1) / 2000.0) for i in range(999)))
print("|h'(p)| = |log2((1-p)/p)| <= log2(rho/(1-rho)) on [1/2,rho] for rho=rho(0.03):",
      all(abs(log2((1 - p) / p)) <= log2(rho_of(0.03) / (1 - rho_of(0.03))) + 1e-12
          for p in [0.5 + (rho_of(0.03) - 0.5) * i / 1000.0 for i in range(1001)]))

# ----------------------------------------------------------------- B
print()
print("=" * 72)
print("B. entropy step: S(k) = sum_{o >= rho k - 2} C(k,o) vs bounds")
print("   bound1 = (k+1) 2^{k h(rho - 2/k)} (valid once rho k - 2 > k/2)")
print("   bound2 = (k+1) (rho/(1-rho))^2 2^{k h(rho)}  (what the lemma uses)")
print("   also: C(k,o) <= 2^{k h(o/k)} for all o (standard)")
print("=" * 72)
ok_std = True
for k in range(1, 401):
    for o in range(0, k + 1):
        if comb(k, o) > 2.0 ** (k * h(o / k)) * (1 + 1e-9):
            ok_std = False
print("C(k,o) <= 2^{k h(o/k)} for all k<=400, o:", ok_std)
for th in (0.03, 0.1):
    r = rho_of(th)
    k0 = k0_of(th)
    worst1 = worst2 = worst_note = 0.0
    fail1 = fail2 = []
    fail1 = []
    fail2 = []
    fail_note = []
    for k in range(1, 401):
        omin = max(0, math.ceil(r * k - 2 - 1e-12))
        S = sum(comb(k, o) for o in range(omin, k + 1))
        b2 = (k + 1) * (r / (1 - r)) ** 2 * 2.0 ** (k * h(r))
        note_bound = (k + 1) * 2.0 ** (k * h(r))        # the note's "(L1)" claim: ratio below k+1
        if k >= k0:
            b1 = (k + 1) * 2.0 ** (k * h(r - 2.0 / k))
            worst1 = max(worst1, S / b1)
            worst2 = max(worst2, S / b2)
            if S > b1:
                fail1.append(k)
            if S > b2:
                fail2.append(k)
        if S > note_bound:
            fail_note.append(k)
        worst_note = max(worst_note, S / note_bound)
    print("theta=%.2f (rho=%.5f, k_0=%d): max S/bound1 over k in [k_0,400] = %.4f, failures %s"
          % (th, r, k0, worst1, fail1[:5]))
    print("            max S/bound2 = %.4f, failures %s" % (worst2, fail2[:5]))
    print("            note's (L1) claim S <= (k+1) 2^{k h(rho)} for ALL k<=400: max ratio S/((k+1)2^{kh}) = %.4f, failures %s"
          % (worst_note, fail_note[:8]))

# ----------------------------------------------------------------- C
print()
print("=" * 72)
print("C. no-dip sets F_b(X,theta) = {y in [1,X]: T_b^i(y) >= y X^{-theta}, 0<=i<=k}, k=floor(log2 X)")
print("=" * 72)
XS = [2 ** 10, 2 ** 12, 2 ** 14, 2 ** 16, 2 ** 18, 2 ** 20]
THETAS = [0.03, 0.1]


def count_F(X, thetas, b):
    """returns dict theta -> (count, count_minus, count_plus) where the +- variants
    perturb the threshold by a relative 1e-9 to expose floating-point fragility."""
    k = int(math.floor(log2(X) + 1e-12))
    thr = {th: X ** (-th) for th in thetas}
    thr_lo = {th: thr[th] * (1 - 1e-9) for th in thetas}
    thr_hi = {th: thr[th] * (1 + 1e-9) for th in thetas}
    tmin = min(thr_lo.values())
    cnt = {th: [0, 0, 0] for th in thetas}
    for y in range(1, X + 1):
        x = y
        m = 1.0
        for _ in range(k):
            x = x // 2 if x % 2 == 0 else (3 * x + b) // 2
            rr = x / y
            if rr < m:
                m = rr
                if m < tmin:
                    break
        for th in thetas:
            if m >= thr[th]:
                cnt[th][0] += 1
            if m >= thr_lo[th]:
                cnt[th][1] += 1
            if m >= thr_hi[th]:
                cnt[th][2] += 1
    return cnt


results = {}
for b in (1, -1):
    for X in XS:
        c = count_F(X, THETAS, b)
        for th in THETAS:
            results[(b, X, th)] = c[th]
for th in THETAS:
    r = rho_of(th)
    print("theta=%.2f  rho=%.5f  h(rho)=%.5f  k_0=%d (lemma asserted only for X >= 2^%d)" % (th, r, h(r), k0_of(th), k0_of(th)))
    print("  %-8s %-6s %-9s %-9s %-14s %-9s %-9s %-9s" % ("X", "k", "#F_{+1}", "#F_{-1}", "bound", "ratio+", "logF/logX", "slope+"))
    prev = None
    for X in XS:
        k = int(log2(X))
        cp = results[(1, X, th)]
        cm = results[(-1, X, th)]
        B = lemma_bound(X, th, 1)
        slope = ""
        if prev is not None:
            slope = "%.4f" % (log2(cp[0] / prev) / 2.0)
        boundary = ""
        if cp[1] != cp[0] or cp[2] != cp[0] or cm[1] != cm[0] or cm[2] != cm[0]:
            boundary = "  [X^(-theta) = %s is exact: strict-inequality counts %d / %d]" % (
                Fraction(X ** (-th)).limit_denominator(1 << 30), cp[2], cm[2])
        print("  2^%-5d %-6d %-9d %-9d %-14.1f %-9.5f %-9.4f %-9s%s"
              % (k, k, cp[0], cm[0], B, cp[0] / B, log2(cp[0]) / k, slope, boundary))
        prev = cp[0]
    print("  sheet difference |#F_{+1} - #F_{-1}| max over X:",
          max(abs(results[(1, X, th)][0] - results[(-1, X, th)][0]) for X in XS))
    # least-squares slope of log2 #F against log2 X over 2^14..2^20
    pts = [(log2(X), log2(results[(1, X, th)][0])) for X in XS if X >= 2 ** 14]
    mx = sum(p[0] for p in pts) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    lsq = sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)
    print("  lemma exponent h(rho)=%.4f ; the note's 'growth exponent' is log F/log X at 2^20 = %.4f ; least-squares slope 2^14..2^20 = %.4f ; local slope 2^18->2^20 = %.4f"
          % (h(r), log2(results[(1, 2 ** 20, th)][0]) / 20.0, lsq,
             log2(results[(1, 2 ** 20, th)][0] / results[(1, 2 ** 18, th)][0]) / 2.0))

# ----------------------------------------------------------------- D
print()
print("=" * 72)
print("D. Proposition 6 identities with exact fractions")
print("=" * 72)


def syracuse(n, b, L):
    """odd iterates m_0=n, m_1, ..., m_L and cumulative halvings d_0=0, d_1, ..., d_L for q=3."""
    ms = [Fraction(n)]
    ds = [0]
    m = Fraction(n)
    d = 0
    for _ in range(L):
        t = 3 * m + b
        # 2-adic valuation of the numerator (denominator odd)
        num = t.numerator
        if num == 0:
            return None
        v = 0
        while num % 2 == 0:
            num //= 2
            v += 1
        m = t / (2 ** v)
        d += v
        ms.append(m)
        ds.append(d)
    return ms, ds


LMAX = 60
bad_minus = 0
bad_plus = 0
factor_min = Fraction(2)
factor_max = Fraction(0)
checked = 0
for n in range(1, 301, 2):
    # minus sheet: m_L 2^{d_L} / 3^L = n prod_{l<L} (1 - 1/(3 m_l))
    ms, ds = syracuse(n, -1, LMAX)
    prod = Fraction(1)
    for L in range(0, LMAX + 1):
        lhs = ms[L] * Fraction(2) ** ds[L] / Fraction(3) ** L
        rhs = n * prod
        if lhs != rhs:
            bad_minus += 1
        # also R_L(d) = n (1 - prod) with R_L = sum_{l<L} 2^{d_l}/3^{l+1}
        RL = sum(Fraction(2) ** ds[l] / Fraction(3) ** (l + 1) for l in range(L))
        if RL != n * (1 - prod):
            bad_minus += 1
        if L < LMAX:
            f = 1 - 1 / (3 * ms[L])
            factor_min = min(factor_min, f)
            factor_max = max(factor_max, f)
            prod *= f
        checked += 1
    # plus sheet: m_j 2^{d_j}/3^j = n + (1/3) sum_{i<j} 2^{d_i}/3^i   (Cor. 3's identity)
    ms, ds = syracuse(n, +1, LMAX)
    for j in range(0, LMAX + 1):
        lhs = ms[j] * Fraction(2) ** ds[j] / Fraction(3) ** j
        rhs = n + Fraction(1, 3) * sum(Fraction(2) ** ds[i] / Fraction(3) ** i for i in range(j))
        if lhs != rhs:
            bad_plus += 1
        if j > 0:
            # 2^{-Delta_j} = 3^j / 2^{d_j}
            pass
print("minus sheet (3n-1), odd n<=300, L<=%d: identity m_L 2^{d_L}/3^L = n prod_{l<L}(1-1/(3m_l)) and R_L(d)=n(1-prod): failures = %d of %d checks"
      % (LMAX, bad_minus, checked))
print("  every factor 1-1/(3 m_l) in [%s, %s] = [%.6f, %.6f] (all in (0,1))" % (factor_min, factor_max, float(factor_min), float(factor_max)))
print("plus sheet (3n+1), odd n<=300, j<=%d: identity m_j 2^{d_j}/3^j = n + (1/3) sum_{i<j} 2^{d_i}/3^i: failures = %d" % (LMAX, bad_plus))
# Prop 6(6) 2-adic statement: 3^L n - B_L (minus sheet) and 3^L n + B_L (plus sheet) divisible by 2^{d_L}
bad2 = 0
for n in range(1, 301, 2):
    for b in (-1, 1):
        ms, ds = syracuse(n, b, 40)
        BL = sum(3 ** (40 - 1 - l) * 2 ** ds[l] for l in range(40))
        if (3 ** 40 * n + b * BL) % (2 ** ds[40]) != 0:
            bad2 += 1
print("Prop 6(6): 3^L n + b B_L = 0 mod 2^{d_L} (L=40, both sheets, odd n<=300): failures = %d" % bad2)
# c_L -> 0 for eventually periodic minus-sheet orbits (all n<=300 enter a cycle); floats via logs
def c_log2(n, L):
    m = n
    d = 0
    for _ in range(L):
        t = 3 * m - 1
        v = (t & -t).bit_length() - 1
        m = t >> v
        d += v
    return log2(m) + d - L * ALPHA


for L in (400, 3000, 20000):
    worst = max(c_log2(n, L) for n in range(1, 301, 2))
    print("minus sheet, odd n<=300: max_n log2 c_L at L=%d is %.2f  (c_L = m_L 2^{d_L}/3^L -> 0 since every orbit enters a cycle)" % (L, worst))
print("  (the slow case is the 3n-1 cycle 17,25,37,55,41,61,91: factor per period = %.6f)" %
      float(Fraction(1) * (1 - Fraction(1, 51)) * (1 - Fraction(1, 75)) * (1 - Fraction(1, 111)) * (1 - Fraction(1, 165)) * (1 - Fraction(1, 123)) * (1 - Fraction(1, 183)) * (1 - Fraction(1, 273))))
# rational examples for Corollary 2's factor
print("rational examples for the factor 1 - 1/(3m):")
for m in (Fraction(1, 3), Fraction(1, 5), Fraction(-1, 5), Fraction(3, 11), Fraction(-3, 7), Fraction(1), Fraction(-1)):
    print("   m=%6s : 1-1/(3m) = %s" % (m, 1 - 1 / (3 * m)))
print("  T_{-1} orbit of 1/3: (3*(1/3)-1)/2 = 0 -> word ends in 0s, not a halving word (factor 0 excluded)")
ms, ds = syracuse(Fraction(1, 5), -1, 6)
print("  T_{-1} odd iterates of 1/5:", [str(m) for m in ms], " (eventually periodic; a negative factor occurs)")

# ----------------------------------------------------------------- E
print()
print("=" * 72)
print("E. section 1.1 sign changes, 1.2 Terras bijection, 1.5 dichotomy/pigeonhole on real segments")
print("=" * 72)


def orbit_until_repeat(x, b, cap=100000):
    seen = set()
    orb = []
    while x not in seen and x != 0 and len(orb) < cap:
        seen.add(x)
        orb.append(x)
        x = T(x, b)
    return orb


# 1.1: sign changes
viol = 0
maxchg = {}
for b in range(-31, 32, 2):
    for x0 in range(-500, 501):
        if x0 == 0:
            continue
        orb = orbit_until_repeat(x0, b)
        chg = 0
        for i in range(len(orb) - 1):
            if orb[i] * orb[i + 1] < 0:
                chg += 1
                if not (orb[i] % 2 != 0 and abs(orb[i]) <= abs(b) / 3):
                    viol += 1
        maxchg[b] = max(maxchg.get(b, 0), chg)
print("sign changes only at odd |x|<=|b|/3 (b in -31..31 odd, |x0|<=500, until repeat): violations =", viol)
print("max #sign changes vs |b|/3+1:", ["b=%d:%d<=%.2f" % (b, maxchg[b], abs(b) / 3 + 1) for b in sorted(maxchg) if maxchg[b] > 0][:12])
print("bound |S| <= |b|/3+1 respected for all b:", all(maxchg[b] <= abs(b) / 3 + 1 for b in maxchg))

# 1.2: Terras bijection for constant odd b
def word(y, b, k):
    w = []
    for _ in range(k):
        w.append(y & 1)
        y = T(y, b)
    return tuple(w)


allbij = True
for b in (1, -1, 3, -3, 5, 7, 9, 11, -15, 21):
    for k in range(1, 13):
        words = {word(y, b, k) for y in range(2 ** k)}
        if len(words) != 2 ** k:
            allbij = False
            print("  NOT a bijection: b=%d k=%d distinct words=%d" % (b, k, len(words)))
print("Terras bijection Z/2^k -> {0,1}^k for constant odd b (b in a sample, k<=12):", allbij)

# 1.3 carry bound check: |beta_k| <= |b| ((3/2)^k - 1) exhaustively over words (via residues)
worst_ratio = 0.0
for b in (1, -1, 5, -7):
    k = 12
    for y in range(2 ** k):
        x = y
        o = 0
        for _ in range(k):
            if x & 1:
                o += 1
            x = T(x, b)
        beta = x - Fraction(3 ** o * y, 2 ** k)
        worst_ratio = max(worst_ratio, float(abs(beta)) / (abs(b) * ((1.5) ** k - 1)))
print("carry bound: max |beta_k| / (|b|((3/2)^k-1)) over residues mod 2^12, b in {1,-1,5,-7}: %.4f (<= 1)" % worst_ratio)


# 1.5: dichotomy on real finite segments
def classify(seg, X, theta):
    """seg: positive segment y_0..y_M (distinct). returns (N, nE, nD, nND, k, N_low, F)"""
    M = len(seg) - 1
    k = int(math.floor(log2(X) + 1e-12))
    thr = X ** (-theta)
    nE = nD = nND = 0
    landing = {}
    for i, y in enumerate(seg):
        if y > X:
            continue
        if i > M - k:
            nE += 1
            continue
        dip = None
        for s in range(1, k + 1):
            if seg[i + s] < y * thr:
                dip = s
                break
        if dip is None:
            nND += 1
        else:
            nD += 1
            landing[i + dip] = landing.get(i + dip, 0) + 1
            assert seg[i + dip] < X ** (1 - theta)
    N = sum(1 for y in seg if y <= X)
    Nlow = sum(1 for y in seg if y <= X ** (1 - theta))
    maxserve = max(landing.values()) if landing else 0
    return N, nE, nD, nND, k, Nlow, maxserve


def F_count_direct(X, theta, b):
    k = int(math.floor(log2(X) + 1e-12))
    thr = X ** (-theta)
    c = 0
    for y in range(1, X + 1):
        x = y
        ok = True
        for _ in range(k):
            x = T(x, b)
            if x < y * thr:
                ok = False
                break
        if ok:
            c += 1
    return c


# longest positive transients
def longest_transient(b, lo, hi):
    best = (0, None)
    for n in range(lo, hi):
        L = len(orbit_until_repeat(n, b))
        if L > best[0]:
            best = (L, n)
    return best


segs = []
segs.append((1, orbit_until_repeat(27, 1)))
L, n = longest_transient(1, 1, 20000)
segs.append((1, orbit_until_repeat(n, 1)))
L, n = longest_transient(-1, 1, 20000)
segs.append((-1, orbit_until_repeat(n, -1)))
L, n = longest_transient(5, 1, 20000)
segs.append((5, orbit_until_repeat(n, 5)))
L, n = longest_transient(-7, 1, 20000)
segs.append((-7, orbit_until_repeat(n, -7)))
allok = True
for b, orb in segs:
    # split into maximal one-signed segments; use the positive ones as T_b segments,
    # negative ones as T_{-b} segments of -x
    pieces = []
    cur = [orb[0]]
    for x in orb[1:]:
        if x * cur[-1] > 0:
            cur.append(x)
        else:
            pieces.append(cur)
            cur = [x]
    pieces.append(cur)
    for piece in pieces:
        if piece[0] > 0:
            bb, seg = b, piece
        else:
            bb, seg = -b, [-x for x in piece]
        if len(seg) < 8:
            continue
        for X in (2 ** 6, 2 ** 8, 2 ** 10, 2 ** 12):
            for th in (0.03, 0.1):
                N, nE, nD, nND, k, Nlow, maxserve = classify(seg, X, th)
                F = F_count_direct(X, th, bb)
                okD = nD <= k * Nlow
                okND = nND <= F
                okpart = (N == nE + nD + nND)
                okE = nE <= k
                okserve = maxserve <= k
                if not (okD and okND and okpart and okE and okserve):
                    allok = False
                    print("  VIOLATION b=%d start=%d X=%d theta=%.2f: N=%d E=%d D=%d ND=%d k=%d N(X^(1-th))=%d F=%d maxserve=%d"
                          % (bb, seg[0], X, th, N, nE, nD, nND, k, Nlow, F, maxserve))
        N, nE, nD, nND, k, Nlow, maxserve = classify(seg, 2 ** 12, 0.1)
        print("  b=%3d segment start=%7d len=%4d  X=2^12 th=0.1: N=%d = E %d + D %d + ND %d ; k N(X^(1-th)) = %d ; #F=%d ; max dippers per landing = %d"
              % (bb, seg[0], len(seg), N, nE, nD, nND, k * Nlow, F_count_direct(2 ** 12, 0.1, bb), maxserve))
print("partition N=E+D+ND, E<=k, D<=k N(X^(1-theta)), ND<=#F, <=k dippers per landing point, on all tested segments:", allok)
print("note: #(E) equals k exactly on every segment above (indices M-k+1..M), i.e. 'at most k', not 'fewer than k'; (R) uses k, so harmless.")

# ----------------------------------------------------------------- F
print()
print("=" * 72)
print("F. Corollary 4: level-2 sign strategies (THM-4474): sigma(1 mod 4), sigma(3 mod 4) in {+1,-1}")
print("=" * 72)


def Tsig(x, s1, s3):
    if x % 2 == 0:
        return x // 2
    return (3 * x + (s1 if x % 4 == 1 else s3)) // 2


def word_sig(y, s1, s3, k):
    w = []
    for _ in range(k):
        w.append(y & 1)
        y = Tsig(y, s1, s3)
    return tuple(w)


for (s1, s3) in ((1, 1), (-1, -1), (1, -1), (-1, 1)):
    print("strategy sigma(1)=%+d sigma(3)=%+d:" % (s1, s3))
    for k in (2, 3, 4, 8, 12):
        fib = {}
        for y in range(2 ** k):
            w = word_sig(y, s1, s3, k)
            fib[w] = fib.get(w, 0) + 1
        maxf = max(fib.values())
        print("   k=%2d: distinct words %5d of %5d, max fibre (classes mod 2^k per word) %d, word (1,1,...) attained: %s"
              % (k, len(fib), 2 ** k, maxf, (tuple([1] * k) in fib)))
    # parity of T(x) for odd x
    par = {x % 4: Tsig(x, s1, s3) % 2 for x in (1, 3, 5, 7, 9, 11)}
    print("   parity of T(x) for odd x: x=1 mod 4 -> %d, x=3 mod 4 -> %d" % (Tsig(1, s1, s3) % 2, Tsig(3, s1, s3) % 2),
          "(consistent over 1,3,5,7,9,11:", len({(x % 4, Tsig(x, s1, s3) % 2) for x in (1, 3, 5, 7, 9, 11)}) == 2, ")")

# counting lemma for the all-odd strategy sigma(1)=-1, sigma(3)=+1
print("all-odd strategy sigma(1)=-1, sigma(3)=+1: T(x) odd and >= x for every odd x >= 1, so every odd y<=X is in F(X,theta).")
for X in (2 ** 12, 2 ** 16, 2 ** 20):
    k = int(log2(X))
    for th in (0.03, 0.1):
        thr = X ** (-th)
        cnt = 0
        for y in range(1, X + 1):
            x = y
            ok = True
            for _ in range(k):
                x = Tsig(x, -1, 1)
                if x < y * thr:
                    ok = False
                    break
            if ok:
                cnt += 1
        print("   X=2^%d theta=%.2f: #F = %d (= X/2 + %d), lemma bound with |b|=1: %.3e, ratio #F/X = %.3f (bound/X = %.3f -> 0)"
              % (k, th, cnt, cnt - X // 2, lemma_bound(X, th, 1), cnt / X, lemma_bound(X, th, 1) / X))
for th in (0.03, 0.1):
    k = 1
    while 2.0 ** (k - 1) <= lemma_bound(2.0 ** k, th, 1):
        k += 1
    print("   theta=%.2f: X/2 exceeds the lemma bound from k = %d on (X = 2^%d): the counting lemma is FALSE for this strategy for all X >= 2^%d"
          % (th, k, k, k))
print("all orbits of that strategy are eventually increasing odd sequences (x -> (3x+-1)/2), so N(X)=O(log X): the THEOREM's conclusion holds there trivially,")
print("but sections 1.2 and 1.4 do not survive, contrary to the proof of the (withdrawn) Corollary 4.")

# ----------------------------------------------------------------- G
print()
print("=" * 72)
print("G. replacement corollaries")
print("=" * 72)

# G1: new Corollary 4 -- T_b is injective on the set of periodic points (unique periodic preimage)
print("G1. periodic points of T_b with |p| <= X0: is T_b injective on them?  (search from every |x| <= X0, step cap 20000)")
for b in (1, -1, 5, -5, 7, -7, 11, 13, -13, 17, 23):
    X0 = 20000
    periodic = set()
    for x0 in range(-X0, X0 + 1):
        if x0 == 0:
            continue
        seen = {}
        x = x0
        steps = 0
        while x not in seen and x != 0 and steps < 20000 and abs(x) < 10 ** 15:
            seen[x] = steps
            x = T(x, b)
            steps += 1
        if x in seen:
            # the cycle is the tail from seen[x]
            cyc = [y for y, i in seen.items() if i >= seen[x]]
            periodic.update(cyc)
    preimages = {}
    for p in periodic:
        q = T(p, b)
        assert q in periodic
        preimages[q] = preimages.get(q, 0) + 1
    inj = max(preimages.values()) == 1 if preimages else True
    ncyc = 0
    rest = set(periodic)
    while rest:
        p = next(iter(rest))
        x = p
        while True:
            rest.discard(x)
            x = T(x, b)
            if x == p:
                break
        ncyc += 1
    small = sorted(q for q in periodic if abs(q) <= X0)
    print("   b=%3d: %2d cycles, %4d periodic points (|p|<=%d: %d, largest |p| = %d); T_b injective on them: %s; count(|p|<=X0)/X0^h* = %.4f"
          % (b, ncyc, len(periodic), X0, len(small), max(abs(q) for q in periodic) if periodic else 0, inj,
             len(small) / X0 ** hstar))

# G2: Corollary 6's elementary inequality log(1-x) >= -x - x^2 on (0, 1/3]
worst = min(math.log(1 - x) + x + x * x for x in [i / 30000.0 for i in range(1, 10001)])
print("G2. min over x in (0,1/3] of log(1-x) + x + x^2 = %.6f (>= 0 required by Cor. 6): %s" % (worst, worst >= 0))

# G3: exponent ladder -- formal plus-sheet orbit m_j = 2^{-Delta_j} (n + (1/3) sum_{i<j} 2^{Delta_i}) with Delta_j = -a log2 j
print("G3. exponent ladder: Delta_j = -a log2 j, m_j(n=1) = 2^(-Delta_j)(1 + (1/3) sum_{i<j} 2^(Delta_i)); log m_j / log j at j = 10^6 vs max(1,a)")
for a in (0.0, 0.5, 1.0, 1.03, 1.05, 1.5, 2.0):
    s = 1.0  # i = 0 term: Delta_0 = 0
    J = 10 ** 6
    for i in range(1, J):
        s += i ** (-a)
    mj = J ** a * (1.0 + s / 3.0)
    print("   a=%.2f: log m_j/log j = %.4f  (max(1,a) = %.2f)" % (a, math.log(mj) / math.log(J), max(1.0, a)))

# G4: level-3 strategies: the parity word of length K is a function of y mod 2^(K+1); count distinct words and fibres
print("G4. all 16 level-3 strategies sigma:{1,3,5,7}->{+-1}; words of length K=10 over y < 2^(K+1)")


def Tsig3(x, sig):
    if x % 2 == 0:
        return x // 2
    return (3 * x + sig[x % 8]) // 2


K = 10
for mask in range(16):
    sig = {1: 1 if mask & 1 else -1, 3: 1 if mask & 2 else -1, 5: 1 if mask & 4 else -1, 7: 1 if mask & 8 else -1}
    fib = {}
    fn_of_mod2K = {}
    welldef = True
    for y in range(2 ** (K + 1)):
        x = y
        w = []
        for _ in range(K):
            w.append(x & 1)
            x = Tsig3(x, sig)
        w = tuple(w)
        fib[w] = fib.get(w, 0) + 1
        r = y % (2 ** K)
        if r in fn_of_mod2K and fn_of_mod2K[r] != w:
            welldef = False
        fn_of_mod2K[r] = w
    const = len(set(sig.values())) == 1
    print("   sigma=(%+d,%+d,%+d,%+d)%s: distinct words %4d of %4d, max fibre %4d (2 = bijective), word determined by y mod 2^K: %s"
          % (sig[1], sig[3], sig[5], sig[7], " CONST" if const else "      ", len(fib), 2 ** K, max(fib.values()), welldef))
print("   => the parity word map is a bijection exactly for the two constant shifts; every non-constant level-3 strategy loses words.")
print()
print("done.")
