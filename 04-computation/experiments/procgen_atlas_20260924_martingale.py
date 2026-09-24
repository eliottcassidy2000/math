#!/usr/bin/env python3
"""
procgen_atlas_20260924_martingale.py

Implication-atlas lane (session collatz-procgen-20260922, 2026-09-24):
the "fair map" / deterministic-martingale nodes.

Sections (printed in this order):
  F1  pair-sum classification for f_{q,r}(n) = n/2 (n even), (qn+r)/2 (n odd), q,r odd:
      which maps preserve the sum of every pair of a perfect matching of Z (resp. of the
      positive integers) into pairs; the consecutive pairings {2i-1,2i}, {2i,2i+1}.
  F2  arithmetic-mean and geometric-mean step factors of qx+1 and of Collatz's 1932
      permutation; the Lundberg exponent theta(q) of the log-walk; theta(3) = 1 exactly.
  F3  the Haar martingale M_j = 3^{a_j}/2^j on Z_2: exact P_k(lam) = Haar{max_{j<=k} M_j >= lam}
      (exact rationals for k <= 60, floats with a Doob-controlled truncation for k = 4000),
      Doob's bound P <= 1/lam, and the Cramer-Lundberg product lam * P.
  F4  Doob's convergence theorem => Terras: Haar{tau > k} -> 0 (exact values).
  F5  deterministic check: on n in [2^40, 2^40 + 2^20) the fraction of n with
      max_{j<=20} T^j(n) >= lam*n equals P_20(lam) exactly; contrast with n <= 2^20.
  F6  the martingale is exact only in the leading coefficient: sum_{n<=X} T^k(n) - sum n.

Everything in F1, F4, F5, F6 and the k <= 60 part of F3 is exact integer/rational arithmetic.
"""
import math
from fractions import Fraction as Fr
import sys

def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2

# ---------------------------------------------------------------- F1
def F1():
    hdr("F1  pair-sum classification of f_{q,r}(n) = n/2 (even), (qn+r)/2 (odd); q,r odd")
    def f(q, r, n):
        return n // 2 if n % 2 == 0 else (q * n + r) // 2
    # (a) consecutive pairings: P1 = {2i-1, 2i}, P0 = {2i, 2i+1}; exhaustive over |q|,|r| <= 51
    hits = []
    for q in range(-51, 52, 2):
        for r in range(-51, 52, 2):
            ok1 = all(f(q, r, 2*i-1) + f(q, r, 2*i) == 4*i - 1 for i in range(-300, 301))
            ok0 = all(f(q, r, 2*i) + f(q, r, 2*i+1) == 4*i + 1 for i in range(-300, 301))
            if ok1: hits.append(("{2i-1,2i}", q, r))
            if ok0: hits.append(("{2i,2i+1}", q, r))
    print("consecutive pairings preserving every pair sum (|q|,|r|<=51, checked on |i|<=300):")
    for h in hits:
        print("   ", h)
    assert hits == [("{2i-1,2i}", 3, 1), ("{2i,2i+1}", 3, -1)] or sorted(hits) == sorted([("{2i-1,2i}", 3, 1), ("{2i,2i+1}", 3, -1)])
    print("  symbolic: f(2i-1)+f(2i) = (q+1)i + (r-q)/2 ; f(2i)+f(2i+1) = (q+1)i + (q+r)/2;")
    print("  equality with 4i-1 resp. 4i+1 for all i forces q=3 and r=1 resp. r=-1.  PROVED.")
    # (b) arbitrary matchings.  Pair types: even-even {2i,2j} needs i+j=0; odd-odd {n,m} needs
    # (q-2)(n+m) = -2r; even-odd {2i,n} needs 2i = (q-2)n + r.
    print("  arbitrary matchings: EE pairs need i+j=0, OO pairs need (q-2)(n+m) = -2r, EO pairs are {n,(q-2)n+r}.")
    print("  On Z the reflection pairs make degenerate matchings possible, e.g. 5x+3: {0,-1}, {n,-2-n}, {2i,-2i}:")
    f53 = lambda n: f(5, 3, n)
    assert f53(0) + f53(-1) == -1 and all(f53(n) + f53(-2 - n) == -2 for n in range(-99, 100, 2)) \
        and all(f53(2 * i) + f53(-2 * i) == 0 for i in range(1, 100))
    print("    (checked).  On the positive integers there are no EE pairs and only finitely many OO pairs, so all but")
    print("    finitely many evens must be (q-2)n + r with n odd: this forces |q-2| = 1, and q = 1 gives only finitely")
    print("    many positive partners.  Census of unmatched integers in [X/2, X], X = 20000, for the EO matching:")
    X = 20000
    for q in (1, 3, 5, 7, 9):
        for r in (-3, -1, 1, 3, 5):
            evens_hit = set((q - 2) * n + r for n in range(1, 2 * X, 2))
            un_e = sum(1 for e in range(X // 2 + (X // 2) % 2, X + 1, 2) if e not in evens_hit)
            if q == 3 and r in (1, -1):
                assert un_e == 0
            if q in (3,) or r == 1:
                print(f"    q={q}, r={r:+d}: evens in [X/2,X] with no odd partner: {un_e} of {X // 4 + 1}")
    print("  => up to finitely many exceptions, a sum-preserving matching of the positive integers exists iff q = 3")
    print("     (partner n + r); it is exact on N iff (q,r) = (3,1) with pairs {2i-1,2i}, and on N_0 for (3,-1) with {2i,2i+1}.")
    for r in (1, -1, 3, 5, -3):
        # positive integers: odds n>=1 matched to n+r must be >=1 and cover all evens >= 2
        matched_evens = set(n + r for n in range(1, 2001, 2) if n + r >= 1)
        uncovered = [e for e in range(2, 40, 2) if e not in matched_evens]
        bad_odds = [n for n in range(1, 40, 2) if n + r < 1]
        print(f"    3x{r:+d} on N: evens uncovered <40: {uncovered[:6]}  odds without partner <40: {bad_odds}")
        if r == 1:
            assert not uncovered and not bad_odds
            assert all(T(n) + T(n + 1) == 2 * n + 1 for n in range(1, 200001, 2))
    print("  CHECK T(2i-1)+T(2i) = 4i-1 for i <= 100000: PASS")

# ---------------------------------------------------------------- F2
def F2():
    hdr("F2  arithmetic-mean vs geometric-mean step factors; Lundberg exponents")
    print("  qx+1 shortcut map, parity fair coin: AM factor (q/2+1/2)/2 = (q+1)/4 ; GM factor sqrt(q)/2")
    for q in (1, 3, 5, 7, 9):
        am = (q + 1) / 4
        gm = math.sqrt(q) / 2
        # Lundberg exponent: ((q/2)^t + (1/2)^t)/2 = 1, t>0  (exists iff log-drift < 0 and q/2 > 1)
        th = None
        if q > 2 and gm < 1:
            lo, hi = 1e-9, 60.0
            g = lambda t: ((q / 2) ** t + 0.5 ** t) / 2 - 1
            for _ in range(200):
                mid = (lo + hi) / 2
                if g(mid) < 0: lo = mid
                else: hi = mid
            th = (lo + hi) / 2
        print(f"    q={q}: AM={am:.4f}  GM={gm:.6f}  log-drift={math.log(gm):+.6f}  Lundberg theta={th}")
    print("  q=3 is the only odd q with AM factor 1 (a martingale), and theta(3)=1 exactly because")
    print("  ((3/2)^1+(1/2)^1)/2 = 1.  Drift log(sqrt3/2) = -0.143841 = the AM-GM gap.")
    # Collatz 1932 permutation g: residues mod 3 uniform: factors 2/3, 4/3, 4/3 ; inverse U: 3/2,3/2,3/4,3/4 (mod 4)
    amg = (2/3 + 4/3 + 4/3) / 3
    gmg = (2/3 * 4/3 * 4/3) ** (1/3)
    amu = (3/2 + 3/2 + 3/4 + 3/4) / 4
    gmu = (3/2 * 3/2 * 3/4 * 3/4) ** (1/4)
    print(f"  Collatz permutation g: AM={amg:.6f} GM={gmg:.6f} (log {math.log(gmg):+.6f}) ; "
          f"inverse U: AM={amu:.6f} GM={gmu:.6f} (log {math.log(gmu):+.6f})")
    print("  both directions expand in AM and GM: g is NOT fair; infinite orbits are the expected behaviour.")
    print("  Duality (elementary; the image-density identity is Kohl, JSC 2008 Remark 1.2 [P]): for an rcwa map with")
    print("  branches n -> (a n + b)/c on the classes mod m, the Haar mean of the inverse ratio c/|a| is its image density,")
    print("  which equals 1 for every permutation; by Jensen its log-drift E log(|a|/c) >= 0, with equality iff |a|/c is")
    print("  constant.  T is the mirror case: the Haar mean of the FORWARD ratio a/c is 1 (the pair-sum property), so its")
    print("  log-drift is < 0.  Check:")
    from fractions import Fraction as Fr_
    for name, br in (("T (3x+1)", [(1, 2), (3, 2)]), ("g (Collatz 1932)", [(2, 3), (4, 3), (4, 3)]),
                     ("U = g^-1", [(3, 2), (3, 2), (3, 4), (3, 4)]), ("5x+1", [(1, 2), (5, 2)])):
        fwd = sum(Fr_(a, c) for a, c in br) / len(br)
        inv = sum(Fr_(c, a) for a, c in br) / len(br)
        drift = sum(math.log(a / c) for a, c in br) / len(br)
        print(f"    {name:18s}: Haar mean of a/c = {str(fwd):>5s}, of c/a (image density) = {str(inv):>5s}, log-drift = {drift:+.6f}")

# ---------------------------------------------------------------- F3
LOG2, LOG3 = math.log(2), math.log(3)

def haar_max_exact(k, lam):
    """Exact Haar{max_{0<=j<=k} 3^{a_j}/2^j >= lam} for rational lam > 1 (M_0 = 1 < lam)."""
    lam = Fr(lam)
    u, v = lam.numerator, lam.denominator
    alive = {0: Fr(1)}  # a -> mass, at step j, among paths with max < lam so far
    killed = Fr(0)
    for j in range(1, k + 1):
        new = {}
        for a, m in alive.items():
            for a2 in (a, a + 1):
                mm = m / 2
                # M_j = 3^{a2}/2^j >= u/v  <=>  3^{a2} v >= u 2^j
                if 3 ** a2 * v >= u * 2 ** j:
                    killed += mm
                else:
                    new[a2] = new.get(a2, 0) + mm
        alive = new
    return killed

def haar_max_float(k, lam, cut=1e-14):
    """Float DP with truncation of states whose M_j < cut*lam; by Doob the truncated mass m
    can later reach lam with probability <= M_j/lam, so the error is <= sum m*M_j/lam."""
    L = math.log(lam)
    alive = {0: 1.0}
    killed = 0.0
    err = 0.0
    for j in range(1, k + 1):
        new = {}
        for a, m in alive.items():
            for a2 in (a, a + 1):
                mm = m * 0.5
                x = a2 * LOG3 - j * LOG2
                if x >= L - 1e-15:
                    killed += mm
                elif x < L + math.log(cut):
                    err += mm * math.exp(x - L)
                else:
                    new[a2] = new.get(a2, 0.0) + mm
        alive = new
    # remaining alive mass at step k: its future contribution <= sum m*M_k/lam
    tail = sum(m * math.exp(a * LOG3 - k * LOG2 - L) for a, m in alive.items())
    return killed, err + tail

def F3():
    hdr("F3  the Haar martingale M_j = 3^{a_j}/2^j: exact maximal probabilities vs Doob")
    print("  Under Haar measure on Z_2 the T-parity bits are iid fair (Terras/Everett bijection;")
    print("  Lagarias 1985 Thm K), so E[M_{j+1}|F_j] = M_j((3/2)+(1/2))/2 = M_j: a nonnegative martingale.")
    print("  Doob's maximal inequality: P_k(lam) := Haar{max_{j<=k} M_j >= lam} <= 1/lam.")
    print()
    print("  exact rationals, k = 60:")
    for lam in (Fr(3, 2), Fr(2), Fr(3), Fr(4), Fr(8), Fr(16)):
        P = haar_max_exact(60, lam)
        print(f"    lam={str(lam):>4}: P_60 = {float(P):.12f}   lam*P = {float(lam*P):.6f}   Doob bound 1/lam = {float(1/lam):.6f}")
        assert P <= 1 / lam
    print()
    print("  floats, k = 4000, with Doob-controlled truncation error (limit k -> infinity):")
    rows = []
    for lam in (2.0, 4.0, 10.0, 100.0, 1000.0, 1e4, 1e5, 1e6):
        P, e = haar_max_float(4000, lam)
        rows.append((lam, P, e))
        print(f"    lam={lam:>9.0f}: P_inf ~ {P:.10f} (+ at most {e:.2e})   lam*P = {lam*P:.6f}")
        assert P <= 1 / lam + 1e-12
    print("  lam*P_inf(lam) stays in a bounded band below 1: Pareto(1) tail of sup_j M_j (Cramer-Lundberg")
    print("  with exponent theta = 1, which is exactly the martingale property); the band oscillation")
    print("  reflects the two-step lattice log(3/2), log 2 at moderate lam.")

# ---------------------------------------------------------------- F4
def F4():
    hdr("F4  Doob's convergence theorem => Terras: Haar{tau > k} -> 0, tau = first j with M_j < 1")
    print("  M_j -> M_inf a.s. (nonnegative martingale); |log M_{j+1} - log M_j| >= log(3/2) forces")
    print("  M_inf = 0 a.s., hence tau < infinity a.s., hence (Terras 1976 / Everett 1977) the natural density")
    print("  of {n : sigma(n) < infinity} is lim_k Haar{tau <= k} = 1.  Exact values of Haar{tau > k}:")
    alive = {0: Fr(1)}
    for j in range(1, 201):
        new = {}
        for a, m in alive.items():
            for a2 in (a, a + 1):
                if 3 ** a2 >= 2 ** j:   # M_j >= 1: still not stopped
                    new[a2] = new.get(a2, 0) + m / 2
        alive = new
        if j in (1, 2, 5, 10, 20, 50, 100, 200):
            s = sum(alive.values())
            print(f"    k={j:4d}: Haar{{tau > k}} = {float(s):.6e}")
    print("  (decay is exponential, rate given by the large-deviation cost of staying above slope log_3 2).")

# ---------------------------------------------------------------- F5
def F5():
    hdr("F5  deterministic Doob: density of {n : max_{j<=k} T^j(n) >= lam n} equals P_k(lam) <= 1/lam")
    k = 20
    lams = (Fr(2), Fr(4), Fr(8))
    Pk = {lam: haar_max_exact(k, lam) for lam in lams}
    base = 2 ** 40
    cnt = {lam: 0 for lam in lams}
    cnt_small = {lam: 0 for lam in lams}
    for n0 in range(2 ** k):
        n = base + n0
        x = n; mx = n
        for j in range(k):
            x = T(x)
            if x > mx: mx = x
        for lam in lams:
            if mx * lam.denominator >= lam.numerator * n:
                cnt[lam] += 1
        n = n0 + 1
        x = n; mx = n
        for j in range(k):
            x = T(x)
            if x > mx: mx = x
        for lam in lams:
            if mx * lam.denominator >= lam.numerator * n:
                cnt_small[lam] += 1
    for lam in lams:
        exact = Pk[lam] * 2 ** k
        print(f"  lam={lam}: #n in [2^40, 2^40+2^20) = {cnt[lam]}  vs  2^20 * P_20 = {exact}  "
              f"(equal: {Fr(cnt[lam]) == exact})   ;  #n in [1,2^20] = {cnt_small[lam]}  ; Doob 2^20/lam = {2**k/lam}")
        assert Fr(cnt[lam]) == exact
    print("  PROVED (elementary): for every k and lam > 1, d{n : max_{j<=k} T^j(n) >= lam n} = P_k(lam) <= 1/lam,")
    print("  since T^j(n)/n = M_j(n) + R_j/(2^j n) with 0 < R_j/2^j <= 3^{a_j}/4, and the k-prefix of n is Haar-uniform.")

# ---------------------------------------------------------------- F6
def F6():
    hdr("F6  exactness only in the leading coefficient: sum_{n<=X} T^k(n) - sum_{n<=X} n")
    for k in range(1, 7):
        X = 2 ** (k + 12)
        s = 0
        for n in range(1, X + 1):
            x = n
            for _ in range(k):
                x = T(x)
            s += x
        base = X * (X + 1) // 2
        print(f"  k={k}: X=2^{k+12}: sum T^k - sum n = {s - base}  ( = {Fr(s - base, X)} per element )")
    print("  k=1: exactly 0 (the pair-sum identity on {2i-1,2i}); k>=2: a linear excess ~ c_k X from the +1's.")

if __name__ == "__main__":
    F1(); F2(); F3(); F4(); F5(); F6()
    print("\nDONE martingale")
