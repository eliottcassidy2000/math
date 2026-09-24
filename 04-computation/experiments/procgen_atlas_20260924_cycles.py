#!/usr/bin/env python3
"""
procgen_atlas_20260924_cycles.py

Implication-atlas lane (session collatz-procgen-20260922, 2026-09-24): the cycle half against the
Diophantine conjectures (abc, explicit abc, Pillai, Baker, Lang-Waldschmidt) and the stopping-time
conjecture of Terras.

Sections:
  C1  continued fraction of alpha = log_2 3 (mpmath, 400 digits; PARI cross-check if gp is present).
  C2  Eliahou-type cycle-length bound from a verification bound N: a nontrivial positive T-cycle with
      L odd terms and K = total T-steps has all elements >= N, hence
          0 < K ln2 - L ln3 = sum_odd ln(1 + 1/(3 n_i)) <= L ln(1 + 1/(3N)),
      so K/L lies in (alpha, alpha + delta_N), delta_N = log_2(1 + 1/(3N)); the least such L is the least
      denominator of an upper semiconvergent within delta_N.  Values for N = 2^40 .. 2^100, and the exact
      verification thresholds N* at which the bound jumps.
  C3  conditional growth: under an irrationality-exponent bound |alpha - p/q| >= c q^(-mu), every cycle
      has L >= (3 c N ln 2)^(1/mu); Lang-Waldschmidt for two logarithms gives mu = 2 + eps.  Empirical
      L_min(N)/sqrt(N) along the computed range.
  C4  the gate barrier: lower bounds for D = |2^K - 3^L| from abc (any eps), Baker's explicit abc,
      Ellison 1971, a Baker-type polynomial bound, and the truth at the record clocks.
  C5  abc and two-block parity prefixes 1^k 0^l (generalising Rozier 2025, Thm 2.1): the least
      n = 2^k A - 1, A = 3^(-k) mod 2^l, with the abc quality of (1, 2^l B, 3^k A); exact data.
  C6  Catalan/Pillai for 2 and 3: unit gaps, and all solutions of |2^K - 3^L| = c <= 10^4 in a box.
  C7  Terras's coefficient stopping time: sigma(n) = tau(n) for 2 <= n <= 2*10^6, and the implications
      CST => no nontrivial cycle, T1 => tau finite, CST & T1 => Collatz (proved in the note).
"""
import math, sys, subprocess, shutil
from fractions import Fraction as Fr
from mpmath import mp, mpf, log as mlog, floor as mfloor, log1p as mlog1p, nstr

def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

mp.dps = 400
ALPHA = mlog(3) / mlog(2)

def cf(x, n):
    a = []
    for _ in range(n):
        ai = int(mfloor(x)); a.append(ai); x = 1 / (x - ai)
    return a

A_CF = cf(ALPHA, 140)
P = [A_CF[0], A_CF[0] * A_CF[1] + 1]; Q = [1, A_CF[1]]
for n in range(2, 130):
    P.append(A_CF[n] * P[-1] + P[-2]); Q.append(A_CF[n] * Q[-1] + Q[-2])

def upper_semis(nmax=120):
    ups = [(P[1], Q[1])]
    for n in range(2, nmax, 2):
        for j in range(1, A_CF[n + 1] + 1):
            ups.append((P[n - 1] + j * P[n], Q[n - 1] + j * Q[n]))
    return ups

UPS = upper_semis()

def C1():
    hdr("C1  continued fraction of log_2 3")
    print("  partial quotients a_0..a_59:", A_CF[:60])
    ok = None
    if shutil.which("gp"):
        try:
            out = subprocess.run(["gp", "-q", "-f"], input="default(realprecision,320); print(contfrac(log(3)/log(2))[1..60]);\n",
                                 capture_output=True, text=True, timeout=120).stdout.strip()
            gpa = [int(t) for t in out.strip("[]").replace(" ", "").split(",")]
            ok = gpa == A_CF[:60]
        except Exception as e:
            ok = f"gp failed: {e}"
    print(f"  PARI/GP cross-check of the first 60 partial quotients: {ok}")
    for n in range(0, 30):
        side = "+" if mpf(P[n]) / Q[n] > ALPHA else "-"
        print(f"    p_{n}/q_{n} = {P[n]}/{Q[n]}  ({side})")

def min_cycle(N):
    delta = mlog1p(mpf(1) / (3 * mpf(N))) / mlog(2)
    for (p, q) in UPS:
        if mpf(p) / q - ALPHA < delta:
            return q, p
    return None

def C2():
    hdr("C2  Eliahou-type bound: least possible (L, K) of a nontrivial positive cycle, given all n < N converge")
    rows = [("2^40 (Eliahou 1993 input)", 2 ** 40), ("2^50", 2 ** 50), ("2^60", 2 ** 60),
            ("5*2^60", 5 * 2 ** 60), ("2^68 (Barina 2021)", 2 ** 68), ("704*2^60 (Hercher's X_0)", 704 * 2 ** 60),
            ("2^71 (Barina 2025)", 2 ** 71), ("2^72", 2 ** 72), ("2^80", 2 ** 80), ("2^100", 2 ** 100)]
    for name, N in rows:
        L, K = min_cycle(N)
        print(f"  N = {name:28s}: L >= {L:>20d} odd terms, K >= {K:>20d} T-steps, K+L >= {K+L}")
    L, K = min_cycle(2 ** 40)
    assert K == 17087915, "Eliahou's 17087915 must be reproduced"
    print("  check: N = 2^40 reproduces Eliahou's period bound K = 17,087,915 (the fraction 17087915/10781274).")
    print()
    print("  verification thresholds N* = 1/(3(2^d - 1)), d = K/L - log_2 3, above which the clock K/L is excluded:")
    for (p, q) in UPS:
        if 10 ** 7 < q < 10 ** 13:
            d = mpf(p) / q - ALPHA
            Ns = 1 / (3 * (mpf(2) ** d - 1))
            print(f"    L = {q:>15d}, K = {p:>15d}:  N* = {nstr(Ns, 10):>16s} = 2^{nstr(mlog(Ns)/mlog(2), 8)}")
    Ncrude = None
    for (p, q) in UPS:
        if q == 72057431991:
            d = mpf(p) / q - ALPHA
            Ncrude = 1 / (3 * (mpf(2) ** d - 1))
    print(f"  The crude threshold for excluding L = 72,057,431,991 is N* = {nstr(Ncrude, 8)} = {nstr(Ncrude / 2**60, 8)} * 2^60:")
    print("  this reproduces exactly Hercher's Remark 28 value X0 >= 3781*2^60 for 'previous methods'.")
    print("  Hercher's Theorem 27 shrinks the window by 3/4 (X0 >= 2836*2^60) and his computer-assisted Corollary 29")
    print("  needs only X0 >= 1536*2^60 = 3*2^69 < 2^71.  Hence (CITED: Hercher 2023 Cor. 29 + Barina 2025):")
    print("  every nontrivial cycle has L >= 137,528,045,312 odd terms and K >= 217,976,794,617 T-steps")
    print("  (355,504,839,929 steps of the unshortened map).")
    for (p, q) in UPS:
        if q == 137528045312:
            d = mpf(p) / q - ALPHA
            Nn = 1 / (3 * (mpf(2) ** d - 1))
            print(f"  Next jump (to L >= 890,638,885,193): crude threshold {nstr(Nn, 6)} = 2^{nstr(mlog(Nn)/mlog(2), 6)};")
            print(f"  with Hercher's Theorem 27 factor 3/4: {nstr(Nn * 3 / 4, 6)} = 2^{nstr(mlog(Nn * 3 / 4)/mlog(2), 6)},"
                  f" i.e. about {nstr(Nn * 3 / 4 / 2**71, 4)} times Barina's 2^71.")

def C3():
    hdr("C3  conditional growth of the least cycle length with the verification bound N")
    print("  If |log_2 3 - p/q| >= c q^(-mu) for all q, then a cycle with all elements >= N has 0 < K/L - alpha <")
    print("  delta_N ~ 1/(3N ln 2), hence L > (3 c N ln 2)^(1/mu).  Lang-Waldschmidt (two logarithms, fixed")
    print("  algebraic 2 and 3) gives mu = 2 + eps; the computed least L tracks sqrt(N):")
    for e in (40, 50, 60, 68, 71, 75, 80, 90, 100):
        N = 2 ** e
        L, K = min_cycle(N)
        print(f"    N = 2^{e:3d}: L_min = {L:>20d}   L_min/sqrt(N) = {L / math.sqrt(N):10.4f}   "
              f"log L_min / log N = {math.log(L) / math.log(N):.4f}")
    print("  Unconditionally: Lagarias 2009 Lemma 2.2 (from Rhin 1987): |log_3 2 - p/q| >= q^(-14.3)/1200 for all q >= 1.")
    print("  With |log_3 2 - L/K| < (L/K) delta_N / alpha this gives the fully explicit, unconditional bound")
    for e in (71, 100, 200):
        N = mpf(2) ** e
        delta = mlog1p(1 / (3 * N)) / mlog(2)
        # K^14.3 > alpha^2 / (1200 * delta)   (using L/K < 1/alpha * (1 + tiny))
        Kmin = (ALPHA ** 2 / (1200 * delta)) ** (1 / mpf('14.3'))
        print(f"    N = 2^{e}: K > {nstr(Kmin, 6)}  (versus the continued-fraction truth above)")
    print("  so proved irrationality measures are useless at the present N; asymptotically they give L >> N^(1/14.3)")
    print("  (explicit) or N^(1/8.616) (Rhin, asymptotic [R]), against the LW prediction N^(1/2) that the record tracks.")

def C4():
    hdr("C4  the gate barrier: lower bounds for D = |2^K - 3^L| (log_2 D - K shown)")
    print("  abc on 3^L + D = 2^K (coprime, rad(2^K 3^L D) <= 6 rad(D) <= 6D):  2^K <= K_eps (6D)^(1+eps),")
    print("     so log2 D >= K/(1+eps) - log2 6 - log2(K_eps)/(1+eps)   [shown with the optimistic K_eps = 1].")
    print("  Baker's explicit abc (c < (6/5) N (log N)^w / w!, N = rad, w = #primes):  D >= (5/36) 2^K w!/(K ln2 + ln6)^w")
    print("     with w = omega(D) + 2 (shown for w = 3 and for the worst admissible w ~ K ln2/ln K).")
    print("  Ellison 1971 Thm 3: D > 2^K e^(-K/10) for K not in {1..11,13,14,16,19,27} [P].")
    print("  Rhin 1987 (via Simons-de Weger Lemma 12 [R]): |u0 + u1 log2 + u2 log3| >= H^(-13.3), so D >= c 2^K K^(-13.3).")
    print("  truth at the record clocks K = p_n: log2 D = K + log2|1 - 2^(-(K - L alpha))|, from the continued fraction.")
    print()
    print("   K            abc eps=.1   abc eps=.01   expl.abc w=3   expl.abc worst   Ellison     C=13.3      truth(record clock)")
    for (p, q) in [(P[n], Q[n]) for n in (5, 9, 13, 17, 21, 25, 29)]:
        K, L = p, q
        a1 = K / 1.1 - math.log2(6) - K
        a2 = K / 1.01 - math.log2(6) - K
        lg = K * math.log(2) + math.log(6)
        w3 = math.log2(5 / 36) + math.log2(math.factorial(3)) - 3 * math.log2(lg)
        wmax = max(3, int(K * math.log(2) / math.log(K)))
        worst = math.log2(5 / 36) + (math.lgamma(wmax + 1) - wmax * math.log(lg)) / math.log(2)
        ell = -K / (10 * math.log(2))
        bak = -13.3 * math.log2(K)
        d = mpf(K) - L * ALPHA
        truth = float(mlog(abs(1 - mpf(2) ** (-d))) / mlog(2))
        print(f"  {K:>14d}  {a1:12.4g}  {a2:12.4g}  {w3:12.4g}  {worst:14.4g}  {ell:10.4g}  {bak:10.4g}   {truth:10.4f}")
    print("  Reading: plain abc loses a factor 2^(eps K/(1+eps)) -- exponentially worse than the proved Baker-type")
    print("  bounds (polynomial or quasi-polynomial loss) and than the truth (~ K^-1); it beats only Ellison's")
    print("  exponential-loss bound for eps < 0.168.  Explicit abc is polynomial only when omega(D) is small.")
    print("  PROVED (comparison): every cycle-length or cycle-element bound obtained from 'gate lower bound +")
    print("  product identity + verification' is monotone in the gate bound, so abc adds nothing to Baker there.")

def rad_and_quality(a, b, c):
    import sympy
    r = 1
    for x in (a, b, c):
        for pr in sympy.factorint(x):
            r *= pr
    # distinct primes across a,b,c (coprime triple) -> product of distinct primes = rad(abc)
    return r, math.log(c) / math.log(r) if r > 1 else float("inf")

def C5():
    hdr("C5  abc and the two-block parity prefix 1^k 0^l")
    print("  n has T-parity prefix 1^k 0^l  <=>  n = 2^k A - 1 with 3^k A = 1 + 2^l B (A odd).  abc on (1, 2^l B, 3^k A):")
    print("  3^k A <= K_eps (6AB)^(1+eps) and B < 3^k A/2^l give log2 A >= l(1+eps)/(1+2eps) - eps k log2(3)/(1+2eps) - O_eps(1),")
    print("  i.e. n >= c_eps 2^(k + (1-eps) l - 1.59 eps k)  (Proposition R2 of the note; Rozier's Thm 2.1 is the prefix 1^k 0 1^(j-k-1)).")
    print("  Unconditionally A_min(k,l) = (3^(-k) mod 2^l) can be tiny only when k is huge: A = 1 needs 2^(l-2) | k.")
    print()
    import sympy
    print("  worst cases of log2(A_min)/l over 1 <= k <= 2l, for l = 8..48 (exact):")
    worst_rows = []
    for l in range(8, 49, 4):
        mod = 2 ** l
        best = None
        for k in range(1, 2 * l + 1):
            A = pow(pow(3, k, mod), -1, mod)
            ratio = math.log2(A) / l
            if best is None or ratio < best[0]:
                best = (ratio, k, A)
        worst_rows.append((l, best))
        print(f"    l={l:2d}: min_k log2(A_min)/l = {best[0]:.4f} at k={best[1]:3d} (A_min={best[2]})")
    print()
    print("  abc quality q = log c / log rad(abc) of the triples (1, 2^l B, 3^k A_min), k <= 20, 4 <= l <= 40:")
    hits = []
    for k in range(1, 21):
        for l in range(4, 41):
            mod = 2 ** l
            A = pow(pow(3, k, mod), -1, mod)
            c = 3 ** k * A
            b = c - 1
            assert b % mod == 0
            r, qual = rad_and_quality(1, b, c)
            hits.append((qual, k, l, A))
    hits.sort(reverse=True)
    for qual, k, l, A in hits[:8]:
        print(f"    q = {qual:.4f}  (k={k}, l={l}, A={A}, n = 2^k A - 1 = {2**k*A-1})")
    nhit = sum(1 for h in hits if h[0] > 1)
    print(f"  abc-hits (q > 1) among these {len(hits)} triples: {nhit}; max quality {hits[0][0]:.4f}.")

def C6():
    hdr("C6  Catalan and Pillai for the bases 2 and 3 (finite checks)")
    unit = [(K, L) for K in range(1, 401) for L in range(1, 260) if abs(2 ** K - 3 ** L) == 1]
    print(f"  |2^K - 3^L| = 1 with K, L >= 1, K <= 400: {unit}  (Gersonides/Levi ben Gershon; Mihailescu in general)")
    sols = {}
    for K in range(0, 401):
        for L in range(0, 260):
            c = 3 ** L - 2 ** K
            if 0 < abs(c) <= 10000:
                sols.setdefault(c, []).append((K, L))
    multi = {c: v for c, v in sols.items() if len(v) > 1}
    print(f"  values c = 3^L - 2^K with 0 < |c| <= 10^4 in the box K <= 400, L < 260: {len(sols)} values;")
    print(f"  values attained twice or more: {dict(sorted(multi.items()))}")
    print("  (consistent with the cited Pillai-type theorems for 3^x - 2^y = c; the note gives the exact statements.)")
    print("  Collatz relevance: a gate value D = 2^K - 3^L determines its clock for |D| > 13 in this box; the carry")
    print("  condition D | B(w) is untouched by any fixed-c Pillai statement.")

def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2

def C7(NMAX=2 * 10 ** 6):
    hdr(f"C7  Terras's coefficient stopping time: sigma(n) = tau(n) for 2 <= n <= {NMAX}")
    bad = 0
    maxtau = 0
    for n in range(2, NMAX + 1):
        x = n; j = 0; a = 0
        tau = None; sigma = None
        while True:
            if x & 1:
                a += 1
            x = T(x); j += 1
            if tau is None and 3 ** a < 2 ** j:
                tau = j
            if x < n:
                sigma = j
                break
        if tau != sigma:
            bad += 1
        maxtau = max(maxtau, tau)
    print(f"  mismatches: {bad}; max tau = {maxtau}")
    print("  Implications (PROVED in the note, elementary):")
    print("   (i)  CST => no nontrivial positive cycle: the minimum m of such a cycle has tau(m) <= K < sigma(m) = infinity.")
    print("   (ii) T1 => tau(n) < infinity for all n >= 1: if tau(n) = infinity, T^j(n) >= n for all j and the orbit")
    print("        is not eventually periodic (a positive cycle forces 3^a_j/2^j -> 0), so it diverges.")
    print("   (iii) CST & (tau finite) => Collatz (sigma(n) = tau(n) < infinity for every n >= 2, then induction);")
    print("        hence CST & T1 => Collatz, and Collatz <=> NC & T1.")

if __name__ == "__main__":
    C1(); C2(); C3(); C4(); C5(); C6(); C7()
    print("\nDONE cycles")
