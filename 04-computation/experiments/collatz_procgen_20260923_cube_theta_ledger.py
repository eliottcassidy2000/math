#!/usr/bin/env python3
"""
collatz_procgen_20260923_cube_theta_ledger.py

Cube-theta lane (HYP-9127), exact ledgers.  Sections:

  L0  constants and thresholds
  L1  independent check of the Y3 block identity  Phi_2(Y3) = -1 - 512/(3^9-2^10) - (256/3^9)(X-1)
  L2  Liouville ledger for X = sum rho^(k^3): exact valuations and heights of partial sums
  L3  Theorem Q: exact Tschakaloff-Zudilin linear forms for quadratic position families
      (squares, triangular, both one-sided pentagonal branches, 2k^2+k) under several maps
  L4  the cubic analogue of the Zudilin construction (kernel of M_ts = rho^(3ts(t+s))):
      exact exponents -- far below 1
  L5  threshold map: where mu_bar < phi bites; nearest-to-critical cubic instances
  L6  the near-critical cube-swap number: how far the Liouville certificates reach

Usage:  python3 collatz_procgen_20260923_cube_theta_ledger.py [--quick]
All certificates are exact integer computations.
"""
import math
import sys
import time
from fractions import Fraction

import gmpy2
from gmpy2 import mpz

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from collatz_procgen_20260923_cube_theta_core import (  # noqa: E402
    PHI, LOG2_3, L_Y3, M0_Y3, mubar, v2, X_mod, theta_mod, partial_sum_fraction,
    log2abs, cube_word, is_cube, bernstein_3x1_mod, parity_check_3x1)

QUICK = "--quick" in sys.argv


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


# ----------------------------------------------------------------------------- L0
def section_L0():
    hdr("L0  constants")
    mu = mubar(L_Y3, M0_Y3)
    print(f"rho = 2^10/3^9 ; |rho|_2 = 2^-10 ; H(rho) = 3^9 = {M0_Y3}")
    print(f"mu_bar(Y3, 3x+1) = log2(3^9)/10 = {mu:.6f}   (= 9/10 * log2 3)")
    print(f"phi = {PHI:.6f} ; Lemma P / Theorem Q need mu_bar < phi ; unshifted (m=0) needs mu_bar < 3/2")
    print(f"Liouville (partial sums) for sum rho^(k^d), d>=2: works iff mu_bar = 1 (exponent ratio -> 1/mu_bar = {1/mu:.6f})")
    print(f"formal Pade budget per unit degree at W = 2E (diagonal):  2L - log2 M0 = {2*L_Y3 - math.log2(M0_Y3):.4f} bits")
    print(f"formal budget at W = w E: L*w - log2 M0 ; zero at w = mu_bar = {mu:.4f}")


# ----------------------------------------------------------------------------- L1
def section_L1():
    hdr("L1  independent check of the Y3 identity (Bernstein formula + parity iteration)")
    N = 6000 if QUICK else 20000
    nblocks = N // 10 + 2
    word = cube_word(nblocks, swap=is_cube)
    t0 = time.time()
    phi_word = bernstein_3x1_mod(word, N)
    mod = mpz(1) << N
    X = X_mod(N)
    inv = gmpy2.invert(mpz(3 ** 9 - 2 ** 10), mod)
    inv39 = gmpy2.invert(mpz(3 ** 9), mod)
    rhs = (-1 - 512 * inv - 256 * inv39 * (X - 1)) % mod
    print(f"N = {N} bits, word length {len(word)} letters, nonzero cubes used: "
          f"{sum(1 for n in range(nblocks) if is_cube(n))}")
    print(f"Phi(word) == RHS mod 2^N : {phi_word == rhs}   ({time.time()-t0:.2f}s)")
    agree = parity_check_3x1(phi_word, word, N)
    print(f"parity iteration of Phi(word) mod 2^N reproduces the first {agree} letters of Y3 (expected {N})")
    print(f"X mod 2^64 = {int(X % (1 << 64)):#018x}")
    import mpmath
    mpmath.mp.dps = 50
    rr = mpmath.mpf(2) ** 10 / mpmath.mpf(3) ** 9
    XR = mpmath.nsum(lambda k: rr ** (k ** 3), [0, mpmath.inf])
    PhiR = -1 - mpmath.mpf(512) / (3 ** 9 - 2 ** 10) - mpmath.mpf(256) / 3 ** 9 * (XR - 1)
    print(f"real analogue: X_R = sum (2^10/3^9)^(k^3) = {mpmath.nstr(XR, 40)}")
    print(f"               Phi_R(Y3) (same series summed in R) = {mpmath.nstr(PhiR, 40)}")
    # the 'rational' test in the same data: a rational a/b with odd b<2^32 and |a|<2^32 would give b*X - a == 0 mod 2^N
    return X


# ----------------------------------------------------------------------------- L2
def section_L2():
    hdr("L2  Liouville ledger for X: partial sums S_K = A_K/3^(9K^3)")
    Kmax = 9 if QUICK else 12
    Nbig = L_Y3 * (Kmax + 2) ** 3 + 64
    X = X_mod(Nbig)
    mod = mpz(1) << Nbig
    mu = mubar(L_Y3, M0_Y3)
    print(" K |  v2(D X - A)  | 10(K+1)^3 | log2 max(|A|,D) | gain = v2 - log2 H | v2/log2 H | 1/mu_bar")
    first_neg = None
    for K in range(1, Kmax + 1):
        A, D = partial_sum_fraction(K)
        val = v2((D * X - A) % mod)
        lh = max(log2abs(A), log2abs(D))
        print(f"{K:2d} | {val:13d} | {10*(K+1)**3:9d} | {lh:15.1f} | {val - lh:18.1f} | {val/lh:9.4f} | {1/mu:.4f}")
        if val - lh < 0 and first_neg is None:
            first_neg = K
    print(f"gain < 0 from K = {first_neg} on; asymptotically gain = L[(K+1)^3 - mu_bar K^3] + O(log K)"
          " ~ -(mu_bar-1) L K^3: partial sums certify only finitely many heights.")


# ----------------------------------------------------------------------------- L3
def qbinom_products(n):
    """Coefficients C_k(q) of R_n(T) = prod_{j=1}^n (1 - q^j T) as dicts {qdeg: coeff}."""
    C = [{0: 1}]  # polynomial in T with coefficients in Z[q]
    for j in range(1, n + 1):
        new = [dict() for _ in range(len(C) + 1)]
        for k, ck in enumerate(C):
            for d, c in ck.items():
                new[k][d] = new[k].get(d, 0) + c
                new[k + 1][d + j] = new[k + 1].get(d + j, 0) - c
        C = [{d: c for d, c in ck.items() if c != 0} for ck in new]
    return C


def tschakaloff_form(n, m, a2, s, L, M0, Nbig, theta):
    """
    Quadratic family P(l) = a2*l*(l-1)/2 + s*l  (alpha = a2/2, alpha+beta = s), theta = sum_l rho^P(l) mod 2^Nbig.
    q = rho^-a2, z = rho^s.  Zudilin form I_n = A_n theta - B_n with
      A_n = sum_k C_k(q) z^-k q^(k(k-1)/2 + k m),  B_n = sum_k (same coeff) * sum_{l <= k+m} z^l q^(-l(l-1)/2).
    Returns (v2 of the cleared form, log2 height, predicted v2, lo, hi).
    """
    P = lambda l: a2 * l * (l - 1) // 2 + s * l
    C = qbinom_products(n)
    Acoef = {}  # exponent of rho -> integer coefficient
    Bcoef = {}
    for k, ck in enumerate(C):
        for d, c in ck.items():
            e0 = -a2 * (d + k * (k - 1) // 2 + k * m) - s * k
            Acoef[e0] = Acoef.get(e0, 0) + c
            for l in range(0, k + m + 1):
                e = e0 + P(l)
                Bcoef[e] = Bcoef.get(e, 0) + c
    Acoef = {e: c for e, c in Acoef.items() if c}
    Bcoef = {e: c for e, c in Bcoef.items() if c}
    exps = list(Acoef) + list(Bcoef)
    lo = min(0, min(exps))
    hi = max(0, max(exps))
    def clear(coef):
        tot = mpz(0)
        for e, c in coef.items():
            tot += c * (mpz(2) ** (L * (e - lo))) * (mpz(M0) ** (hi - e))
        return tot
    PA = clear(Acoef)
    QB = clear(Bcoef)
    mod = mpz(1) << Nbig
    G = (PA * theta - QB) % mod
    val = v2(G)
    pred = L * (-lo) + L * P(n + m + 1)
    lh = max(log2abs(PA), log2abs(QB))
    return val, lh, pred, lo, hi


def section_L3():
    hdr("L3  Theorem Q: Tschakaloff-Zudilin forms for quadratic position families (exact)")
    families = [
        ("squares k^2", 2, 1),               # P(l) = l^2 : a2 = 2 (alpha=1), s = alpha+beta = 1
        ("triangular k(k+1)/2", 1, 1),       # P(l) = l(l+1)/2 : alpha=1/2, beta=1/2
        ("pentagonal k(3k-1)/2", 3, 1),      # alpha=3/2, beta=-1/2
        ("pentagonal k(3k+1)/2", 3, 2),      # alpha=3/2, beta=+1/2
        ("2k^2+k", 4, 3),                    # alpha=2, beta=1
    ]
    maps = [
        ("3x+1, block 1^9 0 (Y)", 10, 3 ** 9),
        ("5x+1, block 1^6 0^4", 10, 5 ** 6),
        ("5x+1, block 1^7 0^3", 10, 5 ** 7),
        ("5x+1, block 1^8 0^2", 10, 5 ** 8),
    ]
    nmax = 8 if QUICK else 12
    for (mname, L, M0) in maps:
        mu = mubar(L, M0)
        print(f"\n-- map {mname}: rho = 2^{L}/{M0}, mu_bar = {mu:.4f} ({'< phi' if mu < PHI else '>= phi'})")
        for (fname, a2, s) in families:
            P = lambda l, a2=a2, s=s: a2 * l * (l - 1) // 2 + s * l
            # sanity: family values
            vals = [P(l) for l in range(6)]
            Nbig = 0
            # generous precision: predicted v2 at n=nmax, m=round(nmax/phi)
            mm = round(nmax / PHI)
            Nbig = L * (a2 * nmax * (nmax + mm) + (s + a2) * nmax + P(nmax + mm + 1)) + 4096
            theta = theta_mod(Nbig, P, L, M0)
            rows = []
            for n in range(1, nmax + 1):
                for m in (0, round(n / PHI)):
                    val, lh, pred, lo, hi = tschakaloff_form(n, m, a2, s, L, M0, Nbig, theta)
                    rows.append((n, m, val, pred, lh, val - lh))
            ok = all(r[2] == r[3] for r in rows)
            m0 = [r for r in rows if r[1] == 0]
            ms = [r for r in rows if r[1] == round(r[0] / PHI) and r[0] > 0]
            print(f"   {fname:24s} P(0..5)={vals}  v2==predicted: {ok}")
            print("      m=0     margins v2 - log2 H, n=1..: " + " ".join(f"{r[5]:+.0f}" for r in m0))
            print("      m=n/phi margins v2 - log2 H, n=1..: " + " ".join(f"{r[5]:+.0f}" for r in ms))


def section_L3b():
    hdr("L3b Theorem Q forms just above phi: squares, 5x+1 block 1^7 0^3 (mu_bar = 1.6253), larger n")
    L, M0, a2, s = 10, 5 ** 7, 2, 1
    P = lambda l: l * l
    nmax = 16 if QUICK else 40
    mm = round(nmax / PHI) + 2
    Nbig = L * (a2 * nmax * (nmax + mm) + (s + a2) * nmax + P(nmax + mm + 1)) + 4096
    theta = theta_mod(Nbig, P, L, M0)
    line = []
    for n in range(4, nmax + 1, 4):
        best = None
        for m in range(max(0, round(n / PHI) - 2), round(n / PHI) + 3):
            val, lh, pred, lo, hi = tschakaloff_form(n, m, a2, s, L, M0, Nbig, theta)
            assert val == pred
            if best is None or val - lh > best[0]:
                best = (val - lh, m)
        line.append(f"n={n}: {best[0]:+.0f} (m={best[1]})")
    print("   best margin over m in [n/phi-2, n/phi+2]:  " + "; ".join(line))
    x = 1 / PHI
    lead = (x * x + 4 * x + 3) - mubar(L, M0) * (x * x + 2 * x + 2)
    print(f"   leading coefficient (alpha L n^2) * [(x^2+4x+3) - mu_bar (x^2+2x+2)] at x = 1/phi: {lead:+.5f} < 0:")
    print("   the margins are eventually negative, but only after the O(n) terms are overtaken.")


# ----------------------------------------------------------------------------- L4
def section_L4():
    hdr("L4  cubic analogue of the Zudilin construction (exact, tiny n)")
    print("R(t) = sum_s c_s rho^(3ts(t+s)) vanishing at t = 0..n-1 (c = kernel of the n x (n+1) matrix);")
    print("Lambda = sum_s c_s rho^(-s^3) (X - S_(s-1)) = p X - r = sum_{t>=n} R(t) rho^(t^3).")
    rho = Fraction(2 ** 10, 3 ** 9)
    nmax = 4 if QUICK else 6
    Nbig = 200000
    X = X_mod(Nbig)
    mod = mpz(1) << Nbig
    print(" n | v2(P X - R) | log2 max(|P|,|R|) | exponent v2/log2H | (Liouville K=n: v2/log2H)")
    for n in range(1, nmax + 1):
        # matrix over Q
        M = [[rho ** (3 * t * s * (t + s)) for s in range(n + 1)] for t in range(n)]
        # kernel via cofactors: c_s = (-1)^s det(M without column s)
        def det(A):
            A = [row[:] for row in A]
            k = len(A)
            d = Fraction(1)
            for i in range(k):
                piv = next(r for r in range(i, k) if A[r][i] != 0)
                if piv != i:
                    A[i], A[piv] = A[piv], A[i]
                    d = -d
                d *= A[i][i]
                for r in range(i + 1, k):
                    f = A[r][i] / A[i][i]
                    if f:
                        for cidx in range(i, k):
                            A[r][cidx] -= f * A[i][cidx]
            return d
        c = []
        for s in range(n + 1):
            sub = [[M[t][j] for j in range(n + 1) if j != s] for t in range(n)]
            c.append((-1) ** s * det(sub))
        # check vanishing
        for t in range(n):
            assert sum(c[s] * M[t][s] for s in range(n + 1)) == 0
        p = sum(c[s] * rho ** (-(s ** 3)) for s in range(n + 1))
        r = Fraction(0)
        for s in range(n + 1):
            Ssm1 = sum((rho ** (k ** 3) for k in range(0, s)), Fraction(0))
            r += c[s] * rho ** (-(s ** 3)) * Ssm1
        den = math.lcm(p.denominator, r.denominator)
        P = mpz(p.numerator * (den // p.denominator))
        R = mpz(r.numerator * (den // r.denominator))
        g = gmpy2.gcd(P, R)
        P //= g
        R //= g
        val = v2((P * X - R) % mod)
        lh = max(log2abs(P), log2abs(R))
        A, D = partial_sum_fraction(n)
        lvl = v2((D * X - A) % mod)
        llh = max(log2abs(A), log2abs(D))
        print(f"{n:2d} | {val:11d} | {lh:17.1f} | {val/lh:17.4f} | {lvl/llh:.4f}")
    print("The kernel coefficients c_s carry the full tropical spread of rho^(3ts(t+s)); the construction")
    print("is worse than the partial sums at every n (no q-binomial cancellation exists for f(t,s) = ts(t+s)).")


# ----------------------------------------------------------------------------- L5
def section_L5():
    hdr("L5  threshold map")
    print("Square-swap words (Theorem Y / Q) need mu_bar = A log2 m / L < phi = 1.6180.")
    for m in (3, 5, 7, 9):
        crit = 1 / math.log2(m)                      # supercritical: A/L > log_m 2
        lim = PHI / math.log2(m)                    # Lemma P: A/L < phi/log2 m
        print(f"  {m}x+1: supercritical iff A/L > {crit:.4f}; Lemma P covers A/L < {min(lim,1):.4f}"
              f"{'  (all supercritical blocks covered)' if lim >= 1 else ''}")
    print("\nCheapest open square-swap instances (mu_bar in [phi, phi+0.01), L <= 40):")
    for m in (5, 7):
        best = []
        for L in range(1, 41):
            for A in range(1, L + 1):
                mu = A * math.log2(m) / L
                if PHI <= mu < PHI + 0.01 and math.gcd(A, L) == 1:
                    best.append((mu - PHI, A, L, mu))
        best.sort()
        for (_, A, L, mu) in best[:4]:
            print(f"   {m}x+1: A = {A:2d} ones in L = {L:2d}: mu_bar = {mu:.5f}  (phi + {mu-PHI:.5f})")
    print("   e.g. 5x+1 with block 1^7 0^3 has mu_bar = 1.62535 (checked in L3: margins turn negative).")
    print("\nCube-swap words: every supercritical block has mu_bar > 1, and every method on record needs mu_bar <= 1.")
    print("Nearest-to-critical supercritical 3x+1 blocks (L <= 60), mu_bar - 1:")
    rows = []
    for L in range(1, 61):
        A = math.floor(L / math.log2(3)) + 1          # least A with A/L > log_3 2
        if A <= L:
            rows.append((A * LOG2_3 / L - 1, A, L))
    rows.sort()
    for (d, A, L) in rows[:6]:
        print(f"   A = {A:2d}, L = {L:2d}: mu_bar = 1 + {d:.6f}")


# ----------------------------------------------------------------------------- L6
def section_L6():
    hdr("L6  near-critical cube swap: reach of the partial-sum certificates")
    # A=12, L=19 under 3x+1: rho = 2^19/3^12
    for (A, L) in ((12, 19), (9, 10)):
        M0 = 3 ** A
        mu = A * LOG2_3 / L
        best = None
        # gain(K) = L (K+1)^3 - log2 max(M0^{K^3}, |A_K|) - O(log);  use exact dominant terms:
        # |A_K| <= (K+1) max(2^L, M0)^{K^3}  -> log2 H <= mu L K^3 + log2(K+1)
        for K in range(1, 20000):
            g = L * (K + 1) ** 3 - mu * L * K ** 3 - math.log2(K + 1)
            if best is None or g > best[0]:
                best = (g, K)
        print(f"  A={A}, L={L}: mu_bar = {mu:.6f}; max_K gain = {best[0]:.4g} bits at K = {best[1]}"
              f" (every gain < 0 beyond K ~ {3/(mu-1):.0f})")
    print("  So near-critical cube-swap numbers have huge FINITE certificates (no rational of height")
    print("  <= 2^(max gain)) from their partial sums alone, yet no proof: the gain is eventually negative.")
    # exact check for A=12, L=19 at small K
    A, L = 12, 19
    M0 = 3 ** A
    Kmax = 6 if QUICK else 9
    Nbig = L * (Kmax + 2) ** 3 + 64
    Xn = theta_mod(Nbig, lambda k: k ** 3, L, M0)
    mod = mpz(1) << Nbig
    out = []
    for K in range(1, Kmax + 1):
        AK, D = partial_sum_fraction(K, L=L, M0=M0)
        val = v2((D * Xn - AK) % mod)
        lh = max(log2abs(AK), log2abs(D))
        out.append(f"K={K}: v2={val} (pred {L*(K+1)**3}), gain={val-lh:.1f}")
    print("  exact (A=12, L=19): " + "; ".join(out))


# ----------------------------------------------------------------------------- L7
def section_L7():
    hdr("L7  repetition profile of swap words at block level: squares, floor(k^(3/2)), cubes")
    import numpy as np
    print("LPF(j) = max_{a<j} LCE(a,j) of the swap indicator (block level; letters = 10 x blocks up to O(1)).")
    print("Two sources: (i) inside one gap (period 1 block): LPF >= distance to the next swap;")
    print("(ii) runs of equal consecutive gaps (matched swap patterns), extended by one partial gap.")
    print("Dio(w) - 1 = limsup LPF(j)/j; the table gives max LPF(j)/j over dyadic ranges of j.")

    def positions(kind, nmax):
        out, k = [], 1
        while True:
            v = {"square": k * k, "k32": math.isqrt(k * k * k), "cube": k ** 3}[kind]
            if v >= nmax:
                break
            out.append(v)
            k += 1
        return np.array(sorted(set(out)), dtype=np.int64)

    for kind, nmax in (("square", 4 * 10 ** 6), ("k32", 4 * 10 ** 6 if not QUICK else 4 * 10 ** 5), ("cube", 4 * 10 ** 9)):
        pos = positions(kind, nmax)
        g = np.diff(pos)
        best = {}
        def upd(j, l):
            if j <= 0:
                return
            t = int(math.log2(j))
            r = l / j
            if r > best.get(t, 0.0):
                best[t] = r
        # (i) within-gap repetitions: j just after a swap, l = gap - 1
        for a_, gg in zip(pos[:-1], g):
            upd(int(a_) + 1, int(gg) - 1)
        # (ii) matched gap runs
        n = len(g)
        for d in range(1, n):
            eq = (g[:-d] == g[d:])
            if not eq.any():
                continue
            e = np.concatenate(([0], eq.astype(np.int8), [0]))
            starts = np.flatnonzero(np.diff(e) == 1)
            ends = np.flatnonzero(np.diff(e) == -1)
            for st, en in zip(starts, ends):
                j = int(pos[st + d])
                l = int(pos[en + d] - pos[st + d])
                nxt = [int(g[en]) if en < n else 0, int(g[en + d]) if en + d < n else 0]
                upd(j, l + min(nxt))
        ts = sorted(best)
        prof = " ".join(f"2^{t}:{best[t]:.3f}" for t in ts if t >= 6)
        tail = [t for t in ts if t >= 6][-7:]
        slope = (math.log2(best[tail[-1]]) - math.log2(best[tail[0]])) / (tail[-1] - tail[0])
        print(f"  {kind:7s} (blocks < {nmax:.0e}): {prof}")
        print(f"          fitted decay over the last dyadic ranges: LPF(j)/j ~ j^({slope:.2f})")
    print("All three profiles decay to 0: Dio = 1 for each, so Theorem D reaches none of them;")
    print("floor(k^(3/2)) is no easier than Y3 by this measure.")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_cube_theta_ledger.py" + (" --quick" if QUICK else ""))
    section_L0()
    section_L1()
    section_L2()
    section_L3()
    section_L3b()
    section_L4()
    section_L5()
    section_L6()
    section_L7()
    print(f"\n[ledger done in {time.time()-t0:.1f}s]")
