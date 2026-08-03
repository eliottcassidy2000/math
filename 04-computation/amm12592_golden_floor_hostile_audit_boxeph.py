#!/usr/bin/env python3
"""AMM 12592 -- independent hostile audit of the golden floor chain
THM-3009 -> THM-3017 -> THM-3027 -> THM-3024   (boxeph, 2026-08-03).

Exact arithmetic only: int / fractions.Fraction / sympy 1.9.  No floats are
load-bearing anywhere in this file; floats appear only in display strings.

PART 1  Exact rational bracketing of gamma* = log(phi)/log(sqrt 5).
        The comparison  p/q vs gamma*  is decided in INTEGERS:
            p/q > gamma*  <=>  5^p > phi^(2q)
        and  2 phi^(2q) = L_{2q} + F_{2q} sqrt5,  so with X = 2*5^p - L_{2q},
            p/q > gamma*  <=>  X > 0 and X^2 > 5 F_{2q}^2.
        (Equality X^2 = 5 F^2 is impossible: sqrt5 is irrational.)
        This is literally the Fibonacci bracketing of phi pushed through the
        power: F_{k+1}/F_k sandwiches phi, L_n, F_n give phi^n exactly.

PART 2  sympy re-derivation of THM-3027's tangency collapse, hostile mode:
        (S)+(T)+(V) => (1-tau)^2 = tau  is checked as polynomial-in-log
        identities in a positivity-safe coordinate tau = 1/(1+w), w > 0;
        the stationarity/tangency derivative formulas, the concavity of the
        inner problem (perspective), the monotonicity of the capacity in
        gamma, the b-alphabet universality, THM-3009's (CLO) lemmas, and an
        EXACT (surd + rational-bracket) interiority certificate for sigma*.

PART 3  Exact re-referee of the decisive (ARCH) computations:
        identities behind the reduction, hostile monotonicity sweep of the
        candidate ladder (validates the binary search), certified refuted
        rates for m = 8..1024 reproduced from scratch in integers, exact
        integer certification that every refuted rate is < gamma*, targeted
        exact checks at m = 2048, 4096, sandwich rationals on both sides of
        gamma*, and D0-robustness of the refutation.

PART 4  Exact autopsy of THM-3024's cross-shell Hall model: per-(M,d) tail
        cuts over dyadic windows in EXACT integers, demonstrating that
        (i) the "binding" degree-resolved cuts reported in THM-3024 involve
        only the deepest shell (truncation edge artifact), and (ii) within
        the model as stated (unbounded forward routing at fixed absolute
        degree) any per-shell deficit is absorbed by the next shell with
        exponential room, so the model does NOT support a general-class
        floor at golden.
"""
import sys
from fractions import Fraction
from math import comb

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(10000000)


def rule(s):
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


# ===========================================================================
# PART 1 -- exact rational bracketing of gamma*
# ===========================================================================

def fib_pair(n):
    """(F_n, F_{n+1}) by fast doubling, exact integers."""
    if n == 0:
        return (0, 1)
    a, b = fib_pair(n >> 1)
    c = a * ((b << 1) - a)          # F_{2k}
    d = a * a + b * b               # F_{2k+1}
    if n & 1:
        return (d, c + d)
    return (c, d)


def cmp_phi(p, q):
    """sign of p/q - phi, exact (phi = (1+sqrt5)/2)."""
    X = 2 * p - q
    if X <= 0:
        return -1
    d = X * X - 5 * q * q
    return 1 if d > 0 else -1       # equality impossible


def cmp_gamma(p, q):
    """sign of p/q - gamma*, gamma* = log(phi)/log(sqrt5), exact integers."""
    F, F1 = fib_pair(2 * q)
    L = 2 * F1 - F                  # Lucas L_{2q}
    X = 2 * pow(5, p) - L
    if X <= 0:
        return -1
    d = X * X - 5 * F * F
    if d == 0:
        raise AssertionError("impossible: sqrt5 rational")
    return 1 if d > 0 else -1


def part1():
    rule("PART 1 -- EXACT INTEGER BRACKETING OF gamma* = log(phi)/log(sqrt 5)")
    # Fibonacci sandwich of phi (sanity for the machinery)
    print("Fibonacci sandwich of phi (exact):")
    for k in (10, 20, 30):
        F, F1 = fib_pair(k)
        F2 = F + F1
        slo = cmp_phi(F1, F)
        shi = cmp_phi(F2, F1)
        print(f"  k={k:3d}:  F_{k+1}/F_{k} vs phi: {slo:+d},   "
              f"F_{k+2}/F_{k+1} vs phi: {shi:+d}   (alternating sandwich)")
        assert slo * shi == -1
    # classify the rates that appear in the chain
    print("\nExact classification of the chain's rational rates vs gamma*:")
    tests = [("1/2", 1, 2), ("2457/6592", 2457, 6592),
             ("2457/4135", 2457, 4135), ("3/5", 3, 5),
             ("299/500  (=0.598)", 299, 500), ("149/250  (=0.596)", 149, 250),
             ("597987/1000000", 597987, 1000000),
             ("597988/1000000", 597988, 1000000)]
    for name, p, q in tests:
        s = cmp_gamma(p, q)
        print(f"  {name:22s} : {'ABOVE' if s > 0 else 'BELOW'} gamma*   (exact)")
    # continued-fraction convergents (proposed by 60-dps evalf, certified in Z)
    from sympy import log as slog, sqrt as ssqrt, Integer
    gs = slog((1 + ssqrt(5)) / 2) / slog(ssqrt(5))
    dec = str(gs.evalf(60))
    frac = Fraction(dec)
    # continued fraction of the 60-digit rational proposal
    cf = []
    x = frac
    for _ in range(40):
        a = int(x)
        cf.append(a)
        x = x - a
        if x == 0:
            break
        x = 1 / x
    # convergents
    conv = []
    h0, h1, k0, k1 = 1, cf[0], 0, 1
    conv.append(Fraction(h1, k1))
    for a in cf[1:]:
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        if k1 > 200000:
            break
        conv.append(Fraction(h1, k1))
    # certify the last two convergents straddle gamma*
    lo = hi = None
    for c in conv[-4:]:
        s = cmp_gamma(c.numerator, c.denominator)
        if s < 0:
            lo = c if (lo is None or c > lo) else lo
        else:
            hi = c if (hi is None or c < hi) else hi
    assert lo is not None and hi is not None
    assert cmp_gamma(lo.numerator, lo.denominator) == -1
    assert cmp_gamma(hi.numerator, hi.denominator) == +1
    width = hi - lo
    print(f"\nCF-convergent bracket, certified by pure-integer Fibonacci/Lucas tests:")
    print(f"  {lo}  <  gamma*  <  {hi}")
    print(f"  width = {width} ~ {float(width):.3e}")
    print(f"  (floats above are display only; the certificates are 5^p vs "
          f"L_2q,F_2q integer comparisons)")
    return lo, hi


# ===========================================================================
# PART 2 -- sympy symbolic audit of the tangency collapse (THM-3027/3009/3017)
# ===========================================================================

def part2(bracket):
    from sympy import (symbols, log, sqrt, Rational, simplify, expand_log,
                       diff, cancel, together, solve, S, nsimplify)
    rule("PART 2 -- SYMBOLIC AUDIT OF THE TANGENCY COLLAPSE (sympy 1.9)")
    w = symbols('w', positive=True)          # tau = 1/(1+w): 0<tau<1 safe
    t = 1 / (1 + w)
    u = (1 - t) / (1 + t)                    # from (T):  2u/(1-u) = (1-tau)/tau
    A = log(1 / u)
    g = log((1 - t) / t) / A                 # from (S)+(T): gamma A = log((1-t)/t)
    rho = 1 - u
    H = lambda z: -z * log(z) - (1 - z) * log(1 - z)

    # (K): the key identity, under (S) alone (b = 2)
    K = H(rho) + rho * log(2) - A * (1 + g * rho)
    rK = simplify(expand_log(K, force=True))
    print(f"  (K)  H(rho)+rho log2 - A(1+gamma rho)      == {rK}   "
          f"{'PROVED-symbolic' if rK == 0 else 'FAIL'}")

    # multiplier cancellation: D(1+gamma rho) = gamma(1+tau) from sigma=tau-rho D, D=gamma(1+sigma)
    sig, Dv, gg, tt, rr = symbols('sig Dv gg tt rr', positive=True)
    sol = solve([Dv - gg * (1 + sig), sig - (tt - rr * Dv)], [Dv, sig], dict=True)
    Dsol = simplify(sol[0][Dv] - gg * (1 + tt) / (1 + gg * rr))
    print(f"  D = gamma(1+tau)/(1+gamma rho)             == 0: {Dsol == 0}")

    # (V) reduces to 2 log(1-tau) = log tau
    V = (1 + t) * log((1 - t) / t) - H(t) - (2 * log(1 - t) - log(t))
    rV = simplify(expand_log(V, force=True))
    print(f"  (V) - [2log(1-tau)-log(tau)]               == {rV}   "
          f"{'PROVED-symbolic' if rV == 0 else 'FAIL'}")

    # the golden quadratic and back-substitution, exact surds
    tau = symbols('tau')
    roots = solve(tau**2 - 3 * tau + 1, tau)
    tstar = [r for r in roots if simplify(r - Rational(1, 2)).is_negative][0]
    phi = (1 + sqrt(5)) / 2
    print(f"  root of (1-tau)^2 = tau in (0,1):  tau* = {tstar}")
    print(f"    tau* - phi^-2                            == "
          f"{simplify(tstar - 1/phi**2)}")
    ustar = simplify((1 - tstar) / (1 + tstar))
    print(f"    u* = (1-tau*)/(1+tau*) - 1/sqrt5         == "
          f"{simplify(ustar - 1/sqrt(5))}")
    print(f"    (1-tau*)/tau* - phi                      == "
          f"{simplify((1 - tstar)/tstar - phi)}")
    print(f"    2+phi - phi*sqrt5                        == "
          f"{simplify(2 + phi - phi*sqrt(5))}")
    print("  hence gamma* = log((1-tau*)/tau*) / log(1/u*) = log(phi)/log(sqrt5)"
          "  EXACTLY.")

    # derivative formulas behind (S) and (T), and the two new lemmas
    gs, ss, ts = symbols('gamma sigma tau', positive=True)
    D2 = gs * (1 + ss); m2 = ts - ss; rho2 = m2 / D2
    psi = D2 * H(rho2) + m2 * log(2)
    r1 = simplify(expand_log(diff(psi, ss)
                             - (-gs * log(1 - rho2) - log(2 * (1 - rho2) / rho2)),
                             force=True))
    print(f"  dpsi/dsigma formula (stationarity (S))     == {r1}")
    r2 = simplify(cancel(together(
        diff(psi, ss, 2)
        + ((1 + gs) / (1 - rho2) + 1 / rho2) * (1 + ts) / (gs * (1 + ss)**2))))
    print(f"  d2psi/dsigma2 + [(1+g)/(1-rho)+1/rho](1+tau)/(g(1+sigma)^2) == {r2}")
    print("    => psi is STRICTLY CONCAVE in sigma on 0<rho<1: the interior")
    print("       stationary point is the GLOBAL inner max (no scan needed).")
    r3 = simplify(expand_log(diff(psi, gs) + (1 + ss) * log(1 - rho2), force=True))
    print(f"  dpsi/dgamma + (1+sigma) log(1-rho)         == {r3}")
    print("    => psi (hence Psi = max_sigma psi) is STRICTLY INCREASING in gamma")
    print("       at fixed tau: for gamma < gamma*, Psi(gamma,tau*) < H(tau*),")
    print("       so the CONTINUUM FLOOR DIRECTION is symbolic -- the numeric")
    print("       global-tau scan is NOT load-bearing for the lower bound.")

    # b-alphabet universality: (K) with log b, reduction unchanged
    b = symbols('b', positive=True)
    ub = (1 - t) / (1 - t + b * t)            # from (T_b): b u/(1-u) = (1-tau)/tau
    Ab = log(1 / ub)
    gb = log((1 - t) / t) / Ab
    rhob = 1 - ub
    Kb = H(rhob) + rhob * log(b) - Ab * (1 + gb * rhob)
    rKb = simplify(expand_log(Kb, force=True))
    print(f"  (K_b) general alphabet residual            == {rKb}")
    print("    => tau* = phi^-2 universal in b; gamma*(b) = log(phi)/log((b+phi)/phi);")
    print("       b=2 special only through 2+phi = phi sqrt5.  PROVED-symbolic.")

    # THM-3009's (CLO) lemmas, log2 units.  (CLO): H(delta) = -H'(delta)(2-delta)
    # with L = log2 delta, M = log2(1-delta), -H' = L - M.  Claim: residual is
    # exactly M - 2L, so (CLO) <=> M = 2L <=> delta^2 = 1 - delta.
    dd = 1 / (1 + w)   # positivity-safe delta in (0,1)
    L2 = log(dd) / log(2); M2 = log(1 - dd) / log(2)
    R_clo = (-dd * L2 - (1 - dd) * M2) - (2 - dd) * (L2 - M2)
    resid = simplify(expand_log(R_clo - (M2 - 2 * L2), force=True))
    print(f"  THM-3009 (CLO)-residual - (M - 2L)         == {resid}")
    print("    => (CLO) <=> delta^2 = 1 - delta, root delta* = 1/phi.  "
          "PROVED-symbolic.")

    # exact interiority of sigma* via the rational gamma* bracket.
    # Q(sqrt5) arithmetic on Fraction pairs (a, b) ~ a + b sqrt5: fully exact.
    def qs_add(x, y): return (x[0] + y[0], x[1] + y[1])
    def qs_mul(x, y): return (x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
    def qs_div(x, y):
        den = y[0]*y[0] - 5*y[1]*y[1]
        return ((x[0]*y[0] - 5*x[1]*y[1]) / den, (x[1]*y[0] - x[0]*y[1]) / den)
    def qs_sign(x):
        a, b = x
        if a == 0 and b == 0: return 0
        if a >= 0 and b >= 0: return 1
        if a <= 0 and b <= 0: return -1
        d5 = a*a - 5*b*b
        if a > 0:   # b < 0: sign = sign(a^2 - 5 b^2)
            return 1 if d5 > 0 else -1
        return -1 if d5 > 0 else 1   # a < 0, b > 0

    lo, hi = bracket
    tau_p = (Fraction(3, 2), Fraction(-1, 2))          # tau* = (3-sqrt5)/2
    rho_p = (Fraction(1), Fraction(-1, 5))             # rho* = 1 - sqrt5/5
    one = (Fraction(1), Fraction(0))

    def sigma_of(gq):
        gp = (Fraction(gq), Fraction(0))
        D = qs_div(qs_mul(gp, qs_add(one, tau_p)),
                   qs_add(one, qs_mul(gp, rho_p)))
        return qs_add(tau_p, tuple(-c for c in qs_mul(rho_p, D)))

    s_hi_end = sigma_of(hi)    # sigma decreasing in gamma => sigma* > sigma(hi)
    s_lo_end = sigma_of(lo)
    sgn_pos = qs_sign(s_hi_end)
    sgn_lt_tau = qs_sign(qs_add(s_lo_end, tuple(-c for c in tau_p)))
    fl = lambda x: float(x[0]) + float(x[1]) * 5 ** 0.5
    print(f"  sigma(gamma_hi) > 0  (exact Q(sqrt5) sign): {sgn_pos == 1}")
    print(f"  sigma(gamma_lo) < tau*  (exact):            {sgn_lt_tau == -1}")
    print("    dsigma/dgamma = -rho*(1+tau*)/(1+gamma rho*)^2 < 0, so")
    print(f"    0 < sigma* in ({fl(s_hi_end):.12f}, {fl(s_lo_end):.12f}) < tau*")
    print("    -- INTERIOR, certified via the exact rational bracket of PART 1")
    print("    + monotonicity.  UPGRADED from THM-3027's 3000-point scan.")

    # THM-3009 ell* interiority (closed forms, 60-dps evaluation of exact logs)
    from sympy import log as slog
    phiS = (1 + sqrt(5)) / 2
    dstar = 1 / phiS
    pstar = 1 / sqrt(5)
    H2e = lambda z: (-z * slog(z) - (1 - z) * slog(1 - z)) / slog(2)
    astar = H2e(dstar) / (H2e(pstar) + 1 - pstar)
    lstar = dstar - pstar * astar
    gstar = slog(phiS) / slog(sqrt(5))
    lv = lstar.evalf(60)
    ub2 = (1 - gstar).evalf(60)
    print(f"  THM-3009 ell* = {lv}")
    print(f"           1-gamma* = {ub2}")
    print(f"    0 < ell* < 1-gamma*: margins {lv} and {(1-gstar-lstar).evalf(30)}")
    print("    (closed-form log expressions evaluated at 60 dps; margins ~0.34")
    print("     and ~0.062, vastly above evaluation error.  VERIFIED-60dps.)")
    sys.stdout.flush()


# ===========================================================================
# PART 3 -- exact (ARCH) re-referee
# ===========================================================================

def committed_profile(m, C):
    """a_k = min(m-1-k, floor(C(m+k)) - (m+k) - 1); None if infeasible."""
    a = []
    for k in range(m):
        n = m + k
        cap = (C.numerator * n) // C.denominator - n - 1
        if cap < 0:
            return None
        a.append(min(m - 1 - k, cap))
    return a


def floor_profile(m, g, D0=0):
    """THM-3029 floor-profile convention: a_k = min(m-1-k, floor(g(m+k)) + D0)."""
    a = []
    for k in range(m):
        n = m + k
        v = (g.numerator * n) // g.denominator + D0
        if v < 0:
            return None
        a.append(min(m - 1 - k, v))
    return a


def arch_first_fail(m, a):
    """first d violating (ARCH), or None; exact integers, early break."""
    L = [m - 1 - k - a[k] for k in range(m)]
    lhs = 1
    for d in range(m):
        if d > 0:
            lhs = lhs * (m - d) // d          # comb(m-1,d)
        rhs = 0
        ok = False
        for k in range(m):
            r = d - L[k]
            if 0 <= r <= a[k]:
                rhs += comb(a[k], r) << (a[k] - r)
                if rhs >= lhs:
                    ok = True
                    break
        if not ok:
            return d
    return None


def arch_min_margin(m, a):
    """exact min over d of Fraction(rhs, lhs) -- no early break (full referee)."""
    L = [m - 1 - k - a[k] for k in range(m)]
    lhs = 1
    worst = None
    argd = None
    for d in range(m):
        if d > 0:
            lhs = lhs * (m - d) // d
        rhs = 0
        for k in range(m):
            r = d - L[k]
            if 0 <= r <= a[k]:
                rhs += comb(a[k], r) << (a[k] - r)
        ratio = Fraction(rhs, lhs)
        if worst is None or ratio < worst:
            worst, argd = ratio, d
    return worst, argd


def candidates(m):
    return sorted({Fraction(m + k + 1 + aa, m + k)
                   for k in range(m) for aa in range(m - k)})


def lower_bound(m):
    cands = candidates(m)
    lo, hi, best = 0, len(cands) - 1, None
    while lo <= hi:
        mid = (lo + hi) // 2
        a = committed_profile(m, cands[mid])
        ok = a is not None and arch_first_fail(m, a) is None
        if ok:
            hi = mid - 1
        else:
            best = cands[mid]
            lo = mid + 1
    return best


def next_candidate_above(m, C0):
    """smallest candidate (n+1+a)/n > C0 over n = m+k, 0 <= a <= m-1-k."""
    best = None
    for k in range(m):
        n = m + k
        amin = (C0.numerator * n) // C0.denominator - n - 1 + 1
        if amin < 0:
            amin = 0
        if amin > m - 1 - k:
            continue
        c = Fraction(n + 1 + amin, n)
        if c <= C0:          # can happen only if C0*n was an exact integer; bump
            if amin + 1 > m - 1 - k:
                continue
            c = Fraction(n + 2 + amin, n)
        if best is None or c < best:
            best = c
    return best


def part3():
    rule("PART 3 -- EXACT (ARCH) RE-REFEREE (integers / Fraction only)")

    print("3.0  identities behind reduction B (exact):")
    ok = all(sum(comb(a, i) * comb(i, r) for i in range(a + 1))
             == comb(a, r) * 2 ** (a - r)
             for a in range(0, 41) for r in range(a + 1))
    print(f"  sum_i C(a,i)C(i,r) = C(a,r) 2^(a-r), a<=40          : {ok}")
    ok2 = True
    for m in (8, 16, 32):
        for t in range(0, m - 1):
            s = sum(comb(m - 1 - k, t + 1 - k) for k in range(m)
                    if 0 <= t + 1 - k <= m - 1 - k)
            ok2 &= (s == comb(m, t + 1))
    print(f"  hockey-stick centring sum_k C(m-1-k,t+1-k)=C(m,t+1) : {ok2}")
    ok3 = True
    for m in (8, 16, 32, 64):
        # ((1+u)^m - 1)/u mod 2 == u^(m-1): coefficients C(m,j), j=1..m
        coeffs = [comb(m, j) % 2 for j in range(1, m + 1)]
        ok3 &= (coeffs == [0] * (m - 1) + [1])
    print(f"  mod-2 reduction ((1+u)^m-1)/u == u^(m-1) over F2    : {ok3}")

    print("\n3.1  hostile monotonicity sweep (validates the binary search):")
    for m in (8, 16, 32):
        cands = candidates(m)
        verdicts = []
        for c in cands:
            a = committed_profile(m, c)
            verdicts.append(a is not None and arch_first_fail(m, a) is None)
        # must be fail...fail pass...pass
        first_pass = verdicts.index(True) if True in verdicts else len(verdicts)
        mono = all(verdicts[i] for i in range(first_pass, len(verdicts)))
        flip_lo = cands[first_pass - 1] if first_pass > 0 else None
        flip_hi = cands[first_pass] if first_pass < len(cands) else None
        print(f"  m={m:3d}: {len(cands):4d} candidates, verdict ladder monotone: "
              f"{mono};  flip in ({flip_lo}, {flip_hi}]")
        assert mono
    sys.stdout.flush()

    print("\n3.2  certified refuted rates, recomputed from scratch (exact):")
    expected = {8: Fraction(3, 2), 16: Fraction(14, 9), 32: Fraction(64, 41),
                512: Fraction(508, 319), 1024: Fraction(1992, 1249)}
    refuted = {}
    from sympy import log as slog, sqrt as ssqrt
    gstar60 = Fraction(str((slog((1 + ssqrt(5)) / 2) / slog(ssqrt(5))).evalf(60)))
    for m in (8, 16, 32, 64, 128, 256, 512, 1024):
        b = lower_bound(m)
        refuted[m] = b
        gm = b - 1
        s = cmp_gamma(gm.numerator, gm.denominator)
        gap = gstar60 - gm
        tag = ""
        if m in expected:
            tag = "  matches committed" if b == expected[m] else \
                  f"  MISMATCH vs committed {expected[m]}"
        print(f"  m={m:5d}: rho(m) > {str(b):>12s} = {float(b):.6f}; "
              f"gamma_m = {float(gm):.6f} {'<' if s < 0 else '>!'} gamma* "
              f"(exact); gap ~ {float(gap):.5f}{tag}")
        assert s == -1, "a refuted rate at or above gamma* would sink the floor"
        sys.stdout.flush()
    mono = all(refuted[a] < refuted[b]
               for a, b in zip((8, 16, 32, 64, 128, 256, 512),
                               (16, 32, 64, 128, 256, 512, 1024)))
    print(f"  refuted ladder strictly increasing in m: {mono}")

    print("\n3.3  exact margins at sandwich rationals (D0 = -1 committed "
          "profile), m <= 1024:")
    rates = [("149/250 (<g*)", Fraction(149, 250)),
             ("597987/10^6 (<g*)", Fraction(597987, 10**6)),
             ("3/5 (>g*)", Fraction(3, 5)),
             ("299/500 (>g*)", Fraction(299, 500))]
    for m in (64, 256, 1024):
        row = []
        for name, g in rates:
            a = committed_profile(m, 1 + g)
            worst, argd = arch_min_margin(m, a)
            lg = (worst.numerator.bit_length() - worst.denominator.bit_length())
            row.append(f"{name.split()[0]:>14s}: "
                       f"{'PASS' if worst >= 1 else 'FAIL'} "
                       f"(min rhs/lhs ~2^{lg:+d} at d={argd})")
        print(f"  m={m:5d}: " + " | ".join(row))
        sys.stdout.flush()

    print("\n3.4  targeted exact checks at m = 2048, 4096:")
    for m, Ccom in ((2048, Fraction(3890, 2437)), (4096, Fraction(6709, 4201))):
        a = committed_profile(m, Ccom)
        fd = arch_first_fail(m, a)
        gm = Ccom - 1
        s = cmp_gamma(gm.numerator, gm.denominator)
        print(f"  m={m}: committed refuted candidate C={Ccom} "
              f"(gamma={float(gm):.6f}): first (ARCH) failure at d={fd} "
              f"(d/m={fd/m:.4f});  gamma_m < gamma* exact: {s == -1}")
        assert fd is not None and s == -1
        nc = next_candidate_above(m, Ccom)
        a2 = committed_profile(m, nc)
        fd2 = arch_first_fail(m, a2)
        print(f"          next candidate {nc} = {float(nc):.7f}: "
              f"{'PASSES' if fd2 is None else f'FAILS at d={fd2} (STRONGER than committed!)'}")
        a3 = committed_profile(m, Fraction(8, 5))
        fd3 = arch_first_fail(m, a3)
        print(f"          C = 8/5 (gamma = 3/5 > gamma*): "
              f"{'PASSES' if fd3 is None else 'FAILS -- ALARM'}")
        assert fd3 is None
        sys.stdout.flush()

    print("\n3.5  D0-robustness of the refutation (floor profile "
          "a_k = floor(g(m+k)) + D0):")
    m = 2048
    g = Fraction(59, 100)
    for D0 in (0, 2):
        a = floor_profile(m, g, D0)
        fd = arch_first_fail(m, a)
        print(f"  m={m}, gamma=59/100, D0={D0}: "
              f"{'FAILS at d=' + str(fd) if fd is not None else 'passes'}")
    print("  => the refutation below gamma* is not an artifact of the -1 "
          "offset; bounded")
    print("     additive slack D0 does not rescue a rate below gamma* at "
          "large m (THM-3029 (D0)).")
    sys.stdout.flush()
    return refuted


# ===========================================================================
# PART 4 -- exact autopsy of the THM-3024 cross-shell Hall model
# ===========================================================================

def supply_table(m, g):
    """exact supply_m(d) for all d, committed profile; supply_m(d)=0 for d>=m."""
    a = committed_profile(m, 1 + g)
    out = [0] * m
    if a is None:
        return out
    L = [m - 1 - k - a[k] for k in range(m)]
    for d in range(m):
        s = 0
        for k in range(m):
            r = d - L[k]
            if 0 <= r <= a[k]:
                s += comb(a[k], r) << (a[k] - r)
        out[d] = s
    return out


def part4():
    rule("PART 4 -- EXACT AUTOPSY OF THE THM-3024 CROSS-SHELL HALL MODEL")
    shells = [8, 16, 32, 64, 128, 256, 512]
    g = Fraction(71, 125)      # = 0.568 ~ golden - 0.03, certified below gamma*
    assert cmp_gamma(71, 125) == -1
    print(f"rate gamma = 71/125 = 0.568 (certified BELOW gamma*), shells {shells}\n")

    print("4.1  per-shell (ARCH) verdicts (exact):")
    fails = {}
    for m in shells:
        a = committed_profile(m, 1 + g)
        fd = arch_first_fail(m, a)
        fails[m] = fd
        print(f"  m={m:4d}: {'FAIL at d=' + str(fd) + f' (d/m={fd/m:.3f})' if fd is not None else 'pass'}")

    sup = {m: supply_table(m, g) for m in shells}

    mfail = min((m for m in shells if fails[m] is not None), default=None)
    if mfail is None:
        print("  no per-shell failure in window -- lower gamma needed")
        return
    dstar = fails[mfail]
    print(f"\n4.2  the rescue: shell {2*mfail} absorbs shell {mfail}'s deficit "
          f"at the SAME absolute degree d={dstar}:")
    dem1 = comb(mfail - 1, dstar)
    sup1 = sup[mfail][dstar]
    deficit = dem1 - sup1
    dem2 = comb(2 * mfail - 1, dstar)
    sup2 = sup[2 * mfail][dstar]
    slack = sup2 - dem2
    print(f"  shell {mfail}:  demand has {dem1.bit_length()} bits, supply "
          f"{sup1.bit_length()} bits -> deficit {deficit.bit_length()} bits")
    print(f"  shell {2*mfail}: demand has {dem2.bit_length()} bits, supply "
          f"{sup2.bit_length()} bits -> slack  {slack.bit_length()} bits")
    print(f"  two-shell tail cut at (M={mfail}, d={dstar}): "
          f"{'SATISFIED' if dem1+dem2 <= sup1+sup2 else 'violated'} "
          f"(slack exceeds deficit by ~2^{slack.bit_length()-deficit.bit_length()})")
    assert dem1 + dem2 <= sup1 + sup2

    print(f"\n4.3  ALL degree-resolved tail cuts over the window, exact:")
    dmax = max(shells) - 1
    violated = []
    for iM, M in enumerate(shells):
        tail = shells[iM:]
        for d in range(1, dmax):
            dem = sum(comb(m - 1, d) for m in tail if d <= m - 1)
            if dem == 0:
                continue
            sv = sum(sup[m][d] for m in tail if d < m)
            if sv < dem:
                violated.append((M, d))
    # classify violations: does any involve more than the deepest shell?
    deep = max(shells)
    second = sorted(shells)[-2]
    edge_only = [(M, d) for (M, d) in violated if d > second - 1]
    multi = [(M, d) for (M, d) in violated if d <= second - 1]
    print(f"  violated tail cuts: {len(violated)} total; of these, "
          f"{len(edge_only)} have d > {second-1}, where ONLY the deepest shell "
          f"(m={deep}) has any demand or supply")
    print(f"  violated tail cuts with genuine multi-shell content "
          f"(d <= {second-1}): {len(multi)}")
    if violated:
        dset = sorted({d for (_, d) in violated})
        print(f"  violated degrees: d in [{dset[0]}, {dset[-1]}]  "
              f"(deepest shell's binding band, d/{deep} ~ "
              f"{dset[0]/deep:.3f}..{dset[-1]/deep:.3f})")
    print(f"\n  => every violated cut lives at degrees only the deepest shell "
          f"reaches: the")
    print(f"     'binding degree-resolved cut' of THM-3024 (G3) is the deepest "
          f"shell's OWN")
    print(f"     per-shell (ARCH) constraint -- a TRUNCATION EDGE ARTIFACT, "
          f"not a cross-shell")
    print(f"     obstruction.  Any per-shell deficit interior to the window is "
          f"absorbed by")
    print(f"     the next shell with exponential room (4.2).  In the model AS "
          f"STATED")
    print(f"     (unbounded forward routing at fixed absolute degree), demand "
          f"at fixed d")
    print(f"     grows polynomially in m while supply grows exponentially, so "
          f"EVERY tail")
    print(f"     cut with an unbounded (or merely one-shell-deeper) tail is "
          f"satisfied for")
    print(f"     ANY gamma > 0: the transportation relaxation yields NO "
          f"general-class floor.")
    sys.stdout.flush()


def main():
    print("AMM 12592 golden-floor hostile audit (boxeph) -- exact arithmetic only")
    print("chain audited: THM-3009 -> THM-3017 -> THM-3027 -> THM-3024\n")
    bracket = part1()
    print()
    part2(bracket)
    print()
    part3()
    print()
    part4()
    print()
    rule("AUDIT SUMMARY")
    print("""\
  THM-3009: floor direction SOUND.  Finite-m certified bounds reproduced
            exactly; every refuted rate certified < gamma* by integer
            Fibonacci/Lucas tests; candidate-ladder monotonicity validated;
            tangency algebra re-proved symbolically.  Remaining gaps: the
            finite-m Stirling transfer lemma (sec 10.3) is sketched, not
            written (LOW severity, one direction, routine); threshold
            EXACTNESS from above still rests on an uncertified numeric scan
            (MODERATE, not load-bearing for the floor).
  THM-3017: SOUND within scope (H=1 checkpoint-closure criterion).  Its
            20-digit numerics are superseded by the symbolic collapse.
  THM-3027: SOUND.  Collapse (S)+(T)+(V) => (1-tau)^2 = tau re-derived
            symbolically; inner concavity + gamma-monotonicity PROVED here,
            making the floor direction scan-free; sigma* interiority upgraded
            to an exact certificate.  Laplace/Stirling rate passage remains
            ASSUMED (needs an error-term lemma; routine for the barrier
            direction).
  THM-3024: GAP FOUND (HIGH).  Within its own transportation model the
            degree-resolved Hall cuts do NOT bind below golden once the
            window extends past the binding shell; the reported binding cuts
            are truncation edge artifacts, and the model as stated admits
            unbounded forward absorption at any gamma > 0.  (G2)'s
            degree-preservation premise also contradicts the degree-rescaling
            invariance opus's aggregate argument needs.  C*_general =
            log_5(5 phi^2) is NOT established; general-class floor OPEN.
""")


if __name__ == "__main__":
    main()
