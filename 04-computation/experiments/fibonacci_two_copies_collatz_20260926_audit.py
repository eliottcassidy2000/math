#!/usr/bin/env python3
"""fibonacci_two_copies_collatz_20260926_audit.py -- independent audit of
05-knowledge/results/fibonacci_two_copies_collatz_20260926.md (opus S12, session fibonacci-two-copies-20260926)
and its scripts fibonacci_two_copies_20260926.py / collatz_cycle_shadows_20260926.py (auditor: fable, 2026-09-27).

A. Real Binet function F(x) = (phi^x - cos(pi x) phi^-x)/sqrt5: the four identities checked SYMBOLICALLY (sympy)
   and at 50 digits (mpmath); the interpolant-family claim tested against explicit counterexamples; zeros; F(0).
B. Shadow proposition: exact-word residue classes mod 2^(A+1) by brute force (and THM-4512's E_w / C_w), fixed
   point x_w, the U^(mp) formula on both signs; negative odd cycles with |n| <= 10^6 by a full basin scan; the
   '120 primitive no-descent words up to rotation' recounted three ways (the script's rule, all primitive growth
   necklaces, Moebius) and the necklaces the script's rule misses; integer fixed points for p <= 12; the 27 and
   (4^k+17)/3 families; the D(k) dynamic programme against the S11 table and the periodic-word fraction; B7: what the
   codex lane's 'three parameterized rule nodes' and its 'decreasing parameter' family N(k,h) actually are.
C. NegaFibonacci: own implementation (interval greedy), existence/uniqueness by subset enumeration, class
   densities at 10^5 (and 10^6) for both signs, Zeckendorf classes, the unit-class intersection against phi^-4,
   exact Wythoff characterisations of both unit classes (these PROVE the independence), R_3 = floor(phi^2 Z),
   the 'shifted index set' claim.
Usage: python3 fibonacci_two_copies_collatz_20260926_audit.py > ../../05-knowledge/results/fibonacci_two_copies_collatz_20260926_audit.out
"""
import math, sys, time, random, bisect
from fractions import Fraction
from math import isqrt
import sympy as sp
import mpmath as mp

T0 = time.time()
LOG23 = math.log2(3)
PHI = (1 + 5 ** 0.5) / 2
FAILS = []


def hdr(s):
    print("\n== %s ==" % s)


def chk(name, ok, extra=""):
    print(" [%s] %s%s" % ("PASS" if ok else "FAIL", name, (": " + extra) if extra else ""))
    if not ok:
        FAILS.append(name)


def fib(n):
    if n >= 0:
        a, b = 0, 1
        for _ in range(n):
            a, b = b, a + b
        return a
    return (-1) ** (n + 1) * fib(-n)


# ----------------------------------------------------------------------------------------------------------
# A. the real Binet function
# ----------------------------------------------------------------------------------------------------------
def partA():
    hdr("A1. symbolic identities (sympy; a = phi^x, c = cos(pi x), s = sin(pi x); cos(pi(x+1)) = -c, phi^(x+1) = phi a)")
    phi = (1 + sp.sqrt(5)) / 2
    a, c, s = sp.symbols('a c s', positive=True)
    r5 = sp.sqrt(5)
    F0 = (a - c / a) / r5                              # F(x)
    F1 = (phi * a + c / (phi * a)) / r5                # F(x+1)
    F2 = (phi ** 2 * a - c / (phi ** 2 * a)) / r5      # F(x+2)
    Fm = (1 / a - c * a) / r5                          # F(-x)
    rec = sp.simplify(sp.expand(F2 - F1 - F0))
    chk("recurrence F(x+2) = F(x+1) + F(x) for all real x", rec == 0, "residual = %s" % rec)
    tw = sp.simplify(sp.expand((Fm - (-c * F0 + (1 / a) * s ** 2 / r5)).subs(s ** 2, 1 - c ** 2)))
    chk("twist F(-x) = -cos(pi x) F(x) + phi^-x sin^2(pi x)/sqrt5", tw == 0, "residual = %s" % tw)
    cas = sp.simplify(sp.expand(F1 ** 2 - F0 * F2 - c))
    chk("smooth Cassini F(x+1)^2 - F(x) F(x+2) = cos(pi x)", cas == 0, "residual = %s (uses phi^2 + phi^-2 = 3)" % cas)
    n = sp.symbols('n', integer=True)
    # F(-n) = (phi^-n - (-1)^n phi^n)/sqrt5 versus (-1)^(n+1) F_n = (-1)^(n+1)(phi^n - (-1)^n phi^-n)/sqrt5
    lhs = (phi ** (-n) - (-1) ** n * phi ** n) / r5
    rhs = (-1) ** (n + 1) * (phi ** n - (-1) ** n * phi ** (-n)) / r5
    chk("F(-n) = F_(-n) = (-1)^(n+1) F_n (symbolic in integer n)", sp.simplify(sp.expand(lhs - rhs)) == 0)

    hdr("A2. numeric spot checks at 50 digits")
    mp.mp.dps = 50
    ph = (1 + mp.sqrt(5)) / 2
    F = lambda x: (ph ** x - mp.cos(mp.pi * x) * ph ** (-x)) / mp.sqrt(5)
    e_int = max(abs(F(mp.mpf(k)) - fib(k)) for k in range(-60, 61))
    chk("F(n) = F_n for -60 <= n <= 60 (max error %.1e)" % e_int, e_int < mp.mpf(10) ** -30)
    rng = random.Random(7)
    xs = [mp.mpf(rng.uniform(-12, 12)) for _ in range(40)]
    e1 = max(abs(F(x + 2) - F(x + 1) - F(x)) for x in xs)
    e2 = max(abs(F(-x) - (-mp.cos(mp.pi * x) * F(x) + ph ** (-x) * mp.sin(mp.pi * x) ** 2 / mp.sqrt(5))) for x in xs)
    e3 = max(abs(F(x + 1) ** 2 - F(x) * F(x + 2) - mp.cos(mp.pi * x)) for x in xs)
    chk("recurrence / twist / Cassini at 40 random x in [-12, 12]: %.1e %.1e %.1e" % (e1, e2, e3), max(e1, e2, e3) < mp.mpf(10) ** -35)
    chk("F(0) = F(2) - F(1) = 0 is forced by the recurrence at x = 0 for every interpolant", fib(2) - fib(1) == 0)

    hdr("A3. the interpolant-family claim ('every function satisfying the recurrence with F(n) = F_n is F + r(x) phi^-x sin(pi x), r 1-periodic')")
    # general solution on R: P(x) phi^x + E(x) phi^-x with P 1-periodic and E 1-ANTIperiodic (E(x+1) = -E(x));
    # values at integers only pin P(n) = 1/sqrt5 and E(n) = -(-1)^n/sqrt5 AT INTEGERS.
    G0 = F0 + s ** 2 * a
    G1 = F1 + s ** 2 * phi * a
    G2 = F2 + s ** 2 * phi ** 2 * a
    chk("counterexample G = F + sin^2(pi x) phi^x satisfies the recurrence (sin^2 is 1-periodic)", sp.simplify(sp.expand(G2 - G1 - G0)) == 0)
    G = lambda x: F(x) + mp.sin(mp.pi * x) ** 2 * ph ** x
    chk("G(n) = F_n at all integers -60..60", max(abs(G(mp.mpf(k)) - fib(k)) for k in range(-60, 61)) < mp.mpf(10) ** -30)
    # is G - F = r(x) phi^-x sin(pi x) with r 1-periodic?  r(x) = sin(pi x) phi^(2x), r(x+1)/r(x) = -phi^2
    r = lambda x: (G(x) - F(x)) / (ph ** (-x) * mp.sin(mp.pi * x))
    x0 = mp.mpf('0.3')
    chk("G - F is NOT r phi^-x sin(pi x) with r 1-periodic: r(1.3)/r(0.3) = %s (would be 1)" % mp.nstr(r(x0 + 1) / r(x0), 8), abs(r(x0 + 1) / r(x0) + ph ** 2) < mp.mpf(10) ** -30)
    # a second counterexample INSIDE the claimed family (non-constant r) with smaller oscillation than Binet's:
    H0 = (a - c ** 3 / a) / r5
    H1 = (phi * a + c ** 3 / (phi * a)) / r5
    H2 = (phi ** 2 * a - c ** 3 / (phi ** 2 * a)) / r5
    chk("H = (phi^x - cos^3(pi x) phi^-x)/sqrt5 satisfies the recurrence (cos^3 is antiperiodic)", sp.simplify(sp.expand(H2 - H1 - H0)) == 0)
    H = lambda x: (ph ** x - mp.cos(mp.pi * x) ** 3 * ph ** (-x)) / mp.sqrt(5)
    chk("H(n) = F_n at all integers -60..60", max(abs(H(mp.mpf(k)) - fib(k)) for k in range(-60, 61)) < mp.mpf(10) ** -30)
    # H - F = cos(pi x) sin^2(pi x) phi^-x / sqrt5 = [sin(2 pi x)/(2 sqrt5)] phi^-x sin(pi x): r = sin(2 pi x)/(2 sqrt 5), 1-periodic
    mp.mp.dps = 20
    normF = mp.quad(lambda x: abs(F(x)) * ph ** x, [-30, -29, -28, -27, -26, -25])  # mean |F| phi^x sqrt... over 5 periods
    normH = mp.quad(lambda x: abs(H(x)) * ph ** x, [-30, -29, -28, -27, -26, -25])
    print("  mean oscillation on [-30,-25] (int |F| phi^x dx vs int |H| phi^x dx): %s vs %s  (limits 5*2/(pi sqrt5) = %s, 5*4/(3 pi sqrt5) = %s)"
          % (mp.nstr(normF, 8), mp.nstr(normH, 8), mp.nstr(10 / (mp.pi * mp.sqrt(5)), 8), mp.nstr(20 / (3 * mp.pi * mp.sqrt(5)), 8)))
    chk("H oscillates LESS than Binet's F on the negative axis (mean |.|): Binet is not the unique minimal-oscillation interpolant", normH < normF)
    print("  sup of |.| phi^x sqrt5 on [-30,-29]: F -> 1 and H -> 1 (both attained at the integer, where the values are forced): 'amplitude' does not single out Binet either")
    print("  verdict A3: the claimed family is INCOMPLETE (missing a(x) phi^x, a 1-periodic vanishing at integers) and 'unique minimal oscillation' is only shown for CONSTANT r.")
    print("  correct statement: interpolants = F + a(x) phi^x + r(x) phi^-x sin(pi x), a, r 1-periodic, a(n) = 0; with a = 0 and r constant the family is the note's,")
    print("  and both restrictions follow from an extra hypothesis such as 'entire of exponential type <= pi' (then periodic coefficients are constants), which the note does not state.")

    hdr("A4. zeros of F on the negative axis (mpmath findroot, 30 digits)")
    mp.mp.dps = 30
    zs = [mp.mpf(0)]
    zs.append(mp.findroot(F, mp.mpf('-0.18')))
    for k in range(1, 12):
        zs.append(mp.findroot(F, -(k + mp.mpf('0.5'))))
    print("  zeros:", [mp.nstr(z, 10) for z in zs])
    note = [0, -0.1838, -1.5708, -2.4704, -3.5109, -4.4958, -5.5016, -6.4994, -7.5002]
    chk("note's zero list (4 decimals) reproduced", all(abs(float(zs[i]) - note[i]) < 6e-5 for i in range(len(note))))
    # count sign changes on [-12, 0.5] with step 1e-3 (independent of findroot)
    g = [F(mp.mpf(i) / 1000) for i in range(-12000, 501)]
    sc = sum(1 for i in range(len(g) - 1) if g[i] * g[i + 1] < 0)
    print("  sign changes of F on [-12, 0.5] at step 1e-3: %d (plus the exact zero at 0): 13 zeros in [-12, 0] as the script found" % sc)
    chk("zeros converge to the half-integers -(k+1/2): |z_k + k + 1/2| for k = 6..11 = %s" % [mp.nstr(abs(zs[k + 1] + k + mp.mpf('0.5')), 3) for k in range(6, 12)],
        all(abs(zs[k + 1] + k + mp.mpf('0.5')) < 0.01 for k in range(6, 12)))
    print("  (the zero -1.570776... is NOT -pi/2 = -1.570796...; a coincidence to 4 decimals, not claimed by the note)")
    # amplitude of F_r for constant r
    print("  F_r = F + r phi^-x sin(pi x)/sqrt5: on the negative axis F_r ~ -phi^|x| (cos(pi x) - r sin(pi x))/sqrt5, envelope phi^|x| sqrt(1+r^2)/sqrt5, minimal at r = 0 among CONSTANT r: correct as far as it goes.")


# ----------------------------------------------------------------------------------------------------------
# B. shadow proposition and the cycle computations
# ----------------------------------------------------------------------------------------------------------
def U(n):
    m = 3 * n + 1
    v = (m & -m).bit_length() - 1
    return m >> v, v


def Ufrac(x):
    m = 3 * x + 1
    v = 0
    while m.numerator % 2 == 0:
        m /= 2
        v += 1
    return m, v


def word_of(n, p):
    w = []
    for _ in range(p):
        n, v = U(n)
        w.append(v)
    return tuple(w), n


def affine(w):
    p = len(w)
    A = sum(w)
    S = 0
    At = 0
    for t in range(p):
        S += 3 ** (p - 1 - t) * 2 ** At
        At += w[t]
    return p, A, S


def res2(x, K):
    """residue mod 2^K of a rational with odd denominator"""
    return (x.numerator * pow(x.denominator, -1, 2 ** K)) % 2 ** K


def nodescent(w):
    At = 0
    for j, v in enumerate(w, 1):
        At += v
        if 2 ** At > 3 ** j:
            return False
    return True


def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def primitive(w):
    return all(w != w[i:] + w[:i] for i in range(1, len(w)))


def compositions(A, p):
    if p == 1:
        yield (A,)
        return
    for v in range(1, A - p + 2):
        for rest in compositions(A - v, p - 1):
            yield (v,) + rest


def partB():
    hdr("B1. shadow proposition, first statement: odd n with EXACT valuation word w form ONE class mod 2^(A+1) (brute force)")
    rng = random.Random(3)
    words = [(1,), (1, 2), (1, 1, 2), (1, 1, 1, 2, 1, 1, 4), (2, 1), (1, 3), (2, 2, 1), (3, 1, 1, 1, 1), (1, 1, 1, 1, 4), (4,), (1, 1, 5, 1)]
    while len(words) < 17:
        p = rng.randint(2, 6)
        w = tuple(rng.randint(1, 4) for _ in range(p))
        if w not in words:
            words.append(w)
    allok = True
    for w in words:
        p, A, S = affine(w)
        N = 2 ** (A + 3)
        exact = [n for n in range(-N + 1, N, 2) if word_of(n, p)[0] == w]
        # coarse: first p-1 exact and last >= v_p
        def coarse(n):
            ww, _ = word_of(n, p)
            return ww[:-1] == w[:-1] and ww[-1] >= w[-1]
        crs = [n for n in range(-N + 1, N, 2) if coarse(n)]
        M = 2 ** (A + 1)
        r = exact[0] % M
        xw = Fraction(S, 2 ** A - 3 ** p)
        ok = (len(exact) == 2 * N // M and all(n % M == r for n in exact)
              and r == res2(xw, A + 1) and r == ((2 ** A - S) * pow(3, -p, M)) % M
              and len(crs) == 2 * N // 2 ** A and all(n % 2 ** A == crs[0] % 2 ** A for n in crs)
              and crs[0] % 2 ** A == (-S * pow(3, -p, 2 ** A)) % 2 ** A)
        # x_w itself has word w^3 (as a 2-adic number: parity of a rational with odd denominator is the numerator's)
        y = xw
        ww = []
        for _ in range(3 * p):
            y, v = Ufrac(y)
            ww.append(v)
        ok &= tuple(ww) == w * 3 and y == xw
        allok &= ok
        print("  w = %-22s A = %2d S = %6d x_w = %10s  exact class = %d mod 2^%d (%d members in [-2^%d, 2^%d)), coarse class mod 2^%d has %d: %s"
              % (w, A, S, xw, r, A + 1, len(exact), A + 3, A + 3, A, len(crs), "ok" if ok else "MISMATCH"))
    chk("exact word <=> one class mod 2^(A+1) = x_w mod 2^(A+1) = THM-4512's E_w; coarse (last valuation >= v_p) <=> one class mod 2^A = C_w; x_w has word w^inf", allok)

    hdr("B2. U^(mp)(n) = (3^p/2^A)^m (n - x_w) + x_w for n = x_w mod 2^(mA+1), both signs, exact rationals")
    allok = True
    for w in words:
        p, A, S = affine(w)
        xw = Fraction(S, 2 ** A - 3 ** p)
        for m in (1, 2, 3):
            K = m * A + 1
            r = res2(xw, K)
            for _ in range(12):
                n = r + 2 ** K * rng.randint(-2000, 2000)
                if n == 0:
                    continue
                ww, end = word_of(n, m * p)
                allok &= ww == w * m and Fraction(end) == Fraction(3 ** p, 2 ** A) ** m * (n - xw) + xw
                # and one MORE period is generally not followed when v2(n - x) < (m+1)A + 1
    chk("shadow formula on 17 words x m = 1,2,3 x 12 random members of each sign", allok)
    # 'exactly floor((v2(n-x)-1)/A) periods'
    ok = True
    for w in words[:6]:
        p, A, S = affine(w)
        xw = Fraction(S, 2 ** A - 3 ** p)
        for _ in range(200):
            m = rng.randint(0, 4)
            K = m * A + 1 + rng.randint(0, A - 1)
            n = res2(xw, K) + 2 ** K * (2 * rng.randint(-500, 500) + 1)   # n = x_w mod 2^K; v2(n - x_w) >= K, computed exactly below
            if n == xw:
                continue                                                   # the cycle point itself has word w^inf (v2 = infinity)
            v2 = 0
            d = n - xw
            while d.numerator % 2 == 0:
                d /= 2
                v2 += 1
            assert v2 >= K
            periods = 0
            while periods <= (v2 - 1) // A + 1 and word_of(n, (periods + 1) * p)[0] == w * (periods + 1):
                periods += 1
            ok &= periods == (v2 - 1) // A
    chk("number of exactly-shadowed periods = floor((v2(n - x_w) - 1)/A) (200 random n per word, 6 words)", ok)
    chk("sign: x_w < 0 iff 3^p > 2^A (S_w > 0 always, its t = 0 term is 3^(p-1))", all((Fraction(affine(w)[2], 2 ** affine(w)[1] - 3 ** len(w)) < 0) == (3 ** len(w) > 2 ** affine(w)[1]) for w in words))

    hdr("B3. negative odd cycles of U with |n| <= 10^6: full basin scan")
    limit = 10 ** 6
    cyc_of = {}
    cycles = []
    maxabs = 0
    for n in range(-1, -limit - 1, -2):
        if n in cyc_of:
            continue
        path = []
        pset = set()
        m = n
        while m not in cyc_of and m not in pset:
            pset.add(m)
            path.append(m)
            m, _ = U(m)
            maxabs = max(maxabs, abs(m))
            if abs(m) > 10 ** 15:
                raise RuntimeError("orbit escaped: %d" % n)
        if m in cyc_of:
            cid = cyc_of[m]
        else:
            cid = len(cycles)
            cycles.append(path[path.index(m):])
        for q in path:
            cyc_of[q] = cid
    print("  every odd n in [-10^6, -1] reaches one of %d cycles; largest |value| met on the way: %d" % (len(cycles), maxabs))
    for cyc in cycles:
        x = min(cyc, key=abs)
        i = cyc.index(x)
        cyc = cyc[i:] + cyc[:i]
        w, back = word_of(x, len(cyc))
        p, A, S = affine(w)
        print("  cycle %s: word %s, p = %d, A = %d, S = %d, 3^p - 2^A = %d, x_w = %s, factor 3^p/2^A = %s = %.5f"
              % (cyc, w, p, A, S, 3 ** p - 2 ** A, Fraction(S, 2 ** A - 3 ** p), Fraction(3 ** p, 2 ** A), 3 ** p / 2 ** A))
    mins = sorted(min(c, key=abs) for c in cycles)
    chk("exactly three negative odd cycles meet |n| <= 10^6 (and every such n lies in their basin): %s" % mins, mins == [-17, -5, -1])
    chk("3^7 - 2^11 = 139 (the owner's 139) and x_w = 2363/(-139) = -17", 3 ** 7 - 2 ** 11 == 139 and Fraction(2363, -139) == -17)
    # the script's stated residues
    chk("script residues: -1 mod 4/16/2^11 = 3, 15, 2047; -5 mod 16/2^10/2^31 = 11, 1019, 2147483643; -17 mod 2^12 = 4079",
        [(-1) % 4, (-1) % 16, (-1) % 2 ** 11, (-5) % 16, (-5) % 2 ** 10, (-5) % 2 ** 31, (-17) % 2 ** 12] == [3, 15, 2047, 11, 1019, 2147483643, 4079])

    hdr("B4. the '120 primitive no-descent words with p <= 8 up to rotation'")
    # (i) the script's rule: enumerate no-descent words (2^A_j < 3^j at every prefix), keep those that are (a) their own
    #     minimal rotation and (b) not a proper power.
    script_count = 0
    script_set = set()
    for p in range(1, 9):
        def rec(prefix, At):
            nonlocal script_count
            j = len(prefix)
            if j == p:
                w = tuple(prefix)
                if w == canon(w) and primitive(w):
                    script_count += 1
                    script_set.add(w)
                return
            v = 1
            while 2 ** (At + v) < 3 ** (j + 1):
                rec(prefix + [v], At + v)
                v += 1
        rec([], 0)
    chk("script's rule reproduced: %d" % script_count, script_count == 120)
    # (ii) ALL primitive growth necklaces (3^p > 2^A), i.e. all primitive negative rational cycles of period p <= 8
    neck = {}
    nd_words = 0
    for p in range(1, 9):
        A = p
        while 2 ** A < 3 ** p:
            for w in compositions(A, p):
                if nodescent(w):
                    nd_words += 1
                if primitive(w):
                    cw = canon(w)
                    if cw not in neck:
                        neck[cw] = []
                    if nodescent(w):
                        neck[cw].append(w)
            A += 1
    # (iii) Moebius: number of primitive necklaces of length p and sum A = (1/p) sum_{d | gcd(p,A)} mu(d) C(A/d - 1, p/d - 1)
    def mu(n):
        return int(sp.mobius(n))
    moeb = 0
    for p in range(1, 9):
        A = p
        while 2 ** A < 3 ** p:
            g = math.gcd(p, A)
            moeb += sum(mu(d) * math.comb(A // d - 1, p // d - 1) for d in range(1, g + 1) if g % d == 0) // p
            A += 1
    print("  no-descent words of length <= 8 (all rotations counted separately): %d" % nd_words)
    print("  primitive growth necklaces (rational negative cycles) with p <= 8: brute force %d, Moebius %d" % (len(neck), moeb))
    print("  necklaces with at least one no-descent rotation (cycle lemma predicts all): %d" % sum(1 for v in neck.values() if v))
    print("  necklaces whose LEX-MIN rotation is no-descent (= the script's count): %d" % sum(1 for cw in neck if nodescent(cw)))
    missed = sorted(cw for cw in neck if not nodescent(cw))
    print("  necklaces MISSED by the script's rule (min rotation has a coefficient descent although another rotation is no-descent):")
    for cw in missed:
        p, A, S = affine(cw)
        print("    min rotation %s (A = %d, 3^p/2^A = %.4f); no-descent rotation(s) %s; x_w = %s" % (cw, A, 3 ** p / 2 ** A, neck[cw], Fraction(S, 2 ** A - 3 ** p)))
    chk("'120' equals the number of primitive growth necklaces with p <= 8", len(neck) == 120,
        "it does not: %d necklaces; 120 = those whose lexicographically minimal rotation happens to be no-descent" % len(neck))
    chk("cycle lemma: every growth necklace has a no-descent rotation", all(v for v in neck.values()))
    chk("Moebius count agrees with brute force", moeb == len(neck))
    ints = sorted((cw, Fraction(affine(cw)[2], 2 ** affine(cw)[1] - 3 ** len(cw))) for cw in neck if Fraction(affine(cw)[2], 2 ** affine(cw)[1] - 3 ** len(cw)).denominator == 1)
    print("  integer fixed points among ALL %d primitive growth necklaces p <= 8: %s" % (len(neck), ints))
    chk("exactly (1) -> -1, (1,2) -> -5, (1,1,1,2,1,1,4) -> -17 (rotation-invariant property, so the missed necklaces cannot add any)",
        sorted(x for _, x in ints) == [-17, -5, -1] and set(cw for cw, _ in ints) == {(1,), (1, 2), (1, 1, 1, 2, 1, 1, 4)})
    # extend the integer-fixed-point census to p <= 12 over all primitive growth necklaces
    ints12 = []
    count12 = {}
    for p in range(1, 13):
        A = p
        cnt = 0
        seen = set()
        while 2 ** A < 3 ** p:
            for w in compositions(A, p):
                cw = canon(w)
                if cw in seen or not primitive(w):
                    continue
                seen.add(cw)
                cnt += 1
                pp, AA, S = affine(cw)
                if S % (3 ** p - 2 ** A) == 0:
                    ints12.append((cw, Fraction(S, 2 ** A - 3 ** p)))
            A += 1
        count12[p] = cnt
    print("  primitive growth necklaces by period p = 1..12: %s (total %d)" % (count12, sum(count12.values())))
    chk("integer negative cycle points with primitive period <= 12: only -1, -5, -17", sorted(x for _, x in ints12) == [-17, -5, -1])
    w112 = (1, 1, 2)
    p, A, S = affine(w112)
    chk("(1,1,2): x = -19/11, factor 27/16 (note's example)", Fraction(S, 2 ** A - 3 ** p) == Fraction(-19, 11) and Fraction(27, 16) == Fraction(3 ** p, 2 ** A))

    hdr("B5. the two families")
    fam = [4 * 8 ** k - 5 for k in range(1, 8)]
    chk("27 family 4*8^k - 5 = %s, n_(k+1) = 8 n_k + 35" % fam, fam[:4] == [27, 251, 2043, 16379] and all(fam[i + 1] == 8 * fam[i] + 35 for i in range(6)))
    chk("v2(n_k + 5) = 3k + 2 exactly (so n_k = -5 mod 2^(3k+2), stronger than the note's mod 2^(3k)); shadowed periods floor((3k+1)/3) = k",
        all((fam[k - 1] + 5) == 2 ** (3 * k + 2) for k in range(1, 8)))
    ok = True
    for k in range(1, 8):
        n = fam[k - 1]
        for j in range(0, k + 1):
            ww, end = word_of(n, 2 * j)
            ok &= end == 4 * 9 ** j * 8 ** (k - j) - 5 and ww == (1, 2) * j
    chk("lane's formula U^(2j)(n_k) = 4 9^j 8^(k-j) - 5 for 0 <= j <= k, word (1,2)^j = shadow formula with x = -5", ok)
    fam2 = [(4 ** k + 17) // 3 for k in range(1, 12)]
    chk("(4^k+17)/3 = %s, n_(k+1) = 4 n_k - 17, integers" % fam2[:8], fam2[:8] == [7, 11, 27, 91, 347, 1371, 5467, 21851]
        and all((4 ** k + 17) % 3 == 0 for k in range(1, 12)) and all(fam2[i + 1] == 4 * fam2[i] - 17 for i in range(10)))
    x = Fraction(17, 3)
    chk("v2(n_k - 17/3) = 2k (n_k - 17/3 = 4^k/3): 2-adic limit 17/3", all(Fraction(fam2[k - 1]) - x == Fraction(4 ** k, 3) for k in range(1, 12)))
    orb = [x]
    for _ in range(7):
        y, v = Ufrac(orb[-1])
        orb.append(y)
    chk("orbit of 17/3 under U: %s" % [str(o) for o in orb], [str(o) for o in orb] == ['17/3', '9', '7', '11', '17', '13', '5', '1'])
    # first-descent times for k = 3..30
    fd = []
    for k in range(3, 31):
        n = fam2[k - 1] if k <= 11 else (4 ** k + 17) // 3
        m = n
        j = 0
        while True:
            m, _ = U(m)
            j += 1
            if m < n:
                break
        fd.append(j)
    chk("first-descent times k = 3..9: %s; k = 10..30 all 6: %s" % (fd[:7], all(t == 6 for t in fd[7:])), fd[:7] == [37, 28, 6, 6, 6, 6, 6] and all(t == 6 for t in fd[7:]))
    # the exact iterates: U^j(n_k) = (3^j/2^(A_j)) (n_k - 17/3) + U^j(17/3) as long as 2k >= A_j + 1 (17/3 is NOT a cycle point:
    # the cycle-shadow formula does not apply; this is the general affine form), reproducing nextforest's 2^(t-1)+9, 3 2^(t-3)+7, ...
    wd = []
    y = x
    for _ in range(6):
        y, v = Ufrac(y)
        wd.append(v)
    ok = True
    for k in range(6, 20):
        n = (4 ** k + 17) // 3
        At = 0
        for j in range(1, 7):
            At += wd[j - 1]
            ww, end = word_of(n, j)
            ok &= Fraction(end) == Fraction(3 ** j, 2 ** At) * (n - x) + orb[j] and ww == tuple(wd[:j])
    chk("word of 17/3 = %s (A_6 = %d, 3^6 < 2^10: a DESCENT word); U^j(n_k) = (3^j/2^A_j)(n_k - 17/3) + U^j(17/3) for k >= 6, j <= 6" % (tuple(wd), sum(wd)), ok)
    print("  hence U^6(n_k) = 243 4^(k-5) + 5 < n_k = (4^k+17)/3 for k >= 6 and the first five iterates 2^(2k-1)+9, 3 2^(2k-3)+7, 9 2^(2k-4)+11, 27 2^(2k-5)+17, 81 2^(2k-7)+13 exceed n_k:")
    print("  exactly the constants 9, 7, 11, 17, 13, 5 of nextforest_20260926_boundary.md section 6 -- the note's '17/3' reading is a correct explanation of that lane's formulas (k = 5 needs the separate check: 2k = 10 < A_6 + 1 = 11).")

    hdr("B6. D(k) of S11 (density among odd integers of exact no-descent words of length k) and the periodic fraction")
    cnt = [dict() for _ in range(61)]
    cnt[0][0] = 1
    for j in range(1, 61):
        for A0, c0 in cnt[j - 1].items():
            v = 1
            while 2 ** (A0 + v) < 3 ** j:
                cnt[j][A0 + v] = cnt[j].get(A0 + v, 0) + c0
                v += 1
    D = lambda k: sum(Fraction(c0, 2 ** A0) for A0, c0 in cnt[k].items())
    s11 = {1: Fraction(1, 2), 2: Fraction(3, 8), 3: Fraction(1, 4), 4: Fraction(13, 64), 5: Fraction(19, 128), 6: Fraction(1, 8), 7: Fraction(113, 1024), 8: Fraction(367, 4096)}
    chk("D(1..8) = %s (matches the S11 table)" % [str(D(k)) for k in range(1, 9)], all(D(k) == s11[k] for k in range(1, 9)))
    s11f = {12: .05212, 16: .03092, 20: .01925, 24: .01314, 28: .00875, 32: .00593, 36: .00427, 40: .00299, 41: .00265, 48: .00156, 52: .00112, 56: .00082, 60: .00062}
    chk("D(12..60) match S11 to 5 decimals: %s" % {k: "%.5f" % float(D(k)) for k in s11f}, all(abs(float(D(k)) - s11f[k]) < 6e-6 for k in s11f))
    chk("D(41) = 1530343662856563/2^59 exactly", D(41) == Fraction(1530343662856563, 2 ** 59))
    print("  periodic no-descent words (w = u^m, m >= 2; u no-descent <=> u^m no-descent) as a fraction of D(k), upper bound summing over all divisors:")
    for k in (4, 8, 12, 16, 24, 32, 40, 48, 60):
        per = sum(sum(Fraction(c0, 2 ** ((k // d) * A0)) for A0, c0 in cnt[d].items()) for d in range(1, k) if k % d == 0)
        print("    k = %2d: D(k) = %.3e, periodic part <= %.3e, ratio <= %.2e" % (k, float(D(k)), float(per), float(per / D(k))))
    print("  the 'aperiodic words dominate by an exponential margin' claim is supported (ratio ~ 2^(-0.4 k)); the note's own evidence (3 2^-K vs D(k)) only concerns the three integer rules.")
    print("  note: the script's 'residual at ~K/1.58 Syracuse steps is of order 10^-2..10^-3' is loose: D(8) = 0.0896 at K = 12 bits.")


def oddpart(x):
    return x >> ((x & -x).bit_length() - 1)


def partB7():
    hdr("B7. what the owner's summary refers to: the codex lane's three rule NODES (entry_20260927_board.md (B1)-(B3)) and its decreasing-parameter family N(k,h) ((B5)-(B6))")
    rng = random.Random(11)
    ok = True
    for k in range(0, 7):
        t = 1
        while 2 ** t * (8 ** (k + 1) - 5) <= 9 ** (k + 1) - 5:
            t += 1
        b0 = (5 * pow(9 ** (k + 1), -1, 2 ** t)) % 2 ** t           # odd, since 5 and the inverse are odd
        for b in (b0, b0 + 2 ** t * 2 * rng.randint(1, 50)):
            n = b * 8 ** (k + 1) - 5
            m = n
            j = 0
            word = []
            while True:
                m, v = U(m)
                j += 1
                word.append(v)
                if m < n:
                    break
            ok &= j == 2 * k + 2 and tuple(word[:2 * k]) == (1, 2) * k and word[2 * k] == 1 and m == oddpart(b * 9 ** (k + 1) - 5)
            print("  k = %d, t = %2d, b = %6d: n = %d, first descent at step %d, word %s" % (k, t, b, n, j, tuple(word)))
    chk("(B2) cylinders: exact first-descent time 2k+2 with word (1,2)^k . 1 . t' = the nodes Repeat12(k), Step, Step: a -5 shadow for k periods, one v = 1 step, one descent step; the -17 word never occurs", ok)
    print("  so the lane's 'three parameterized rules' are three NODES of one certificate grammar for the cylinders b 8^(k+1) - 5, not three growth rules matching the three cycles -1, -5, -17.")
    print("  the owner's 'separate closed family through 27 with a decreasing recursion parameter' is (B5)-(B6): N(k,h) = 2(47 4^h + 7) 8^k / 3^(2k+1) - 5 with 47 4^h = -7 mod 3^(2k+1):")
    ok = True
    fam2set = set((4 ** k + 17) // 3 for k in range(1, 400))
    for k in range(1, 6):
        mod = 3 ** (2 * k + 1)
        h = 0
        pw = 1
        while (47 * pw + 7) % mod != 0:
            h += 1
            pw = pw * 4 % mod
        num = 2 * (47 * 4 ** h + 7) * 8 ** k
        assert num % mod == 0
        N = num // mod - 5
        v2 = ((N + 5) & -(N + 5)).bit_length() - 1
        Nk1 = 2 * (47 * 4 ** h + 7) * 8 ** (k - 1) // 3 ** (2 * k - 1) - 5
        m2 = U(U(N)[0])[0]
        N0 = 2 * (47 * 4 ** h + 7) // 3 - 5
        u0 = U(N0)[0]
        m = N
        j = 0
        fd = None
        first1 = None
        while first1 is None:
            m, v = U(m)
            j += 1
            if fd is None and m < N:
                fd = j
            if m == 1:
                first1 = j
        # v2(N+5) = 3k + v2(a), a = 2(47 4^h + 7)/3^(2k+1): = 3k+1 for h >= 1 (47 4^h + 7 odd) and 3k+2 for h = 0 (54 = 2 * 27)
        ok &= m2 == Nk1 and v2 == 3 * k + (2 if h == 0 else 1) and u0 == 47 and (fd == 2 * k + 1 or N == 27) and first1 == 2 * k + 39 and ((N in fam2set) == (N == 27))
        print("  k = %d, least h = %5d: N = %s (%d bits), v2(N+5) = %d (3k+1, or 3k+2 when h = 0), U^2 N = N(k-1,h): %s, U N(0,h) = 47: %s, first descent at step %s, first 1 at step %d, in the (4^k+17)/3 family: %s"
              % (k, h, N if N < 10 ** 15 else "(large)", N.bit_length(), v2, m2 == Nk1, u0 == 47, fd, first1, N in fam2set))
    chk("N(k,h): a -5 GROWTH shadow for k periods (N = -5 mod 2^(3k+1), 2k rising steps) exiting to 47 by one large-valuation step (first descent 2k+1 except for 27), first 1 at 2k+39; meets the (4^k+17)/3 family only at 27", ok)
    print("  hence the note's item 3 ('the \"decreasing recursion parameter\" family is a descent shadow, not a growth one') attaches the owner's phrase to the wrong family;")
    print("  the decreasing parameter k counts the remaining -5 shadow periods, and the family is a growth shadow with a designed exit, the opposite of what the note says.")


# ----------------------------------------------------------------------------------------------------------
# C. negaFibonacci
# ----------------------------------------------------------------------------------------------------------
FNEG = [0] + [fib(-k) for k in range(1, 61)]   # FNEG[k] = F_(-k)


def nf_range(m):
    """integers representable with indices <= m: m odd: [-(F_m - 1), F_(m+1)], m even: [-(F_(m+1) - 1), F_m]"""
    if m % 2:
        return -(fib(m) - 1), fib(m + 1)
    return -(fib(m + 1) - 1), fib(m)


ODD_HI = [(fib(m + 1), m) for m in range(1, 60, 2)]          # smallest odd m with F_(m+1) >= N  (N > 0)
EVEN_LO = [(fib(m + 1) - 1, m) for m in range(2, 61, 2)]     # smallest even m with F_(m+1) - 1 >= -N (N < 0)
ODD_KEYS = [t[0] for t in ODD_HI]
EVEN_KEYS = [t[0] for t in EVEN_LO]


def nega(N):
    """Knuth's negaFibonacci representation as a tuple of indices k (F_(-k)), no two consecutive"""
    idx = []
    while N != 0:
        if N > 0:
            m = ODD_HI[bisect.bisect_left(ODD_KEYS, N)][1]
        else:
            m = EVEN_LO[bisect.bisect_left(EVEN_KEYS, -N)][1]
        idx.append(m)
        N -= FNEG[m]
    return tuple(reversed(idx))


def nega_lowest(N):
    low = 0
    while N != 0:
        if N > 0:
            m = ODD_HI[bisect.bisect_left(ODD_KEYS, N)][1]
        else:
            m = EVEN_LO[bisect.bisect_left(EVEN_KEYS, -N)][1]
        low = m
        N -= FNEG[m]
    return low


def zeck(N):
    idx = []
    k = 2
    while fib(k + 1) <= N:
        k += 1
    while N > 0:
        while fib(k) > N:
            k -= 1
        idx.append(k)
        N -= fib(k)
        k -= 2
    return tuple(sorted(idx))


FIBS = [fib(i) for i in range(0, 95)]


def zeck_lowest(N):
    k = bisect.bisect_right(FIBS, N) - 1
    low = 0
    while N > 0:
        while FIBS[k] > N:
            k -= 1
        low = k
        N -= FIBS[k]
        k -= 2
    return low


def floor_phi(M):
    """floor(M phi) exactly for integer M (phi = (1+sqrt5)/2)"""
    if M >= 0:
        return (M + isqrt(5 * M * M)) // 2
    return -floor_phi(-M) - 1          # M phi irrational for M != 0


def upperW(M):
    """M >= 1 is in the upper Wythoff sequence {floor(N phi^2)} iff floor((M+1)phi) - floor(M phi) = 1 (iff {M phi} < phi^-2)"""
    return floor_phi(M + 1) - floor_phi(M) == 1


def partC():
    hdr("C1. negaFibonacci: existence and uniqueness by enumeration (all non-consecutive subsets of {1..30})")
    K = 30
    sums = {}
    def rec(i, cur, mx):
        if i > K:
            sums[cur] = sums.get(cur, 0) + 1
            return
        rec(i + 1, cur, mx)
        rec(i + 2, cur + FNEG[i], i)
    rec(1, 0, 0)
    tot = sum(sums.values())
    chk("number of subsets = F_32 = %d" % tot, tot == fib(K + 2))
    lo, hi = nf_range(K)
    chk("sums are exactly the integers of [%d, %d], each once (uniqueness on this range)" % (lo, hi),
        set(sums) == set(range(lo, hi + 1)) and all(v == 1 for v in sums.values()))
    # the range claim level by level (needed to conclude that larger indices cannot reach |N| <= 10^5)
    ok = True
    for m in range(1, 25):
        S = set()
        def rec2(i, cur):
            if i > m:
                S.add(cur)
                return
            rec2(i + 1, cur)
            rec2(i + 2, cur + FNEG[i])
        rec2(1, 0)
        l, h = nf_range(m)
        ok &= S == set(range(l, h + 1))
    chk("range(m) = [-(F_m - 1), F_(m+1)] (m odd) / [-(F_(m+1) - 1), F_m] (m even) for m <= 24 (and m = 30 above)", ok)
    print("  a subset with largest index k >= 31 has |sum| >= F_k - F_(k-2) = F_(k-1) >= F_30 = 832040 > 10^5 (the opposite-sign terms have indices of the other parity <= k-2 and sum to at most F_(k-2)),")
    print("  so the enumeration with indices <= 30 is complete for |N| <= 10^5: existence and uniqueness there are PROVED by this run.")
    X = 10 ** 5
    ok = all(nega(N) in ((),) if N == 0 else True for N in (0,))
    # own greedy agrees with the enumeration for all |N| <= 10^5 and is non-consecutive
    def rec3(i, cur, chosen, out):
        if i > 28:
            if abs(cur) <= X:
                out[cur] = tuple(sorted(chosen))
            return
        rec3(i + 1, cur, chosen, out)
        rec3(i + 2, cur + FNEG[i], chosen + [i], out)
    enum = {}
    rec3(1, 0, [], enum)
    ok = all(nega(N) == enum[N] and all(b - a >= 2 for a, b in zip(nega(N), nega(N)[1:])) for N in range(-X, X + 1))
    chk("interval-greedy representation = enumerated subset for all |N| <= 10^5 (non-consecutive, unique)", ok)
    small = {1: (1,), -1: (2,), 2: (3,), -2: (1, 4), 3: (1, 3), -3: (4,), 4: (2, 5), -4: (2, 4), 5: (5,), 8: (1, 3, 5)}
    chk("note's small cases 1, -1, 2, -2, 3, -3, 4, -4, 5, 8", all(nega(N) == v for N, v in small.items()))
    print("  -8..8:", ["%d=%s" % (N, "+".join("F(-%d)" % k for k in nega(N)) or "0") for N in range(-8, 9)])

    hdr("C2. class densities by lowest index (k = 1: +1, k = 2: -1, odd >= 3, even >= 4), both signs, |N| <= 10^5 and 10^6")
    for X in (10 ** 5, 10 ** 6):
        for sign, rngN in (("positive", range(1, X + 1)), ("negative", range(-X, 0))):
            cnt = {}
            perk = {}
            for N in rngN:
                k = nega_lowest(N)
                perk[k] = perk.get(k, 0) + 1
                c = 'k=1' if k == 1 else 'k=2' if k == 2 else 'odd>=3' if k % 2 else 'even>=4'
                cnt[c] = cnt.get(c, 0) + 1
            dens = {c: cnt[c] / X for c in ('k=1', 'k=2', 'odd>=3', 'even>=4')}
            target = {'k=1': PHI ** -2, 'k=2': PHI ** -3, 'odd>=3': PHI ** -3, 'even>=4': PHI ** -4}
            dev = max(abs(dens[c] - target[c]) for c in dens)
            print("  X = %d %s: %s  max |dens - phi^-(2,3,3,4)| = %.2e; per index k=1..8: %s"
                  % (X, sign, {c: "%.5f" % dens[c] for c in dens}, dev, {k: "%.5f" % (perk.get(k, 0) / X) for k in range(1, 9)}))
            if X == 10 ** 5:
                chk("%s |N| <= 10^5: 4-decimal densities 0.3820/0.2361/0.2361/0.1459 (note)" % sign,
                    ["%.4f" % dens[c] for c in ('k=1', 'k=2', 'odd>=3', 'even>=4')] in (['0.3820', '0.2361', '0.2361', '0.1459'], ['0.3820', '0.2361', '0.2360', '0.1459']))
            chk("%s |N| <= %d: within %s of phi^-2, phi^-3, phi^-3, phi^-4 (per-index density ~ phi^-(k+1))" % (sign, X, "1e-4" if X == 10 ** 5 else "3e-5"), dev < (1e-4 if X == 10 ** 5 else 3e-5))

    hdr("C3. Zeckendorf classes on positives and the two unit classes")
    X = 10 ** 5
    zc = {}
    for N in range(1, X + 1):
        k = zeck_lowest(N)
        c = 'k=2' if k == 2 else 'odd>=3' if k % 2 else 'even>=4'
        zc[c] = zc.get(c, 0) + 1
    print("  Zeckendorf lowest-index classes N <= 10^5: %s (phi^-2 = %.4f, phi^-3 = %.4f)" % ({c: "%.4f" % (v / X) for c, v in sorted(zc.items())}, PHI ** -2, PHI ** -3))
    chk("Zeckendorf classes 0.3820, 0.3820, 0.2361 (S8 Wythoff classes)", ["%.4f" % (zc[c] / X) for c in ('k=2', 'odd>=3', 'even>=4')] == ['0.3820', '0.3820', '0.2361'])
    for X in (10 ** 5, 10 ** 6):
        z1 = n1 = both = 0
        okA = okB = True
        for N in range(1, X + 1):
            a_ = zeck_lowest(N) == 2
            b_ = nega_lowest(N) == 1
            z1 += a_
            n1 += b_
            both += a_ and b_
            okB &= a_ == upperW(N + 1)
            okA &= b_ == (N == 1 or upperW(N - 1))
        print("  X = %d: Zeckendorf ends in F_2: %d (%.5f), negaFibonacci ends in F_(-1): %d (%.5f), both: %d (%.6f); phi^-4 X = %.2f, deviation %+.2f"
              % (X, z1, z1 / X, n1, n1 / X, both, both / X, PHI ** -4 * X, both - PHI ** -4 * X))
        chk("X = %d: Zeckendorf ends in F_2  <=>  N+1 is in the upper Wythoff sequence  <=>  {(N+1) phi} < phi^-2   (exact, isqrt arithmetic)" % X, okB)
        chk("X = %d: negaFibonacci ends in F_(-1)  <=>  N = 1 or N-1 is in the upper Wythoff sequence  <=>  {(N-1) phi} < phi^-2" % X, okA)
        chk("X = %d: intersection count is within 2 of phi^-4 X (bounded discrepancy => the density is EXACTLY phi^-4)" % X, abs(both - PHI ** -4 * X) < 2)
    print("  proof of the independence: in u = {N phi}, class B (Zeckendorf ends in F_2) is u in [phi^-2, 2 phi^-2) and class A (nega ends in F_(-1)) is u in [phi^-1, 1);")
    print("  the intersection [phi^-1, 2 phi^-2) has length 2 phi^-2 - phi^-1 = phi^-2 (2 - phi) = phi^-4 exactly (2 - phi = phi^-2). Equidistribution of N phi mod 1 gives the densities.")
    # R_3 = floor(phi^2 Z): integers whose negaFibonacci representation uses only indices >= 3
    X = 10 ** 5
    R3 = set(N for N in range(-X, X + 1) if N == 0 or nega(N)[0] >= 3)
    W2 = set()
    M = 0
    while True:
        v = floor_phi(M) + M           # floor(M phi^2) = floor(M phi) + M
        if v > X:
            break
        W2.add(v)
        M += 1
    M = -1
    while True:
        v = floor_phi(M) + M
        if v < -X:
            break
        W2.add(v)
        M -= 1
    chk("R_3 (indices >= 3 only) = {floor(N phi^2): N in Z} on [-10^5, 10^5] (shift by two indices = multiplication by phi^2 minus the golden expansion sum {N phi^2})", R3 == W2)
    same = [N for N in range(1, 10 ** 5 + 1) if tuple(k - 1 for k in zeck(N)) == nega(N)]
    chk("'negaFibonacci index set = Zeckendorf set shifted by one' holds for exactly N = 1 among N <= 10^5 (note says 'for no N'): %s" % same, same == [1])


def main():
    partA()
    partB()
    partB7()
    partC()
    hdr("summary")
    print(" failed checks: %d %s" % (len(FAILS), FAILS))
    print(" elapsed %.1f s" % (time.time() - T0))


if __name__ == '__main__':
    main()
