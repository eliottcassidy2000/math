#!/usr/bin/env python3
"""
procgen_rank_20260926_q1.py -- Q1 of the rank lane: infinite banks of 2-adic valuation counters.

Parts (numbering as in the note, section 2):
  A  the product-formula (Liouville) bound v2(n - beta) <= log2((n+1) H(beta)); heights of periodic points
  B  Lemma L (local expansion of a bank potential at a rational point, with the explicit error bound)
  C  Theorem A witnesses: every summable bank is violated along the chain from a preimage z of -1 whose
     charge is too small; the violation grows linearly in the depth with the predicted slope
  D  sufficiency inside shadows: charge c > a*chi on the cycle points plus a periodic h descends along the
     whole shadow of -1, of the -5 cycle and of the -17 cycle; the entry step from an uncharged preimage fails
  E  the forced mass per period, M(p) = sum over expanding primitive necklaces of log(3^a/2^p)
  F  Theorem D: the private-atom bank built from stopping times descends on a finite window
Every printed claim is a check(...).  a = 1 throughout (the rank a log n + Phi scales).
"""
import math
import random
from fractions import Fraction

from procgen_rank_20260926_lib import (KAPPA, LN2, LN3, H_BITS, check, say, v2, height, T, periodic_point,
                                       cycle_points, chi_of_word, is_expanding_word, preimages, Phi, Phi_star,
                                       atom, good_integer, n_primitive_necklaces, orbit_to_one, affine_of_word)


# ----------------------------------------------------------------------------------------------
# A. Liouville bound and heights of periodic points
# ----------------------------------------------------------------------------------------------
def part_A():
    say("== A. product-formula (Liouville) bound and heights of periodic points ==")
    rng = random.Random(20260926)
    worst = -1e9
    cases = 0
    for _ in range(4000):
        D = 2 * rng.randrange(0, 5000) + 1
        u = rng.randrange(-10 ** 6, 10 ** 6)
        beta = Fraction(u, D)
        if beta.denominator == 1 and beta > 0:
            continue
        for _ in range(3):
            n = rng.randrange(1, 10 ** 9)
            if Fraction(n) == beta:
                continue
            val = v2(Fraction(n) - beta)
            worst = max(worst, val - math.log2((n + 1) * height(beta)))
            cases += 1
        for Dd in (10, 30, 60):                    # integers deliberately close to beta
            y = good_integer(beta, Dd)
            val = v2(Fraction(y) - beta)
            worst = max(worst, val - math.log2((y + 1) * height(beta)))
            cases += 1
    check(worst <= 0, "Liouville bound v2(n - beta) <= log2((n+1) H(beta)) on %d (n, beta) pairs, "
          "including 12000 pairs with n 2-adically close to beta (max slack %.3f)" % (cases, worst))
    # heights of periodic points: H(x_w) <= 3^(a-1) 2^p, so log2 H(x_w) < p log2 6
    mx = 0.0
    nw = 0
    for p in range(1, 15):
        for m in range(1 << p):
            w = [(m >> i) & 1 for i in range(p)]
            a = sum(w)
            if a == 0:
                continue
            x = periodic_point(w)
            check_ok = height(x) <= 3 ** (a - 1) * 2 ** p
            if not check_ok:
                check(False, "height bound fails for word %s" % w)
            mx = max(mx, math.log2(height(x)) / p)
            nw += 1
    check(mx < math.log2(6), "H(x_w) <= 3^(a-1) 2^p for all %d words of length <= 14 with a >= 1; "
          "max log2 H / p = %.4f < log2 6 = %.4f" % (nw, mx, math.log2(6)))


# ----------------------------------------------------------------------------------------------
# B. Lemma L: local expansion
# ----------------------------------------------------------------------------------------------
def lemma_L_bound(zeta, D, bank, B=2):
    s = 0.0
    for beta, c in bank:
        if beta != zeta and v2(beta - zeta) >= D:
            s += abs(c) * (1 + math.log2(height(beta)))
    return (B + math.log2(height(zeta))) * s


def part_B():
    say("== B. Lemma L (local expansion of a bank at a rational point) ==")
    rng = random.Random(7)
    zetas = [Fraction(-1), Fraction(-2), Fraction(-5, 3), Fraction(-13, 9), Fraction(-5), Fraction(-17),
             periodic_point([1, 1, 1, 1, 0])]
    n_checks = 0
    max_ratio = 0.0
    for zeta in zetas:
        # a bank: an atom at zeta, the -1/-5/-17 cycle points, and "near atoms" zeta + 2^m/3 of height ~2^m
        bank = [(zeta, 0.7)]
        for w in ([1], [1, 1, 0]):
            for x in cycle_points(w):
                if x != zeta:
                    bank.append((x, rng.uniform(-1, 1)))
        for m in range(1, 241):
            bank.append((zeta + Fraction(2 ** m, 3), rng.uniform(-1, 1) / m ** 3))
        ps = Phi_star(zeta, bank)
        cz = atom(zeta, bank)
        for D in (8, 16, 32, 48, 64, 96, 160):
            y = good_integer(zeta, D)
            err = Phi(y, bank) - cz * D - ps
            bd = lemma_L_bound(zeta, D, bank)
            if abs(err) > bd + 1e-9:
                check(False, "Lemma L error bound fails at zeta=%s D=%d: |err|=%g > %g" % (zeta, D, abs(err), bd))
            if bd > 0:
                max_ratio = max(max_ratio, abs(err) / bd)
            n_checks += 1
        errs = []
        for D in (40, 80, 160):
            yb = good_integer(zeta, D)
            errs.append((abs(Phi(yb, bank) - cz * D - ps), lemma_L_bound(zeta, D, bank)))
        check(all(e <= b + 1e-12 for e, b in errs) and errs[0][1] > errs[1][1] > errs[2][1] > 0,
              "Lemma L at zeta = %s (signed bank with 240 near atoms of height ~2^m): |Phi(y) - c D - Phi*| = "
              "%.1e, %.1e, %.1e <= bound %.1e, %.1e, %.1e at D = 40, 80, 160" %
              (zeta, errs[0][0], errs[1][0], errs[2][0], errs[0][1], errs[1][1], errs[2][1]))
    check(n_checks == 49 and max_ratio <= 1.0,
          "Lemma L error bound respected in all %d (zeta, D) cases (max |err|/bound = %.3f)" % (n_checks, max_ratio))


# ----------------------------------------------------------------------------------------------
# C. Theorem A witnesses
# ----------------------------------------------------------------------------------------------
def R(n, bank):
    return math.log(n) + Phi(n, bank)


def backward_tree(root, depth):
    """points z with T^j z = root, j <= depth, as (z, j)"""
    out = [(Fraction(root), 0)]
    frontier = [Fraction(root)]
    for j in range(1, depth + 1):
        nxt = []
        for x in frontier:
            for y in preimages(x):
                if y != x and y not in [p for p, _ in out]:
                    nxt.append(y)
                    out.append((y, j))
        frontier = nxt
    return out


def predicted_violation(z, j, bank):
    """walk the chain z, Tz, ..., T^j z = -1 and return (step index i, slope) of the first charge increase,
    or ('shadow', kappa - c(-1)) if the charges are non-increasing and the -1 charge is below kappa"""
    chain = [Fraction(z)]
    for _ in range(j):
        chain.append(T(chain[-1]))
    assert chain[-1] == -1
    cs = [atom(x, bank) for x in chain]
    for i in range(j):
        if cs[i + 1] > cs[i] + 1e-12:
            return ("step", i, cs[i + 1] - cs[i], chain)
    if cs[-1] < KAPPA - 1e-12:
        return ("shadow", j, KAPPA - cs[-1], chain)
    return (None, None, 0.0, chain)


def increment_at(z, D, i, bank):
    """R(T^(i+1) y) - R(T^i y) for the good integer y at depth D from z"""
    y = good_integer(z, D)
    a = y
    for _ in range(i):
        a = T(a)
    b = T(a)
    return R(b, bank) - R(a, bank), a


def part_C():
    say("== C. Theorem A: forced charges on the backward orbit of -1; explicit violations ==")
    tree = backward_tree(-1, 6)
    check(len(tree) >= 20 and all(p < 0 for p, _ in tree),
          "the backward tree of -1 to depth 6 has %d points, all negative rationals" % len(tree))
    banks = {}
    banks["B1 {-1: 0.9 kappa}"] = [(Fraction(-1), 0.9 * KAPPA)]
    banks["B2 {-2^j: 1.2 kappa 2^-j, j<=30}"] = [(Fraction(-(2 ** j)), 1.2 * KAPPA / 2 ** j) for j in range(31)]
    banks["B3 {-2^j: 1.2 kappa, j<=20}"] = [(Fraction(-(2 ** j)), 1.2 * KAPPA) for j in range(21)]
    rng = random.Random(11)
    b4 = []
    for p in range(1, 9):
        for m in range(1 << p):
            w = [(m >> t) & 1 for t in range(p)]
            if is_expanding_word(w) and periodic_point(w).denominator % 2 == 1:
                x = periodic_point(w)
                b4.append((x, rng.uniform(0.5, 2.5) * KAPPA / p ** 2))
    for zz, _ in tree:
        if zz not in [b for b, _ in b4]:
            b4.append((zz, rng.uniform(0.2, 2.0) * KAPPA))
    # merge duplicates
    merged = {}
    for b, c in b4:
        merged[b] = merged.get(b, 0.0) + c
    banks["B4 random (periodic points p<=8, tree of -1 to depth 6)"] = list(merged.items())
    total_viol = 0
    for name, bank in banks.items():
        found = None
        for z, j in sorted(tree, key=lambda t: t[1]):
            pv = predicted_violation(z, j, bank)
            if pv[0] is not None:
                found = (z, j, pv)
                break
        if found is None:
            # B3 charges the even chain only: use the uncharged odd preimage -5/3 of -2
            z = Fraction(-5, 3)
            found = (z, 2, predicted_violation(z, 2, bank))
        z, j, (kind, i, slope, chain) = found
        if kind == "step":
            d1, _ = increment_at(z, 60, i, bank)
            d2, _ = increment_at(z, 120, i, bank)
            d3, _ = increment_at(z, 240, i, bank)
            ok = d1 > 0 and abs((d2 - d1) - slope * 60) < 1e-6 and abs((d3 - d2) - slope * 120) < 1e-6
            check(ok, "%s: the preimage z = %s of -1 (j = %d) has a charge increase at chain step %d (%s -> %s); "
                  "R(T^(i+1) y) - R(T^i y) = %.3f, %.3f, %.3f at depths 60, 120, 240: linear with slope %.5f = "
                  "c(T^(i+1) z) - c(T^i z)" % (name, z, j, i, chain[i], chain[i + 1], d1, d2, d3, slope))
        else:
            # in-shadow violation: along the shadow of -1 every step exceeds 0 by about kappa - c(-1)
            y = good_integer(Fraction(-1), 80)
            incs = []
            a = y
            for _ in range(40):
                b = T(a)
                incs.append(R(b, bank) - R(a, bank))
                a = b
            ok = min(incs) > 0 and abs(min(incs) - slope) < 1e-9 + 1.0 / y ** 0.5
            check(ok, "%s: charge at -1 below kappa, so every step inside the shadow of -1 increases R by "
                  "about kappa - c(-1) = %.5f (min over 40 steps at depth 80: %.5f)" % (name, slope, min(incs)))
        total_viol += 1
    check(total_viol == len(banks), "every one of the %d test banks has an explicit, predicted violation" % len(banks))
    # lookahead version: for the geometric bank B2, no j <= L gives R(T^j y) < R(y) at the entry y
    bank = banks["B2 {-2^j: 1.2 kappa 2^-j, j<=30}"]
    z = Fraction(-(2 ** 3))
    for L in (5, 20, 50):
        D = 4 * L + 40
        y = good_integer(z, D)
        base = R(y, bank)
        a = y
        m = float("inf")
        for _ in range(L):
            a = T(a)
            m = min(m, R(a, bank) - base)
        check(m > 0, "lookahead L = %d: from the good integer y at depth %d from -8 (bank B2), "
              "R(T^j y) - R(y) >= %.3f > 0 for every 1 <= j <= L" % (L, D, m))


# ----------------------------------------------------------------------------------------------
# D. sufficiency inside shadows (the tension on the deep cycle) and failure at the entry
# ----------------------------------------------------------------------------------------------
def T_sign(x, sgn):
    """shortcut map with a constant sign: x/2 or (3x + sgn)/2 (sgn = +1 Collatz, -1 the 3n-1 map)"""
    x = Fraction(x)
    return (3 * x + sgn) / 2 if x.numerator % 2 else x / 2


def cycle_points_sign(word, sgn):
    A, B, C = 1, 1, 0
    for b in word:
        if b:
            A, C = 3 * A, 3 * C + sgn * B
        B *= 2
    x = Fraction(C, B - A)
    pts = [x]
    for _ in range(len(word) - 1):
        pts.append(T_sign(pts[-1], sgn))
    assert T_sign(pts[-1], sgn) == x
    return pts


def preimages_sign(x, sgn):
    x = Fraction(x)
    out = [2 * x]
    y = (2 * x - sgn) / 3
    if y.numerator % 2:
        out.append(y)
    return out


def shadow_rank(word, c, sgn=1):
    """bank = the cycle points with weight c; h on the deep nodes = -c tau_j + g_j (note, Proposition 2.6)"""
    pts = cycle_points_sign(word, sgn)
    p = len(word)
    tau = [sum(v2(pts[j] - pts[i]) for i in range(p) if i != j) for j in range(p)]
    lam = sum(word) * LN3 - p * LN2
    mubar = (lam - p * c) / p
    g = [0.0]
    for j in range(p - 1):
        w = KAPPA if word[j] else -LN2
        g.append(g[-1] + mubar - w + c)
    K = max([v2(pts[i] - pts[j]) for i in range(p) for j in range(p) if i != j] + [0]) + 2
    hmap = {}
    for j in range(p):
        x = pts[j]
        r = (x.numerator * pow(x.denominator, -1, 1 << K)) % (1 << K)
        hmap[r] = -c * tau[j] + g[j]
    bank = [(x, c) for x in pts]

    def Rk(n):
        return math.log(n) + Phi(n, bank) + hmap.get(n % (1 << K), 0.0)
    return Rk, K, mubar, pts, bank


def Tint_sign(n, sgn):
    return (3 * n + sgn) >> 1 if n & 1 else n >> 1


def part_D():
    say("== D. sufficiency inside shadows; failure at the entry ==")
    w17 = []
    x = -17
    for _ in range(11):
        w17.append(x & 1)
        x = T(x)
    for word, name, sgn in (([1], "-1 (3n+1)", 1), ([1, 1, 0], "-5 cycle (3n+1)", 1), (w17, "-17 cycle (3n+1)", 1),
                            ([1, 1, 0], "5 -> 7 -> 10 cycle (3n-1)", -1)):
        chi = chi_of_word(word)
        c = chi + 0.02
        Rk, K, mubar, pts, bank = shadow_rank(word, c, sgn)
        worst = -1e9
        steps = 0
        for D in (K + 30, K + 90, K + 200):
            y = good_integer(pts[0], D)
            a = y
            for t in range(D - K - 1):
                b = Tint_sign(a, sgn)
                worst = max(worst, Rk(b) - Rk(a))
                a = b
                steps += 1
        check(worst < mubar / 2 < 0,
              "%s (p = %d, chi = %.5f): with charge c = chi + 0.02 on the %d cycle points and the periodic "
              "correction h = -c tau + g on the deep nodes mod 2^%d, every one of %d shadow steps has "
              "R(Tn) - R(n) <= %.5f < mubar/2 = %.5f" % (name, len(word), chi, len(pts), K, steps, worst, mubar / 2))
        ent = None
        for j, xj in enumerate(pts):
            for pre in preimages_sign(xj, sgn):
                if pre not in pts:
                    ent = (pre, xj)
                    break
            if ent:
                break
        pre, xj = ent
        incs = []
        for D in (60, 120, 240):
            y = good_integer(pre, D)
            incs.append(Rk(Tint_sign(y, sgn)) - Rk(y))
        check(incs[0] > 0 and abs((incs[1] - incs[0]) - c * 60) < 1e-6 and abs((incs[2] - incs[1]) - c * 120) < 1e-6,
              "%s: the entry step from the uncharged preimage %s of %s increases R by %.2f, %.2f, %.2f at depths "
              "60, 120, 240 (slope exactly c = %.5f): the shadow is paid, the seam is not" %
              (name, pre, xj, incs[0], incs[1], incs[2], c))


# ----------------------------------------------------------------------------------------------
# E. the forced mass per period
# ----------------------------------------------------------------------------------------------
def forced_mass(p):
    """M(p) = sum over expanding primitive necklaces of length p of log(3^a / 2^p) (= sum of the forced
    charges chi over their p points); also the number of expanding periodic points of exact period p"""
    M = 0.0
    npts = 0
    for a in range(p + 1):
        if 3 ** a > 2 ** p:
            k = n_primitive_necklaces(p, a)
            M += k * (a * LN3 - p * LN2)
            npts += k * p
    return M, npts


def part_E():
    say("== E. the forced charges are not summable: mass per period ==")
    say("   p   expanding points P(p)   forced mass M(p)   M(p) / (2^(h p) p^(-3/2))")
    rows = []
    for p in list(range(1, 13)) + [16, 20, 24, 32, 40, 48, 64]:
        M, npts = forced_mass(p)
        ratio = M / (2 ** (H_BITS * p) * p ** -1.5)
        rows.append((p, npts, M, ratio))
        say("  %3d   %22d   %16.6g   %.4f" % (p, npts, M, ratio))
    lo_ok = all(abs(math.log2(r[2]) - H_BITS * r[0]) <= 2 * math.log2(r[0]) + 2 for r in rows if r[0] >= 16)
    check(all(r[2] > 0 for r in rows if r[0] != 2) and rows[1][2] == 0 and rows[-1][2] > 1e15 and lo_ok,
          "M(p) > 0 for every p != 2 (M(2) = 0: no expanding primitive word of length 2) and M(64) = %.3g; "
          "|log2 M(p) - h p| <= 2 log2 p + 2 for p = 16..64, as the binomial bounds M(p) = 2^(h p + O(log p)) predict: "
          "the forced charges are not summable" % rows[-1][2])
    rr = [r[3] for r in rows if r[0] >= 16]
    check(max(rr) / min(rr) < 3.0,
          "M(p) / (2^(h p) p^(-3/2)) stays in [%.3f, %.3f] for p = 16..64 (h = %.5f), consistent with the local-CLT "
          "order 2^(h p) p^(-3/2) (heuristic, not claimed)" % (min(rr), max(rr), H_BITS))
    # the words 1^(L-1) 0 alone: chi -> log(3/2), total forced charge per orbit log lambda_L -> infinity
    vals = [((L - 1) * LN3 - L * LN2) for L in (3, 10, 30, 100)]
    check(all(v2_ > v1 for v1, v2_ in zip(vals, vals[1:])) and vals[-1] > 39,
          "the orbits O_L of 1^(L-1)0 alone need total charge log lambda_L = %s (L = 3, 10, 30, 100)" %
          ", ".join("%.3f" % v for v in vals))


# ----------------------------------------------------------------------------------------------
# F. Theorem D on a window
# ----------------------------------------------------------------------------------------------
def part_F(N0=1000):
    say("== F. Theorem D: the private-atom bank (from stopping times) on a finite window ==")
    # window: all values on the orbits of 2..N0
    vals = set()
    for n in range(2, N0 + 1):
        vals.update(orbit_to_one(n))
    Mx = max(vals)
    check(Mx < 10 ** 6, "the orbits of 2..%d (shortcut map) stay below M = %d" % (N0, Mx))
    import numpy as np
    M = Mx
    # stopping times sigma(m) for m <= M (shortcut map), by memoised iteration
    sig = np.full(M + 1, -1, dtype=np.int64)
    sig[1] = 0
    for m in range(2, M + 1):
        if sig[m] >= 0:
            continue
        path = []
        x = m
        while x > M or sig[x] < 0:
            path.append(x)
            x = (3 * x + 1) >> 1 if x & 1 else x >> 1
        s = int(sig[x])
        for y in reversed(path):
            s += 1
            if y <= M:
                sig[y] = s
    a = 1.0
    C = 2 * a
    delta = 0.05
    K = a + 1.0 + 10 * delta
    m_idx = np.arange(1, M + 1, dtype=np.float64)
    eps = delta / (m_idx * (m_idx + 1))                         # eps[m-1] = delta/(m(m+1))
    F = C * sig[1:].astype(np.float64) - a * np.log(m_idx) + K
    # leak L(n) = sum_{m != n, m <= M} eps_m v2(n - m) (+ a tail <= sum_{m > M} eps_m log2(4 n m) <= tail_bd)
    # computed by the 2-adic tree: sum_m eps_m v2(n - m) = sum_s [W_s(n mod 2^s) - eps_n]
    Lk = np.zeros(M, dtype=np.float64)
    s = 1
    nvals = np.arange(1, M + 1, dtype=np.int64)
    while (1 << s) <= 2 * M:
        mod = 1 << s
        W = np.bincount(nvals % mod, weights=eps, minlength=mod)
        Lk += W[nvals % mod] - eps
        s += 1
    Dn = np.floor((F - Lk) / eps)
    check(bool(np.all(Dn >= np.log2(m_idx) + 2)),
          "construction admissible on [1, %d]: every D_m >= log2 m + 2 (min D_m = %.3g)" % (M, Dn.min()))
    check(bool(np.all(Dn > math.log2(M) + 2)),
          "every D_m exceeds log2 M + 2, so v2(n - z_m) = v2(n - m) for all n != m in the window")
    Phi_w = eps * Dn + Lk            # the exact potential at n <= M of the finite bank {(z_m, eps_m) : m <= M}
    err = Phi_w - F
    check(bool(np.all(err <= 1e-9)) and bool(np.all(err > -eps - 1e-9)),
          "the finite bank {z_m = m + 2^(D_m)/3, eps_m} (m <= %d) has -eps_n < Phi(n) - F(n) <= 0 on the window" % M)
    Rw = a * np.log(m_idx) + Phi_w
    worst = -1e9
    nsteps = 0
    for n in range(2, N0 + 1):
        x = n
        while x != 1:
            y = (3 * x + 1) >> 1 if x & 1 else x >> 1
            d = Rw[y - 1] - Rw[x - 1]
            worst = max(worst, d)
            nsteps += 1
            x = y
    eta = C - delta / 2
    check(worst <= -eta + 1e-9,
          "Theorem D window check: along all %d steps of the orbits of 2..%d, R(Tn) - R(n) <= %.4f "
          "<= -C + delta/2 = %.4f (the finite bank is exact on the window; the infinite construction adds a "
          "tail in [0, delta((3 + log2 M)/M + 1/(M ln 2))])" % (nsteps, N0, worst, -eta))
    # the bank is massively non-height-summable: sum eps_m log2 H(z_m) ~ sum eps_m D_m ~ sum F(m)
    hm = [float((eps[:mm] * Dn[:mm]).sum()) for mm in (1000, 10000, 100000, M)]
    check(hm[0] < hm[1] < hm[2] < hm[3] and hm[3] > 1e6,
          "height moment sum_{m<=M} eps_m D_m = %s for M = 10^3, 10^4, 10^5, %d: it diverges (Theorem C)" %
          (", ".join("%.3g" % v for v in hm), M))
    # the heavy preimage chains: Phi(y) >= kappa (V - 1) - j log 2 at y = 2^j (2^V - 1) in the window
    ok = True
    cnt = 0
    for V in range(3, 17):
        for j in range(0, 18 - V):
            y = (2 ** V - 1) * 2 ** j
            if y <= M:
                ok &= Phi_w[y - 1] + 1e-9 >= KAPPA * (V - 1) - j * LN2
                cnt += 1
    check(ok, "the shadow lower bound Phi(y) >= kappa (V-1) - j log 2 holds at all %d points y = 2^j (2^V - 1) "
          "of the window" % cnt)


def run():
    part_A()
    part_B()
    part_C()
    part_D()
    part_E()
    part_F()


if __name__ == "__main__":
    run()
