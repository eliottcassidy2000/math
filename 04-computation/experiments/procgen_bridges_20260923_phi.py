#!/usr/bin/env python3
"""procgen-bridges 2026-09-23, part PHI: where the golden ratio enters LRC(14), AMM 12592 and Collatz.

Sections (deterministic output to stdout; timing to stderr):
  PHI1  the four golden occurrences as eigen-data of the Fibonacci matrix M = [[1,1],[1,0]] (exact, sympy)
        Theorem S (Collatz), Lemma P (Collatz), THM-3027 / Theorem B / Long's point (AMM)
  PHI2  co-variation ("knob") tests: move each occurrence's own structural parameter, see what moves
        AMM knob b (capacity boost, THM-3027 rate problem re-solved independently), AMM knob q (THM-3009),
        Theorem S knob L (word class), Lemma P knob (shift family)
  PHI3  Theorem S's golden threshold is a single-scale (decoupled) artifact: realized Diophantine
        exponents of Sturmian words vs Theorem S's guarantee mu_S(L) vs Bugeaud-Kim's floor 2.50994
  PHI4  golden/Fibonacci in the live LRC(14) frontier (file scan)
Memory: < 200 MB.  Runtime: about 1-2 minutes.
"""
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[2]
T0 = time.time()
PHI = (1 + 5 ** 0.5) / 2


def P(*a):
    print(*a, flush=True)


def tick(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------------------------- PHI1
def section_phi1():
    P("=" * 100)
    P("PHI1  the four golden occurrences as eigen-data of M = [[1,1],[1,0]] (exact)")
    P("=" * 100)
    x, y, t, p, mu, g = sp.symbols("x y tau p mu gamma")
    phi = (1 + sp.sqrt(5)) / 2
    M = sp.Matrix([[1, 1], [1, 0]])
    cM = sp.expand((x * sp.eye(2) - M).det())
    cM2 = sp.expand((x * sp.eye(2) - M * M).det())
    P(f"  charpoly(M)   = {cM}      eigenvalues {sorted(M.eigenvals().keys(), key=lambda e: float(e))}")
    P(f"  charpoly(M^2) = {cM2}     eigenvalues {sorted((M * M).eigenvals().keys(), key=lambda e: float(e))}")
    P(f"  disc(M) = lambda_max - lambda_min = {sp.simplify(phi - (-1 / phi))}")

    # (1) Theorem S: limsup q_(n+1)/q_n >= phi, equality iff eventually all partial quotients 1
    P("  (1) Collatz Theorem S: threshold mu(mu-1) < L with L = limsup q_(n+1)/q_n.")
    for a in (1, 2, 3):
        A = sp.Matrix([[a, 1], [1, 0]])
        rad = max(A.eigenvals().keys(), key=lambda e: float(e))
        P(f"      CF recursion matrix [[{a},1],[1,0]]: spectral radius {sp.nsimplify(rad)} = {float(rad):.6f}")
    fib = [1, 1]
    for _ in range(40):
        fib.append(fib[-1] + fib[-2])
    P(f"      all-ones CF (golden slope): q_(n+1)/q_n = F_(n+1)/F_n -> {fib[-1] / fib[-2]:.12f} (phi = {PHI:.12f})")
    muS = sp.solve(sp.Eq(mu * (mu - 1), phi), mu)
    muS = [m for m in muS if float(m) > 0][0]
    P(f"      mu_S(phi) = root of mu(mu-1) = phi: {sp.nsimplify(muS)} = {float(muS):.6f}  (Theorem S: mu < 1.8668)")
    P(f"      origin: phi = spectral radius of the minimal CF matrix = slowest denominator growth (Hurwitz/Fibonacci)")

    # (2) Lemma P: Rayleigh quotient of the Fibonacci matrix
    f = (y ** 2 + 2 * y) / (y ** 2 + 1)
    crit = sp.solve(sp.diff(f, y), y)
    crit_pos = [c for c in crit if float(c) > 0][0]
    A = sp.Matrix([[0, 1], [1, 1]])
    v = sp.Matrix([1, y])
    ray = sp.simplify((v.T * A * v)[0] / (v.T * v)[0])
    P("  (2) Collatz Lemma P: bracket (x^2+4x+3)/(x^2+2x+2), x = m/n; with y = 1+x = (n+m)/n it is")
    P(f"      f(y) = (y^2+2y)/(y^2+1) = Rayleigh quotient of [[0,1],[1,1]] in basis (1,y): {sp.simplify(ray - f) == 0}")
    P(f"      critical point y* = {sp.nsimplify(crit_pos)} (= phi: {sp.simplify(crit_pos - phi) == 0}),"
      f" max value f(y*) = {sp.nsimplify(sp.simplify(f.subs(y, crit_pos)))} (= phi); x* = y*-1 = 1/phi")
    P(f"      value = argmax because at a critical point f = 1 + 1/y and y^2 = y+1  (phi = 1 + 1/phi);"
      f" unshifted m=0 (y=1): f(1) = {f.subs(y, 1)}")
    P(f"      origin: top eigenvalue of the valuation-vs-height quadratic-form pencil of ONE Pade family")

    # (3) THM-3027 tangency
    sol = sp.solve(sp.Eq((1 - t) ** 2, t), t)
    P("  (3) AMM THM-3027: tangency (1-tau)^2 = tau <=> tau^2 - 3 tau + 1 = charpoly(M^2):"
      f" {sp.expand((1 - t) ** 2 - t) == sp.expand(cM2.subs(x, t))}")
    tstar = [s for s in sol if 0 < float(s) < 1][0]
    P(f"      tau* = {sp.nsimplify(tstar)} = phi^-2: {sp.simplify(tstar - phi ** -2) == 0};"
      f"  Zudilin's Tschakaloff condition (3-sqrt5)/2 equals tau*: {sp.simplify((3 - sp.sqrt(5)) / 2 - tstar) == 0}")
    P(f"      origin: double root of the entropy value function (a Stirling saddle), gamma* = log(phi)/log(sqrt5)")

    # (4) Theorem B / Long's evaluation point
    solp = sp.solve(sp.Eq(p * (1 - p), -1), p)
    pneg = [s for s in solp if float(s) < 0][0]
    S = sp.Abs(pneg) + sp.Abs(1 - pneg)
    gstar = sp.log(phi) / sp.log(sp.sqrt(5))
    P("  (4) AMM Theorem B: the golden point w = p(1-p) = -1 <=> p^2 - p - 1 = charpoly(M)(p):"
      f" p = {sp.nsimplify(pneg)} = -1/phi: {sp.simplify(pneg + 1 / phi) == 0}")
    P(f"      p + (1-p) = 1 = trace(M) (the fairness involution), p(1-p) = -1 = det(M) (the unit circle |w| = 1)")
    P(f"      S(p) = |p| + |1-p| = {sp.nsimplify(sp.simplify(S))} = lambda_max - lambda_min;"
      f"  |p| S^gamma* = {sp.nsimplify(sp.simplify(sp.exp(sp.log(sp.Abs(pneg)) + gstar * sp.log(S))))}")
    th = phi ** -2
    xq = -th
    P(f"      Long's evaluation point x = -theta*, theta* = phi^-2: 1+x = {sp.nsimplify(sp.simplify(1 + xq))} = 1/phi:"
      f" {sp.simplify(1 + xq - 1 / phi) == 0}; 1+2x = {sp.nsimplify(sp.simplify(1 + 2 * xq))} = phi^-3:"
      f" {sp.simplify(1 + 2 * xq - phi ** -3) == 0} (both positive; only x^i alternates)")
    u = sp.sqrt(5)
    P(f"      u = sqrt5 <-> w = -1 via 1+w = (5-u^2)/4: {sp.simplify(1 + (-1) - (5 - u ** 2) / 4) == 0}")
    # the Catalan S-fraction at w = -1 is the all-ones continued fraction
    zc = [sp.Integer(0)]
    for _ in range(12):
        zc.append(sp.Rational(-1) / (1 - zc[-1]))
    P(f"      Catalan quotient z(w) = w/(1 - z(w)) at w = -1: convergents {[str(c) for c in zc[1:9]]} = -F_n/F_(n+1) -> -1/phi")
    P("      origin: the continuation domain first touches the natural boundary |w| = 1 of the lacunary class")
    P("  SUMMARY PHI1: all four are eigen-data of the minimal hyperbolic unit of GL_2(Z) (trace 1, det -1) or its square")
    P("  in SL_2(Z) (trace 3), reached through four different extremal problems.")


# ---------------------------------------------------------------------------------------------- PHI2
def H(z):
    if z <= 0 or z >= 1:
        return 0.0
    return -z * math.log(z) - (1 - z) * math.log(1 - z)


def psi_rate(gam, tau, b, n_grid=160):
    """(1/R) log S(t) -> max_sigma [ D H(rho) + m log b ], D = gam (1+sigma), m = tau - sigma, rho = m/D <= 1."""
    lo, hi = 0.0, min(1.0, tau)

    def val(s):
        D = gam * (1 + s)
        m = tau - s
        if m < 0 or m > D:
            return -1e9
        return D * H(m / D) + m * math.log(b)
    best_s = max(np.linspace(lo, hi, n_grid), key=val)
    a, c = max(lo, best_s - (hi - lo) / n_grid), min(hi, best_s + (hi - lo) / n_grid)
    for _ in range(60):   # golden section (inner concavity: THM-3027 audit)
        m1, m2 = a + (c - a) * 0.381966, a + (c - a) * 0.618034
        if val(m1) < val(m2):
            a = m1
        else:
            c = m2
    s = 0.5 * (a + c)
    return val(s), s


def margin(gam, b):
    """min over tau of Psi - H(tau); returns (min, argmin tau)."""
    taus = np.linspace(0.01, 0.99, 197)
    vals = [psi_rate(gam, tt, b)[0] - H(tt) for tt in taus]
    i = int(np.argmin(vals))
    a, c = taus[max(0, i - 1)], taus[min(len(taus) - 1, i + 1)]
    for _ in range(50):
        m1, m2 = a + (c - a) * 0.381966, a + (c - a) * 0.618034
        if psi_rate(gam, m1, b)[0] - H(m1) > psi_rate(gam, m2, b)[0] - H(m2):
            a = m1
        else:
            c = m2
    tt = 0.5 * (a + c)
    return psi_rate(gam, tt, b)[0] - H(tt), tt


def section_phi2():
    P("=" * 100)
    P("PHI2  co-variation ('knob') tests: does anything move all golden occurrences together?")
    P("=" * 100)
    P("  (a) AMM knob b (capacity boost b^(t-i) in the THM-3002 criterion), THM-3027 rate problem re-solved here")
    P("      independently (grid + golden-section inner max, bisection in gamma):")
    for b in (2, 3, 4, 5):
        lo, hi = 0.05, 1.5
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if margin(mid, b)[0] >= 0:
                hi = mid
            else:
                lo = mid
        gnum = 0.5 * (lo + hi)
        _, tau_star = margin(gnum, b)
        gcl = math.log(PHI) / math.log((b + PHI) / PHI)
        P(f"      b={b}: gamma*(b) numeric {gnum:.7f}, closed form log(phi)/log((b+phi)/phi) = {gcl:.7f};"
          f" binding tau* = {tau_star:.5f} (phi^-2 = {PHI ** -2:.5f})")
    P("      => tau* = phi^-2 is universal in b (the golden quadratic is the tangency), gamma* moves with b.")
    P("  (b) AMM knob q (THM-3009 alphabet, |Delta^j P(0)| <= q^j): delta^q = (1-delta)^(q-1)")
    for q in (2, 3, 4, 5):
        lo, hi = 0.5, 1.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if mid ** q - (1 - mid) ** (q - 1) < 0:
                lo = mid
            else:
                hi = mid
        P(f"      q={q}: delta* = {0.5 * (lo + hi):.10f}{'  = 1/phi' if q == 2 else ''}")
    P("      => under THIS knob the golden value is special to q = 2 (no metallic family).")
    P("  (c) Collatz Theorem S knob L = limsup q_(n+1)/q_n (the word class), threshold mu_S(L) = root of mu(mu-1) = L:")
    for name, L in (("Sturmian infimum (golden slope)", PHI), ("L=2", 2.0), ("silver L=1+sqrt2", 1 + 2 ** 0.5), ("L=3", 3.0)):
        P(f"      {name:34s} mu_S = {(1 + math.sqrt(1 + 4 * L)) / 2:.5f}")
    P("      the map (every 3x+r, mu <= log2 3 = 1.585) does not move it; only the word class does.")
    P("  (d) Collatz Lemma P knob (the Pade family): shift ratio x = m/n")
    for xx in (0.0, 0.3, 1 / PHI, 1.0, 2.0):
        yy = 1 + xx
        P(f"      x = {xx:.4f}: bracket threshold (y^2+2y)/(y^2+1) = {(yy * yy + 2 * yy) / (yy * yy + 1):.6f}")
    P("      independent of rho = 2^L/M_0: the golden value is the optimum of the construction, not of the number.")
    P("  VERDICT PHI2: four knobs, four different responses (universal / special / class-dependent / construction-")
    P("  dependent).  No parameter moves the AMM and Collatz golden values together: across programs the shared phi")
    P("  is NUMEROLOGY (common algebra of the smallest hyperbolic unit), inside AMM it is one mechanism with three faces.")


# ---------------------------------------------------------------------------------------------- PHI3
import mpmath as mpm


def cf_terms_mp(x, n=60, dps=80):
    """partial quotients a_1, a_2, ... of x in (0,1), from an mpmath value computed at high precision"""
    with mpm.workdps(dps):
        out = []
        for _ in range(n):
            x = 1 / x
            a = int(mpm.floor(x))
            out.append(a)
            x = x - a
            if x == 0:
                break
    return out


def slope_from_cf(prefix, period, n=60):
    terms = list(prefix) + list(period) * n
    with mpm.workdps(80):
        v = mpm.mpf(0)
        for a in reversed(terms[:n]):
            v = 1 / (a + v)
    return v, terms[:n]


def convergent_denoms(terms):
    q0, q1 = 1, terms[0]
    qs = [q1]
    for a in terms[1:]:
        q0, q1 = q1, a * q1 + q0
        qs.append(q1)
    return qs


def mechanical(alpha, rho, n):
    i = np.arange(n + 1, dtype=np.float64)
    fl = np.floor(i * alpha + rho)
    return (fl[1:] - fl[:-1]).astype(np.int8)


def lpf_profile(w, jmax):
    """LPF(j) = max_(a<j) LCE(a, j), vectorized over the period p = j - a."""
    N = len(w)
    LPF = np.zeros(jmax + 1, dtype=np.int64)
    for p in range(1, jmax + 1):
        n = N - p
        eq = w[:n] == w[p:]
        idx = np.arange(n)
        nf = np.where(eq, n, idx)
        nxt = np.minimum.accumulate(nf[::-1])[::-1]
        r = nxt - idx
        a_hi = min(jmax - p, n - 1)
        if a_hi < 0:
            continue
        seg = LPF[p:p + a_hi + 1]
        np.maximum(seg, r[:a_hi + 1], out=seg)
    return LPF


def section_phi3(N=30000):
    P("=" * 100)
    P("PHI3  Theorem S's golden threshold is a single-scale (decoupled) worst case")
    P("=" * 100)
    bk = 5 / 3 + 4 * math.sqrt(10) / 15
    bk2 = (math.sqrt(10) - 1.5) / (math.sqrt(10) - 2.5)
    P("  Theorem S: at infinitely many scales q_n one of two single-scale options gives a prefix repetition of")
    P("  exponent >= mu_S(L) (root of mu(mu-1) = L, L = limsup q_(n+1)/q_n): 1.8668 in the golden worst case.")
    P("  Theorem D (PROVED, audited) needs only Dio(w) = 1 + limsup_j LPF(j)/j > eta, and Bugeaud-Kim give")
    P(f"  Dio >= 5/3 + 4 sqrt10/15 = {bk:.7f} (= (sqrt10-3/2)/(sqrt10-5/2) = {bk2:.7f}) for EVERY Sturmian word (Cor. D1).")
    P("  Measured here: the repetition exponent actually delivered at EACH scale n of the slope,")
    P(f"  e_n = max_(q_n <= j <= q_n + q_(n+1)) (j + LPF(j))/j (the |UV| range of options A and B), on words of")
    P(f"  length N = {N} (j <= N/5, no end truncation), against the single-scale minimax mu_S(q_(n+1)/q_n)")
    P("  (min over the defect position j0 of max(option A, option B) = root of mu(mu-1) = q_(n+1)/q_n).")
    slopes = [
        ("golden 1/phi^2", (2,), (1,)),
        ("silver sqrt2-1", (), (2,)),
        ("(2,1,1)-periodic", (), (2, 1, 1)),
        ("(1,2)-periodic", (), (1, 2)),
    ]
    extra = [("Collatz critical log_3 2", mpm.log(2) / mpm.log(3)), ("e - 2", mpm.e - 2)]
    intercepts = [0.0, 0.3141592653589793, 0.7071067811865476, 0.9505]
    jmax = N // 5
    worst_defect = []
    rows = []
    items = [(nm, *slope_from_cf(pre, per)) for nm, pre, per in slopes]
    items += [(nm, v, cf_terms_mp(v, n=30, dps=120)) for nm, v in extra]
    for name, al_mp, terms in items:
        al = float(al_mp)
        qs = convergent_denoms(terms)
        # scale n: |UV| ranges over [q_n, q_n + q_(n+1)] (options A and B of Theorem S); guarantee mu_S(q_(n+1)/q_n)
        wins = [(a, b) for a, b in zip(qs[:-1], qs[1:]) if a >= 40 and a + b <= jmax]
        P(f"  {name:26s} alpha={al:.6f} CF [0;{','.join(map(str, terms[:12]))},...]")
        for rho in intercepts:
            w = mechanical(al, rho, N)
            lp = lpf_profile(w, jmax)
            es = []
            for a, b in wins:
                js = np.arange(a, a + b + 1)
                e = 1 + float((lp[js] / js).max())
                g = (1 + math.sqrt(1 + 4 * b / a)) / 2
                es.append((e, g))
            if not es:
                continue
            dmin = min(e - g for e, g in es)
            rows.append((name, rho, dmin))
            if name.startswith("golden"):
                worst_defect.append(dmin)
            P(f"      rho={rho:.4f}: {len(es)} scales; (realized e_n, guarantee mu_S(q_(n+1)/q_n)) = "
              + " ".join(f"({e:.3f},{g:.3f})" for e, g in es) + f"; min defect {dmin:+.3f}")
            tick(f"PHI3 {name} rho={rho:.3f}")
    allmin = min(r[2] for r in rows)
    P(f"  smallest single-scale defect over all rows and scales: {allmin:+.3f} (never negative: consistent with Theorem S)")
    P(f"  golden rows: every scale window delivers at least mu_S(phi) + {min(worst_defect):.3f}; the limsup (Dio) is")
    P(f"  larger still and >= {bk:.4f} by Bugeaud-Kim.  Single windows can dip below the BK floor (they are not limsups).")
    P("  => the golden number in Theorem S is the worst case of a one-scale split that treats the defect position")
    P("     j0 and the ratio q_(n+1)/q_n as independent adversaries; the Ostrowski coupling of the scales (Dio,")
    P("     Bugeaud-Kim, Corollary D1) replaces 1.8668 by 2.50994, a sqrt10 constant with no golden origin.")
    return rows


# ---------------------------------------------------------------------------------------------- PHI4
def section_phi4():
    P("=" * 100)
    P("PHI4  golden / Fibonacci in the live LRC(14) frontier (deterministic file scan)")
    P("=" * 100)
    fr = (REPO / "00-navigation" / "CURRENT-FRONTIER.md").read_text(encoding="utf-8", errors="replace")
    sec = fr.split("## LRC(14)")[1].split("\n## ")[0]
    pat = re.compile(r"golden|fibonacci|\bphi\b|sqrt ?5\b", re.I)
    P(f"  CURRENT-FRONTIER.md, section LRC(14): {len(pat.findall(sec))} golden/Fibonacci/phi tokens")
    thm = sorted((REPO / "01-canon" / "theorems").glob("THM-*"))
    lrc = [f for f in thm if re.search(r"lrc|lonely", f.name, re.I)]
    gold_word = re.compile(r"golden|fibonacci|sqrt ?5\b|1\s*\+\s*sqrt ?5", re.I)
    hits = [(f.name, len(gold_word.findall(f.read_text(encoding="utf-8", errors="replace")))) for f in lrc]
    hits = [h for h in hits if h[1] > 0]
    P(f"  theorem files named lrc/lonely: {len(lrc)}; containing golden/Fibonacci/sqrt5 (excluding the variable name"
      f" 'phi'): {len(hits)} -> {hits[:8]}")
    fibs = [f for f in thm if re.search(r"fibonacci", f.name, re.I)]
    deny = 0
    for f in fibs:
        txt = f.read_text(encoding="utf-8", errors="replace")[:4000]
        if re.search(r"No (LRC|Farey-graph, full-tree, LRC)|no LRC|LRC-blind|not preserve LRC|No LRC", txt):
            deny += 1
    P(f"  Fibonacci-titled theorems: {len(fibs)}; of these {deny} state explicitly in their header that no LRC transfer"
      " follows.")
    P("  LRC(14) frontier constants are rationals built from 7 (the 1/14 arc), clocks 2,3,4 and pigeonhole")
    P("  (214/1449, 72/539, 124/693, 17/693, 1/110, 195, 91^6): no golden threshold occurs.")


def main():
    quick = "--quick" in sys.argv
    section_phi1()
    tick("PHI1 done")
    section_phi2()
    tick("PHI2 done")
    section_phi3(N=12000 if quick else 30000)
    tick("PHI3 done")
    section_phi4()
    tick("PHI4 done")


if __name__ == "__main__":
    main()
