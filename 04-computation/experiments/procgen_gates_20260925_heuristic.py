#!/usr/bin/env python3
"""procgen_gates_20260925_heuristic.py -- PART C of the cycle-gate equidistribution lane.

C1  heuristic expected numbers of sporadic integral necklaces (primitive, |G| > 1, a >= 1):
      naive    E = sum L(p,a)/|G|                     (residues mod G uniform, independent);
      refined  E = sum L_elig(p,a)/|G|                (only necklaces whose rational cycle has least
                                                       point >= 1 can be integral: perigee bounds);
    for q = 3 and q = 5, dyadic and q-adic sides separately, with partial sums and tails; the
    eligible counts are exact for p <= 30 (part B dump), Monte Carlo on the transition window
    30 < p <= 300 (validated against the exact counts), and bracketed beyond;
C2  comparison with the census of part A (Poisson probabilities);
C3  the bound hierarchy on representative clocks: perigee (Belaga) << size range << sqrt(C)
    (second-moment barrier) << C, and the main term;
C4  the archimedean limit law: S(h)/C at fixed a/p for growing p (normalised DP) against Monte
    Carlo estimates of E e(h Y) for the limit series Y; near-critical clocks.
stdout = results, stderr = timing/memory.  Checks raise on failure.
"""
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_gates_20260925_core import S_dp, cmin_cmax, gate, lean_malloc, lyndon_count, report_mem  # noqa: E402
from procgen_gates_20260925_stats import elig_class  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ELIG_TSV = os.path.join(HERE, "..", "..", "scratch", "procgen_gates", "partB_elig.tsv")
P_EXACT = 30
P_MC = 1000            # transition-window Monte Carlo up to this period (band next to the critical line)
P_MAX = 4000
MC_SAMPLES = 400
RNG = np.random.default_rng(20260925)


def log_lyndon(p, a):
    if p <= 60:
        return math.log(lyndon_count(p, a))
    return math.lgamma(p + 1) - math.lgamma(a + 1) - math.lgamma(p - a + 1) - math.log(p)


def load_exact_elig():
    d = {}
    with open(ELIG_TSV) as fh:
        for line in fh:
            q, p, a, ne, ew, C = map(int, line.split())
            d[(q, p, a)] = (ne, ew, C)
    return d


def elig_class_fast(p, a, q):
    """same classification as elig_class, with float logs away from the boundaries."""
    l2, lq, lq1 = math.log(2), math.log(q), math.log(q + 1)
    d = p * l2 - a * lq
    if abs(d) < 1e-6:
        return elig_class(p, a, q)
    if d > 0:
        e1 = p * l2 - a * lq1  # 'none' iff 2^p > (q+1)^a
        if abs(e1) < 1e-9 or abs((a - 1) * lq - (p * l2 + math.log1p(-math.exp(-d)))) < 1e-9:
            return elig_class(p, a, q)
        if e1 > 0:
            return "none"
        if (a - 1) * lq >= p * l2 + math.log1p(-math.exp(-d)):
            return "all"
        return "transition"
    return elig_class(p, a, q) if p <= 400 else ("all" if q == 3 else _q_adic_fast(p, a, q, d))


def _q_adic_fast(p, a, q, d):
    l2, lq = math.log(2), math.log(q)
    if p * l2 < a * math.log(q - 1):
        return "none"
    lM = a * lq + math.log1p(-math.exp(d))
    lcmin = a * lq + math.log1p(-math.exp(a * (l2 - lq))) - math.log(q - 2)
    if (a - 1) * lq >= lM or lcmin >= lM:
        return "all"
    return "transition"


def mc_elig_fraction(p, a, q, S=MC_SAMPLES):
    """Monte Carlo share of words whose rational cycle has least point >= 1 (uniform word in W(p,a)).
    Dyadic: forward iteration of x -> x/2, (qx+1)/2 along the word (attracting, multiplier q^a/2^p < 1);
    q-adic: backward iteration y -> 2y, (2y+1)/q (the inverse branches; contracting over a period)."""
    G = gate(p, a, q)
    M = abs(G)
    keys = RNG.random((S, p))
    pos = np.sort(np.argsort(keys, axis=1)[:, :a], axis=1)  # uniform a-subsets
    i = np.arange(a)
    lt = (a - 1 - i)[None, :] * math.log(q) + pos * math.log(2) - math.log(M)
    mx = lt.max(axis=1, keepdims=True)
    x = np.exp(mx[:, 0]) * np.exp(lt - mx).sum(axis=1)  # |x_w| = c_w / M
    W = np.zeros((S, p), dtype=bool)
    np.put_along_axis(W, pos, True, axis=1)
    mn = x.copy()
    if G > 0:
        for t in range(p):
            x = np.where(W[:, t], (q * x + 1.0) / 2.0, x / 2.0)
            np.minimum(mn, x, out=mn)
    else:
        for t in range(p - 1, -1, -1):
            x = np.where(W[:, t], (2.0 * x + 1.0) / q, 2.0 * x)
            np.minimum(mn, x, out=mn)
    return float(np.count_nonzero(mn >= 1.0 - 1e-12)) / S


def c1_c2():
    exact = load_exact_elig()
    print("=" * 110)
    print("C1. Heuristic expected numbers of sporadic integral necklaces (primitive words, a >= 1, |G| > 1)")
    print("=" * 110)
    # validation of the Monte Carlo eligibility against exact counts (q = 3, transition clocks p <= 30)
    errs = []
    for (q, p, a), (ne, ew, C) in sorted(exact.items()):
        if q == 3 and elig_class(p, a, q) == "transition" and p >= 20:
            f_ex = ew / C
            f_mc = mc_elig_fraction(p, a, q, S=2000)
            errs.append(abs(f_ex - f_mc))
    print(f"  Monte Carlo eligibility validated on {len(errs)} transition clocks 20 <= p <= 30 (2000 samples): "
          f"max |exact - MC| = {max(errs):.3f}, mean {np.mean(errs):.4f}")
    assert max(errs) < 0.06
    results = {}
    for q in (3, 5):
        t0 = time.time()
        cum = {s: {"naive": 0.0, "ref_lo": 0.0, "ref_hi": 0.0} for s in (1, -1)}
        marks = {}
        per_p = []
        trans_mc = 0
        trans = []
        unsampled = {1: 0.0, -1: 0.0}
        l2, lq = math.log(2), math.log(q)
        for p in range(1, P_MAX + 1):
            row = {s: [0.0, 0.0, 0.0] for s in (1, -1)}
            lfp = math.lgamma(p + 1)
            for a in range(1, p + 1):
                d = p * l2 - a * lq
                # log|G| in floats (exact big integers near the critical line)
                if abs(d) < 0.01:
                    G = gate(p, a, q)
                    if abs(G) <= 1:
                        continue
                    lg = math.log(abs(G))
                    sgn = 1 if G > 0 else -1
                elif d > 0:
                    lg = p * l2 + math.log1p(-math.exp(-d))
                    sgn = 1
                else:
                    lg = a * lq + math.log1p(-math.exp(d))
                    sgn = -1
                llog = lfp - math.lgamma(a + 1) - math.lgamma(p - a + 1) - math.log(p)
                if llog - lg < -60:
                    continue
                if p <= 60:
                    G = gate(p, a, q)
                    if abs(G) <= 1:
                        continue
                    L = lyndon_count(p, a)
                    if L == 0:
                        continue
                    L_over_G = L / abs(G)
                else:
                    L_over_G = math.exp(llog - lg)
                cl = elig_class_fast(p, a, q)
                if p <= P_EXACT and (q, p, a) in exact:
                    ne = exact[(q, p, a)][0]
                    f = ne / lyndon_count(p, a)
                elif cl == "none":
                    f = 0.0
                elif cl == "all":
                    f = 1.0
                else:
                    trans.append((abs(a - p * l2 / lq), a, sgn, L_over_G))
                    row[sgn][0] += L_over_G
                    continue
                row[sgn][0] += L_over_G
                row[sgn][1] += f * L_over_G
                row[sgn][2] += f * L_over_G
            # transition window: Monte Carlo from the critical line outwards, per side; stop after two
            # consecutive clocks with no eligible sample (the eligible share falls off within ~10 clocks)
            for sgn_ in (1, -1):
                band = sorted(t for t in trans if t[2] == sgn_)
                zeros = 0
                for dist, a, _, LG in band:
                    if p <= P_MC and zeros < 2:
                        f = mc_elig_fraction(p, a, q, S=MC_SAMPLES if p <= 300 else 200)
                        trans_mc += 1
                        zeros = zeros + 1 if f == 0.0 else 0
                        row[sgn_][1] += f * LG
                        row[sgn_][2] += f * LG
                    else:
                        unsampled[sgn_] += LG
                        row[sgn_][2] += LG
            trans = []
            for s_ in (1, -1):
                cum[s_]["naive"] += row[s_][0]
                cum[s_]["ref_lo"] += row[s_][1]
                cum[s_]["ref_hi"] += row[s_][2]
            per_p.append((p, row))
            if p in (11, 20, 30, 40, 100, 300, 1000, 4000):
                marks[p] = {s_: dict(cum[s_]) for s_ in (1, -1)}
        results[q] = (marks, per_p)
        print()
        print(f"  q = {q}: cumulative expectations up to P (dyadic = positive cycles of {q}x+1; "
              f"{q}-adic = positive cycles of {q}x-1)   [{time.time() - t0:.1f}s, {trans_mc} MC clocks]")
        print(f"  {'P':>6} | {'naive dyadic':>12} {'refined dyadic':>15} | {'naive q-adic':>12} {'refined q-adic':>15}")
        for P in sorted(marks):
            m = marks[P]
            print(f"  {P:>6} | {m[1]['naive']:12.4f} {m[1]['ref_lo']:15.4f} | {m[-1]['naive']:12.4f} {m[-1]['ref_lo']:15.4f}")
        print(f"  refined = exact eligible counts (p <= {P_EXACT}), the proven classes 'none'/'all', and Monte Carlo on the")
        print(f"  transition band next to the critical line (p <= {P_MC}; each band stops after two consecutive clocks")
        print(f"  with no eligible sample).  Transition clocks outside the sampled band carry naive weight dyadic")
        print(f"  {unsampled[1]:.4f}, q-adic {unsampled[-1]:.4f}; they are counted as 0 (their Monte Carlo share is 0 at the")
        print(f"  band edge and the eligible share decreases away from the critical line).")
        nd = [(P, marks[P][1]["naive"]) for P in sorted(marks) if P >= 100]
        print("  naive dyadic sum minus ln P: " + ", ".join(f"P={P}: {v - math.log(P):+.3f}" for P, v in nd))
    # per-period contributions for small p (where the mass is)
    print()
    print("  per-period refined contributions (q = 3), p <= 30: dyadic / 3-adic (necklace level)")
    marks3, per3 = results[3]
    line = []
    for p, row in per3[:30]:
        line.append(f"{p}:{row[1][1]:.3f}/{row[-1][1]:.3f}")
    for k in range(0, 30, 10):
        print("   " + "  ".join(line[k:k + 10]))
    print()
    print("=" * 110)
    print("C2. Census (part A, p <= 40) against the expectations (Poisson model, sporadic necklaces)")
    print("=" * 110)
    obs = {(3, 1): 0, (3, -1): 1, (5, 1): 3, (5, -1): 0}
    for q in (3, 5):
        m = results[q][0][40]
        for s, nm in ((1, "dyadic"), (-1, f"{q}-adic")):
            lam_n = m[s]["naive"]
            lam_r = m[s]["ref_lo"]
            k = obs[(q, s)]

            def pois(lam, kk):
                return math.exp(-lam) * lam ** kk / math.factorial(kk)

            def cdf(lam, kk):  # P(X <= kk)
                return sum(pois(lam, j) for j in range(kk + 1))
            print(f"  q = {q} {nm:>8}: observed {k}; naive mean {lam_n:.3f}: P(X<={k}) = {cdf(lam_n, k):.3f}, "
                  f"P(X>={k}) = {1 - cdf(lam_n, k - 1) if k else 1.0:.3f};  refined mean {lam_r:.3f}: "
                  f"P(X<={k}) = {cdf(lam_r, k):.3f}, P(X>={k}) = {1 - cdf(lam_r, k - 1) if k else 1.0:.3f}")
    return results


def c3_bounds():
    print()
    print("=" * 110)
    print("C3. Bound hierarchy (log2 values).  N(p,a) = number of integral words on the clock.")
    print("  size   = log2((c_max - c_min)/|G| + 1)   [injectivity of w -> x_w + range of c_w]")
    print("  perig  = log2(p * (1/|2^(p/a) - 3|/2 + 1))  [distinct odd least points <= 1/|2^(p/a)-3|; Belaga 2003]")
    print("  sqrtC  = log2 sqrt(C)   [no second-moment (Parseval) bound can certify less than ~ sqrt(C - C^2/M)]")
    print("  main   = log2 C/|G|")
    print("=" * 110)
    print(f"  {'(p,a)':>14} {'log2|G|':>9} {'log2 C':>9} {'sqrtC':>9} {'size':>9} {'perig':>7} {'main':>9}")
    clocks = [(19, 12), (27, 17), (46, 29), (65, 41), (84, 53), (149, 94), (306 + 179, 306), (1054, 665),
              (100, 55), (100, 60), (1000, 550), (1000, 620), (24727, 15601), (50508, 31867)]
    for p, a in clocks:
        G = gate(p, a)
        M = abs(G)
        lC = (math.lgamma(p + 1) - math.lgamma(a + 1) - math.lgamma(p - a + 1)) / math.log(2)
        cmn, cmx = cmin_cmax(p, a)
        size = math.log2((cmx - cmn) // M + 1)
        r = 2 ** (p / a)
        perig = math.log2(p * (1.0 / abs(r - 3) / 2 + 1))
        main = lC - math.log2(M)
        print(f"  {str((p, a)):>14} {math.log2(M):9.2f} {lC:9.2f} {lC / 2:9.2f} {size:9.2f} {perig:7.2f} {main:9.2f}")
        assert perig <= size + 1e-9 or p < 30
        assert perig < lC / 2 or p < 30


def c4_archimedean():
    print()
    print("=" * 110)
    print("C4. Archimedean part of the spectrum: S(h)/C at fixed theta = a/p (normalised DP) and the limit series")
    print("  dyadic (theta < log_3 2): x_w ~ Y = sum_n 3^n 2^(-D_n), D_n = position of the (n+1)-th one from the end;")
    print("  3-adic (theta > log_3 2): |x_w| ~ Z = sum_n 2^(sigma_n)/3^(n+1), sigma_n = position of the n-th one;")
    print("  letters iid Bernoulli(theta) in the limit; Monte Carlo with 2*10^5 samples, 160 terms (s.e. ~ 0.002).")
    print("=" * 110)
    thetas = [0.30, 0.45, 0.55, 0.75, 0.90]
    ps = [50, 100, 200, 400, 800]
    for th in thetas:
        vals = []
        for p in ps:
            a = int(round(th * p))
            v = S_dp(p, a, [1, 2, 3], 3, normalize=True)
            vals.append(v)
        # Monte Carlo limit, in batches (memory)
        n_s, n_t, bs = 200000, 160, 20000
        Ys = []
        for _ in range(n_s // bs):
            Bm = RNG.random((bs, n_t)) < th
            rank = np.cumsum(Bm, axis=1) - 1
            if th < math.log(2, 3):
                pos = np.arange(1, n_t + 1)[None, :]  # distance from the end: 1..n_t
                terms = np.where(Bm, np.exp(rank * math.log(3) - pos * math.log(2)), 0.0)
            else:
                pos = np.arange(0, n_t)[None, :]  # position from the start: 0..n_t-1
                terms = np.where(Bm, np.exp(pos * math.log(2) - (rank + 1) * math.log(3)), 0.0)
            Ys.append(terms.sum(axis=1))
        Y = np.concatenate(Ys)
        mc = [complex(np.mean(np.exp(2j * math.pi * h * Y))) for h in (1, 2, 3)]
        print(f"  theta = {th:.2f}:")
        for p, v in zip(ps, vals):
            print(f"     p = {p:>4}:  S(1)/C = {v[0].real:+.4f}{v[0].imag:+.4f}i   S(2)/C = {v[1].real:+.4f}{v[1].imag:+.4f}i"
                  f"   S(3)/C = {v[2].real:+.4f}{v[2].imag:+.4f}i")
        print(f"     limit series (MC): E e(Y) = {mc[0].real:+.4f}{mc[0].imag:+.4f}i   E e(2Y) = {mc[1].real:+.4f}"
              f"{mc[1].imag:+.4f}i   E e(3Y) = {mc[2].real:+.4f}{mc[2].imag:+.4f}i")
        # fraction of words with x < 1 in the limit (dyadic) / near 1 (3-adic)
        fr = np.mod(Y, 1.0)
        if th < math.log(2, 3):
            print(f"     limit law: P(Y < 1) = {np.mean(Y < 1):.3f};  P(Y mod 1 < 0.05) = {np.mean(fr < 0.05):.4f}, "
                  f"P(Y mod 1 > 0.95) = {np.mean(fr > 0.95):.4f} (uniform: 0.05 each)")
        else:
            print(f"     limit law: P(Z < 1.05) = {np.mean(Y < 1.05):.4f} (Z >= 1 always); P(Z mod 1 < 0.05) = "
                  f"{np.mean(fr < 0.05):.4f}, P(Z mod 1 > 0.95) = {np.mean(fr > 0.95):.4f} (uniform: 0.05 each)")
        dev = abs(vals[-1][0] - mc[0])
        print(f"     |S(1)/C at p = 800  -  E e(Y)| = {dev:.4f}")
        assert dev < 0.05, (th, vals[-1][0], mc[0])
    print()
    print("  near-critical clocks (amp = max(2^p,3^a)/|G|): |S(h)|/C for h = 1, 2, 3 (normalised DP)")
    for p, a in [(19, 12), (27, 17), (46, 29), (65, 41), (84, 53), (149, 94), (485, 306), (1054, 665)]:
        v = np.abs(S_dp(p, a, [1, 2, 3], 3, normalize=True))
        ampl = max(2 ** p, 3 ** a) / abs(gate(p, a))
        lC = (math.lgamma(p + 1) - math.lgamma(a + 1) - math.lgamma(p - a + 1)) / math.log(2)
        extra = ""
        if p >= 140:  # certify the double-precision DP against a 40-digit DP where |S|/C is tiny
            from procgen_gates_20260925_core import S_dp_mpmath
            ref = abs(S_dp_mpmath(p, a, 1, 3, dps=40)) / math.comb(p, a)
            assert abs(ref - v[0]) <= 1e-6 * ref, (p, a, ref, v[0])
            extra = f"   [40-digit DP: {ref:.5e}]"
        print(f"     ({p},{a}) amp = {ampl:9.3g}:  |S(1)|/C = {v[0]:.3e}  |S(2)|/C = {v[1]:.3e}  |S(3)|/C = {v[2]:.3e}"
              f"   (1/sqrt(C) = 2^{-lC / 2:.1f}){extra}")


def main():
    lean_malloc()
    t = time.time()
    c1_c2()
    print(f"[C1-C2 {time.time() - t:.1f}s]", file=sys.stderr, flush=True)
    c3_bounds()
    c4_archimedean()
    report_mem("end")
    print(f"[C total {time.time() - t:.1f}s]", file=sys.stderr)


if __name__ == "__main__":
    main()
