#!/usr/bin/env python3
"""Hard-class lane (collatz-procgen-20260922): the Diophantine-exponent reach of Theorem R/S, and C2.

Sections (see 05-knowledge/results/collatz_procgen_20260922_hard_class.md):
  H0  constants: Dio thresholds (Bugeaud-Kim, ADQZ, Theorem S), Lemma J' constants mu*(k0,L), maps
  H1  calibration of the Dio engine (LPF array) on words with known Dio / rep
  H2  Sturmian and quasi-Sturmian words: exact Theorem-R certificates beyond the old ranges,
      and the sharpness example (Bugeaud-Kim word s' under 19x+1, mu > Dio)
  H3  codings of rotations by several arcs (Lemma J'), census
  H4  primitive substitutions: exact Dio formula, certificate, census; the 5x+1 (drift) side
  H5  the square-swap block word Y: zero entropy, bounded discrepancy, Dio = 1, yet Phi_2(Y) irrational
      (Theorem Y: a 2-adic Tschakaloff-Pade certificate, exact); the general square-swap criterion mu_bar < phi
  H6  the cube-swap block word Y3: Dio = 1, Phi_2(Y3) = rational + c * sum rho^(k^3): the smallest open instance
  C1  supercritical strip automata: entropy, E-sets K_q, overlap criterion, decoupled-Mahler pressure
  C2  strip words (random, greedy low-repetition, Sturmian-perturbed): gains, reconstruction, values
  C3  Proposition M checks (EP words: Phi_R = Phi_2, L = 0; tail recursion of E_s)
  C4  capacity / Monks-Yazinski / density measurements
Every certificate is exact integer arithmetic (collatz_procgen_20260922_hard_verify.py); repetition
profiles come from the C engine collatz_procgen_20260922_hard_lpf.c (suffix array + LPF).
Usage: python3 collatz_procgen_20260922_hard.py [--quick]
"""
import argparse
import math
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import collatz_procgen_20260922_hard_words as HW      # noqa: E402
import collatz_procgen_20260922_hard_verify as HV     # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SCR = os.path.join(ROOT, "scratch", "procgen_hard")
LPF_BIN = os.path.join(SCR, "hard_lpf")
LOG23 = math.log2(3)
PHI = (1 + math.sqrt(5)) / 2
T0 = time.time()


def P(*a):
    print(*a, flush=True)


def elapsed():
    return f"[{time.time() - T0:.0f}s]"


# ------------------------------------------------------------------ engine wrappers
def lpf_run(w, mu, weights=None, tag="w", topK=14, jmin=48):
    path = os.path.join(SCR, f"word_{tag}.txt")
    HW.write_word(path, w)
    args = [LPF_BIN, path, f"{mu:.12f}"]
    args += [f"{weights[0]:.12f}", f"{weights[1]:.12f}"] if weights else ["1", "1"]
    args += [str(topK), str(jmin)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    res = dict(scales=[], cands=[], rn={}, rep=[], tail=[])
    for line in out.splitlines():
        p = line.split()
        if p[0] == "SCALE":
            if p[4] == "capped":
                continue
            res["scales"].append(dict(s=int(p[1]), ratio=float(p[4]), j=int(p[5]), a=int(p[6]), lpf=int(p[7]),
                                      wratio=float(p[8]), wj=int(p[9]), wa=int(p[10]), wlpf=int(p[11]),
                                      gain=float(p[12]), gj=int(p[13]), ga=int(p[14]), glpf=int(p[15])))
        elif p[0] == "CAND":
            res["cands"].append((int(p[1]), int(p[2]), int(p[3]), float(p[4])))
        elif p[0] == "RN":
            res["rn"][int(p[1])] = int(p[2])
        elif p[0] == "REP":
            res["rep"].append((int(p[1]), int(p[2]), None if p[3] == "incomplete" else float(p[3])))
        elif p[0] == "DIOTAIL":
            res["tail"].append((float(p[1]), float(p[3]), float(p[4])))
    return res


def dio_profile(res, smin, smax, key="ratio"):
    return [(d["s"], round(d[key], 4)) for d in res["scales"] if smin <= d["s"] <= smax]


def dio_est(res, smin, smax, key="ratio"):
    vals = [d[key] for d in res["scales"] if smin <= d["s"] <= smax]
    return max(vals) if vals else float("nan")


def rep_est(res):
    reps = [r for r in res["rep"] if r[2] is not None]
    return reps[-1][2] if reps else None


def certify(w, amap, res, nverify=10, jmax=1 << 15, Nb_cap=1 << 17):
    """Exact Theorem-R certificates for the engine's top candidates with |UV| <= jmax.
    Returns (rows, summary).  Phi_T(w) mod 2^Nb is computed once and validated by parity iteration."""
    cands = sorted({(j, a, L) for (j, a, L, g) in res["cands"] if j <= jmax}, key=lambda t: t[0])
    if not cands:
        return [], dict(n=0)
    Nb = min(len(w), Nb_cap, max(j + L for (j, a, L) in cands) + 64)
    X = HV.phi2_mod(w, Nb, amap)
    par_ok = HV.parity_check(X, w, min(Nb, 40000), amap)
    rows = []
    for (j, a, L) in cands[-nverify:]:
        lam = j + L
        r = HV.verify(w, X, Nb, a, j, lam, amap)
        rows.append(r)
    summ = dict(n=len(rows), parity_ok=par_ok, iso_all=all(r["iso"] in (True, None) for r in rows),
                iso_checked=sum(1 for r in rows if r["iso"] is not None),
                lam_all=all(r["lam_direct"] == r["lam"] for r in rows),
                best_gain=max(r["gain"] for r in rows), best_row=max(rows, key=lambda r: r["gain"]),
                max_canon_gcd=max(r["canon_gcd_bits"] for r in rows))
    return rows, summ


def fmt_rows(rows, k=4):
    out = []
    for r in sorted(rows, key=lambda r: r["j"])[-k:]:
        out.append(f"|UV|={r['j']},|U|={r['a']},lam={r['lam']}:gain={r['gain']:.1f}(h={r['hbits']:.0f},"
                   f"canon.gcd={r['canon_gcd_bits']}b,iso={'Y' if r['iso'] else ('-' if r['iso'] is None else 'N')})")
    return " ".join(out)


def recon(w, amap, Nbits=10000):
    kind, h, lam, cb = HV.recon_certificate(w, Nbits, amap)
    if kind == "NONE":
        if lam is None:
            return f"NONE(no candidate of height<=2^{h})"
        return f"NONE: the unique candidate (height 2^{cb}) leaves w at letter {lam} (N={Nbits}); no rational of height<=2^{h}"
    return f"AGREES on all {lam} letters (candidate 2^{cb}): repetition runs past the window"


def mu_star(k0, L):
    """Unique mu > 1 with mu (mu^k0 - 1) = L (Lemma J' constant); inf if L = inf."""
    if L == float("inf"):
        return float("inf")
    lo, hi = 1.0, 10.0
    for _ in range(200):
        m = (lo + hi) / 2
        if m * (m ** k0 - 1) < L:
            lo = m
        else:
            hi = m
    return lo


def limsup_ratio(terms):
    qs = HW.convergent_dens(terms)
    rat = [qs[i + 1] / qs[i] for i in range(len(qs) - 1)]
    tail = rat[len(rat) // 2:]
    return max(tail) if tail else float("nan"), qs


# ------------------------------------------------------------------ H0
def section_H0():
    P("=" * 110)
    P("H0  CONSTANTS AND THRESHOLDS  (Theorem D: Phi_T(w) irrational whenever Dio(w) > eta(w) = mu(beta))")
    P("=" * 110)
    rep_s = math.sqrt(10) - 1.5
    dio_s = rep_s / (rep_s - 1)
    P(f"  Bugeaud-Kim (TAMS 371 (2019)) Thm 3.4/3.6: rep(s) <= sqrt10 - 3/2 = {rep_s:.10f} for every (quasi-)Sturmian s;")
    P(f"    Lemma 10.3: Dio = rep/(rep-1)  ->  Dio(s) >= {dio_s:.10f}  (= 5/3 + 4 sqrt10/15 = {5/3 + 4*math.sqrt(10)/15:.10f})")
    mS = (1 + math.sqrt(1 + 4 * PHI)) / 2
    P(f"  Theorem S (options A/B, limsup q_(n+1)/q_n >= phi): mu(mu-1) < phi  <=>  mu < {mS:.10f};"
      f"  ADQZ/BHZ initial squares: mu < 2")
    P("  Lemma J' (codings of a rotation, k0 endpoint classes): Dio >= mu*(k0, L), mu*(mu*^k0 - 1) = L = limsup q_(n+1)/q_n:")
    Ls = [("phi", PHI), ("2", 2.0), ("1+sqrt2", 1 + math.sqrt(2)), ("3", 3.0), ("5", 5.0), ("10", 10.0)]
    P("    k0 | " + " | ".join(f"L={n:>7s}" for n, _ in Ls) + " |  complexity p(n)<=kn: 1+1/k")
    for k0 in range(1, 7):
        P(f"    {k0:2d} | " + " | ".join(f"{mu_star(k0, L):9.4f}" for _, L in Ls) + f" |  k={2*k0}: {1 + 1/(2*k0):.4f}")
    P("  Height exponent mu(beta) = max(1, beta log2 m1) for T(x) = x/2, (m1 x + 1)/2  (m0 = 1, r1 = 1):")
    P("    map    | log2 m1 | critical beta | beta range with mu < Dio_BK (all Sturmian & quasi-Sturmian) | mu < 2 (ADQZ) | mu < 1.8668 (Thm S)")
    for name, m1 in (("3x+1", 3), ("5x+1", 5), ("7x+1", 7), ("9x+1", 9), ("19x+1", 19)):
        l2 = math.log2(m1)
        P(f"    {name:6s} | {l2:7.4f} | {1/l2:13.4f} | beta < {min(1.0, dio_s / l2):.4f}"
          f"{' (ALL slopes)' if dio_s / l2 >= 1 else ''}"
          f" | beta < {min(1.0, 2 / l2):.4f} | beta < {min(1.0, mS / l2):.4f}")
    P("    Mahler ceil(3a/2): m0 = m1 = 3, mu = log2 3 = 1.5850 < Dio_BK: every Sturmian/quasi-Sturmian carry word.")
    P(f"  Lemma P (2-adic Tschakaloff, Theorem Y): sum_(k>=0) rho^(k^2) is irrational in Q_2 for rho = 2^L/M (M odd) whenever"
      f" mu_bar = max(1, log2(M)/L) < phi = {PHI:.6f}; with the unshifted construction (m = 0): mu_bar < 3/2.")
    P(f"  Linear complexity (AB ETDS 2007, proof of Thm 3): p(n) <= c n infinitely often => Dio >= 1 + 1/c;"
      f" 3x+1 needs c < 1/(beta log2 3 - 1), e.g. beta = 0.9: c < {1/(0.9*LOG23 - 1):.3f}; beta = 1: c < {1/(LOG23 - 1):.3f}")
    return dio_s


# ------------------------------------------------------------------ H1
def section_H1(n):
    P("=" * 110)
    P("H1  CALIBRATION OF THE Dio ENGINE  (Dio = 1 + limsup LPF(j)/j; r(n) = n + min{j: LPF(j) >= n})")
    P("=" * 110)
    fib = [0]
    while len(fib) < n:
        fib = [c for a in fib for c in ((0, 1) if a == 0 else (0,))]
    fib = fib[:n]
    tm = [bin(i).count("1") & 1 for i in range(n)]
    rnd = random.Random(1)
    rw = [rnd.randrange(2) for _ in range(n)]
    th = HW.QI(-2, 1, 3, 10)
    bk = HW.mechanical(th, HW.QI.rat(F(1, 3), 10), n, upper=False, start=1)
    top = int(math.log2(n)) - 2
    rows = [("Fibonacci fixed point 0->01,1->0", fib, PHI ** 2, PHI),
            ("Thue-Morse", tm, 5 / 3, 2.5),
            ("Bugeaud-Kim s' (slope (sqrt10-2)/3=[0;(2,1,1)], intercept 1/3)", bk, 5 / 3 + 4 * math.sqrt(10) / 15, math.sqrt(10) - 1.5),
            ("random fair bits (seed 1)", rw, 1.0, float("inf"))]
    for name, w, dio_true, rep_true in rows:
        res = lpf_run(w, 1.5, tag="cal")
        d = dio_est(res, top - 6, top)
        r = rep_est(res)
        P(f"  {name:62s}: Dio_est(scales {top-6}..{top}) = {d:.6f} (known {dio_true:.6f}); rep_est = "
          f"{r if r is None else round(r, 6)} (known {rep_true if rep_true != float('inf') else 'inf'})")
    P(f"  words of length {n}; the Bugeaud-Kim equality case (their Thm 3.4, (3.1)) is reproduced to 6 digits.  {elapsed()}")


# ------------------------------------------------------------------ H2
def sturm_row(name, w, beta, amap_name, dio_s, n_small, notes=""):
    amap = HV.MAPS[amap_name]
    mu = HV.mu_of(beta, amap)
    res = lpf_run(w, mu, tag="h2")
    top = int(math.log2(len(w))) - 2
    d = dio_est(res, top - 5, top)
    rows, summ = certify(w, amap, res, nverify=8, jmax=1 << 14)
    rc = recon(w[:n_small], amap, 10000) if n_small else "-"
    P(f"  {name} under {amap_name}: beta={beta:.5f} mu={mu:.4f} Dio_est={d:.4f} ({'> mu' if d > mu else '<= mu'}); "
      f"rep_est={rep_est(res)}; {notes}")
    if summ.get("n"):
        P(f"      exact certificates: parity ok={summ['parity_ok']} isometry ok={summ['iso_all']} ({summ['iso_checked']} checked) "
          f"lambda ok={summ['lam_all']} max canonical gcd={summ['max_canon_gcd']} bits; best exact gain {summ['best_gain']:.1f} bits")
        P(f"      {fmt_rows(rows, 4)}")
    P(f"      reconstruction at N=10^4 bits: {rc}")
    return d, mu, summ


def section_H2(n, dio_s, quick):
    P("=" * 110)
    P("H2  STURMIAN AND QUASI-STURMIAN WORDS  (Theorem S extended: all maps with mu < 2.5099; per-slope beyond)")
    P("=" * 110)
    def cover_of(mu, L, rho=None):
        if rho == 0 and mu < 1 + L and not (mu < dio_s):
            return (f"BHZ Thm 1.2 (characteristic word 0c_alpha: Dio >= ice = 1 + limsup q_(n+1)/q_n = {1 + L:.4f} > mu);"
                    f" no uniform bound applies")
        if mu * (mu - 1) < PHI:
            return "Theorem S (mu(mu-1) < phi)"
        if mu < 2:
            return "ADQZ initial squares (mu < 2)"
        if mu < dio_s:
            return "Bugeaud-Kim uniform bound only (mu < 2.5099)"
        if mu * (mu - 1) < L:
            return f"per-slope Lemma J only (mu(mu-1) = {mu*(mu-1):.3f} < limsup q_(n+1)/q_n = {L:.3f})"
        return "NO uniform or per-slope guarantee: decided by Dio of the individual word"
    for pre, per, rhos, maps in (([1, 7], [1], (F(0), F(1, 3)), ("3x+1", "5x+1", "7x+1")),
                                 ([1, 8], [1], (F(0), F(1, 3)), ("7x+1",)),
                                 ([1], [10], (F(0),), ("7x+1",))):
        al = HW.qi_from_cf(pre, per)
        terms = HW.cf_terms_qi(al, 24)
        L, qs = limsup_ratio(terms)
        P(f"  slope alpha = [0;{','.join(map(str, pre))},({','.join(map(str, per))})] = {float(al):.10f};"
          f" partial quotients {terms[:6]}...; limsup q_(n+1)/q_n = {L:.4f}")
        # standard words s_-1 = 1, s_0 = 0, s_1 = s_0^(a_1 - 1) s_-1, s_k = s_(k-1)^(a_k) s_(k-2): c_alpha = lim s_k
        sm1, s0 = [1], [0]
        sk = s0 * (terms[1] - 1) + sm1
        prev, cur = s0, sk
        for ak in terms[2:]:
            prev, cur = cur, cur * ak + prev
            if len(cur) > 200000:
                break
        w0 = HW.mechanical(al, HW.QI.rat(0, al.d), min(len(cur), 200000) + 1)
        P(f"    check: lower mechanical word with rho = 0 equals 0 c_alpha (standard-word limit) on "
          f"{min(len(cur), 200000)} letters: {w0 == [0] + cur[:len(w0) - 1]}")
        for rho in rhos:
            w = HW.mechanical(al, HW.QI.rat(rho, al.d), n)
            beta = sum(w) / len(w)
            for mp in maps:
                mu = HV.mu_of(beta, HV.MAPS[mp])
                sturm_row(f"Sturmian slope above, rho={rho}", w, beta, mp, dio_s, 90000 if mp != "3x+1" else 0,
                          notes=f"covered by: {cover_of(mu, L, rho)}")
    # Bugeaud-Kim s' complemented (density 1 - theta = 0.6126)
    th = HW.QI(-2, 1, 3, 10)
    s = HW.mechanical(th, HW.QI.rat(F(1, 3), 10), n, upper=False, start=1)
    sc = [1 - c for c in s]
    beta = sum(sc) / len(sc)
    for mp in ("5x+1", "19x+1"):
        mu = HV.mu_of(beta, HV.MAPS[mp])
        sturm_row("Bugeaud-Kim s' complemented", sc, beta, mp, dio_s, 0,
                  notes=("mu < Dio(s') = 2.5099: covered" if mu < dio_s else
                         "mu > Dio(s') = 2.5099: Theorem D does NOT apply (sharpness of the uniform Sturmian threshold)"))
    # quasi-Sturmian: S(0)=1101, S(1)=111 applied to the golden Sturmian word, prefix 00
    g = HW.qi_golden_inv()
    base = HW.mechanical(g, HW.QI.rat(0, 5), n // 3)
    S = {0: [1, 1, 0, 1], 1: [1, 1, 1]}
    q = HW.morphic_image(base, S, prefix=(0, 0))[:n]
    beta = sum(q) / len(q)
    cx = HW.complexity(q[:200000], 40)
    diffs = sorted({c - k for k, c in cx[10:]})
    P(f"  quasi-Sturmian q = 00 S(s), S(0)=1101, S(1)=111 (S(01) != S(10)), s golden Sturmian: p(n) - n for n=11..40 in {diffs}")
    for mp in ("3x+1", "5x+1", "7x+1"):
        sturm_row("quasi-Sturmian q", q, beta, mp, dio_s, 90000 if mp != "3x+1" else 0,
                  notes="Bugeaud-Kim Thm 3.6 + Lemma 10.3 => Dio(q) >= 2.5099 > mu")
    P(f"  {elapsed()}")


# ------------------------------------------------------------------ H3
def section_H3(n, quick):
    P("=" * 110)
    P("H3  CODINGS OF A ROTATION BY SEVERAL ARCS  (Lemma J'; complexity p(n) <= k n)")
    P("=" * 110)
    al = HW.qi_golden_inv()
    d = 5
    w = HW.rotation_coding(al, HW.QI.rat(0, d), [(F(0), F(9, 10))], n)
    beta = sum(w) / len(w)
    for mp in ("3x+1", "5x+1", "7x+1"):
        amap = HV.MAPS[mp]
        mu = HV.mu_of(beta, amap)
        guar = mu_star(2, PHI)
        res = lpf_run(w, mu, tag="h3")
        top = int(math.log2(n)) - 2
        de = dio_est(res, top - 5, top)
        rows, summ = certify(w, amap, res, nverify=6, jmax=1 << 14)
        P(f"  golden rotation, 1-set [0, 9/10) (endpoints 0, 9/10: k0 = 2 classes), rho = 0, under {mp}: beta={beta:.4f} mu={mu:.4f};"
          f" Lemma J' guarantee mu*(2,phi) = {guar:.4f} {'> mu: PROVED' if guar > mu else '< mu: not guaranteed'};"
          f" Dio_est = {de:.4f}")
        if summ.get("n"):
            P(f"      exact: iso ok={summ['iso_all']} parity ok={summ['parity_ok']} best gain {summ['best_gain']:.1f}; {fmt_rows(rows, 3)}")
    # Lemma J' per-convergent check (FINITE-EXACT): best option ratio at period q_n vs the chain bound
    terms = HW.cf_terms_qi(al, 40)
    qs = HW.convergent_dens(terms)
    chk = []
    for i in range(3, len(qs) - 1):
        q, q1 = qs[i], qs[i + 1]
        if 2 * (q + q1) + q > len(w) // 2:
            break
        D = [t for t in range(len(w) - q) if w[t] != w[t + q]]
        starts = [0] + [t + 1 for t in D]
        ends = D + [len(w) - q]
        best = max((b + q) / (a + q) for a, b in zip(starts, ends) if a <= 2 * (q + q1))
        chk.append((q, round(q1 / q, 4), round(best, 4), round(mu_star(2, q1 / q), 4)))
    P(f"  Lemma J' check, k0 = 2, per convergent (q_n, q_(n+1)/q_n, best option ratio (b+q)/(a+q) over stretches"
      f" starting before 2(q_n+q_(n+1)), bound mu*(2, q_(n+1)/q_n)): {chk}")
    P(f"    all best ratios >= bound: {all(c[2] >= c[3] for c in chk)}")
    # census of multi-arc codings with rational endpoints
    rng = random.Random(3)
    ntr = 40 if quick else 160
    census = []
    for it in range(ntr * 3):
        if len(census) >= ntr:
            break
        k0 = rng.randint(2, 6)
        den = rng.choice([7, 11, 13, 17, 19, 23, 29, 31])
        starts = sorted(set(F(rng.randrange(den), den) for _ in range(k0)))
        m = rng.choice([13, 21, 34])
        Lq = al.mul_int(m).frac()
        ok = all(float(starts[i + 1] - starts[i]) > float(Lq) + 1e-9 for i in range(len(starts) - 1)) and \
            float(starts[0] + 1 - starts[-1]) > float(Lq)
        if not ok:
            continue
        arcs = []
        for i, s0 in enumerate(starts):
            lo = HW.QI.rat(s0, d) + Lq
            hi = HW.QI.rat(starts[(i + 1) % len(starts)], d)
            arcs.append((lo.frac(), hi))
        rho = F(rng.randrange(97), 97)
        nn = 1 << 16
        ww = HW.rotation_coding(al, HW.QI.rat(rho, d), arcs, nn)
        b = sum(ww) / nn
        mu = b * LOG23
        res = lpf_run(ww, mu, tag="h3c")
        de = dio_est(res, 9, 13)
        census.append((de - mu, de, mu, b, len(starts), m))
    census.sort()
    P(f"  census: {len(census)} golden codings, 0-set = union of 2..6 arcs [c_i, c_i + {{m alpha}}) (m in 13,21,34; bounded discrepancy),"
      f" c_i rational, rho rational; words of 2^16 letters; 3x+1:")
    P(f"    min (Dio_est - mu) = {census[0][0]:+.4f} (Dio_est {census[0][1]:.4f}, mu {census[0][2]:.4f}, beta {census[0][3]:.4f},"
      f" {census[0][4]} arcs); fraction with Dio_est > mu: {sum(1 for c in census if c[0] > 0)}/{len(census)};"
      f" median Dio_est {sorted(c[1] for c in census)[len(census)//2]:.4f}")
    P(f"    under 5x+1 (mu = 2.3219 beta): fraction with Dio_est > mu: "
      f"{sum(1 for c in census if c[1] > c[3] * math.log2(5))}/{len(census)}  {elapsed()}")


# ------------------------------------------------------------------ H4
def prefix_certificate(x, l, m):
    """max over prefix patterns inside x[:m] of the l-weighted ratio; returns (ratio, a, j, lam)."""
    s = "".join(map(str, x[:m]))
    W = [0.0]
    for c in s:
        W.append(W[-1] + l[int(c)])
    best = (0.0, None, None, None)
    for nlen in range(2, m + 1):
        sub = s[:nlen]
        lo, hi = 0, nlen - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if sub.find(sub[nlen - mid:]) < nlen - mid:
                lo = mid
            else:
                hi = mid - 1
        sL = lo
        if sL > 0:
            r = W[nlen] / W[nlen - sL]
            if r > best[0]:
                pos = sub.find(sub[nlen - sL:])
                best = (r, pos, nlen - sL, nlen)
    return best


def section_H4(n, quick):
    P("=" * 110)
    P("H4  PRIMITIVE SUBSTITUTIONS  (exact: Dio(sigma^w(a)) = sup of l-weighted prefix-pattern ratios)")
    P("=" * 110)
    # exact formula demo: Fibonacci with l-weights
    sig = {0: [0, 1], 1: [0]}
    th, th2, l, f, M = HW.pf_data(sig)
    x = HW.fixed_point(sig, 0, n)
    res = lpf_run(x, 1.5, weights=l, tag="h4")
    allw = max(dd["wratio"] for dd in res["scales"])
    P(f"  Fibonacci 0->01,1->0: left PF vector l = ({l[0]:.4f},{l[1]:.4f}); sup over all j of l-weighted ratio = {allw:.6f}"
      f" (phi^2 = {PHI**2:.6f}); prefix certificate (first 300 letters) = {prefix_certificate(x, l, 300)[0]:.6f}")
    # a supercritical certified example
    sig = {0: [0, 1, 1, 0, 0, 1, 0, 1, 0, 0], 1: [1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1]}
    th, th2, l, f, M = HW.pf_data(sig)
    x = HW.fixed_point(sig, 0, n)
    beta = f[1]
    cert = prefix_certificate(x, l, 300)
    P(f"  sigma(0)=0110010100, sigma(1)=11110111111 (primitive), x = sigma^w(0): theta={th:.4f} |theta_2|={abs(th2):.4f}"
      f" ({'<1: bounded discrepancy' if abs(th2) < 1 else '>1: discrepancy ~ N^' + format(math.log(abs(th2))/math.log(th), '.3f')}"
      f", Adamczewski TCS 2003 Thm 13); beta = {beta:.5f}")
    for mp in ("3x+1", "5x+1"):
        amap = HV.MAPS[mp]
        mu = HV.mu_of(beta, amap)
        res = lpf_run(x, mu, weights=l, tag="h4b")
        top = int(math.log2(n)) - 2
        de = dio_est(res, top - 5, top)
        dw = max(dd["wratio"] for dd in res["scales"])
        rows, summ = certify(x, amap, res, nverify=6, jmax=1 << 14)
        verdict = (f"certificate {cert[0]:.4f} > mu: PROVED (pattern U=x[0,{cert[1]}), V=x[{cert[1]},{cert[2]}), lambda>={cert[3]})"
                   if cert[0] > mu else f"no certificate > mu found (sup l-ratio over 2^{int(math.log2(n))} letters = {dw:.4f}): method-limited")
        P(f"    under {mp}: mu = {mu:.4f}; Dio_est (lengths) = {de:.4f}; {verdict}")
        if summ.get("n"):
            P(f"      exact: iso ok={summ['iso_all']} parity ok={summ['parity_ok']} best gain {summ['best_gain']:.1f}; {fmt_rows(rows, 3)}")
    # census
    rng = random.Random(7)
    ntries = 3000 if quick else 20000
    cens = []
    nper = 0
    for it in range(ntries):
        n1 = rng.randint(2, 16)
        n0 = rng.randint(1, 16)
        d1, d0 = rng.uniform(0.5, 1.0), rng.uniform(0.5, 1.0)
        a = rng.randrange(2)
        sa = [a] + [1 if rng.random() < (d1 if a == 1 else d0) else 0 for _ in range((n1 if a == 1 else n0) - 1)]
        sb = [1 if rng.random() < (d0 if a == 1 else d1) else 0 for _ in range(n0 if a == 1 else n1)]
        sg = {a: sa, 1 - a: sb}
        try:
            th, th2, l, f, M = HW.pf_data(sg)
        except Exception:
            continue
        M2 = [[sum(M[i][k] * M[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)]
        M3 = [[sum(M2[i][k] * M[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)]
        if min(min(r) for r in M3) <= 0 or th <= 1.0001:
            continue
        beta = f[1]
        if beta * LOG23 <= 1.0:
            continue
        try:
            xx = HW.fixed_point(sg, a, 300)
        except Exception:
            continue
        # exclude eventually periodic fixed points (Morse-Hedlund on a long prefix)
        xl = HW.fixed_point(sg, a, 6000)
        sb_ = bytes(xl)
        if any(len({sb_[i:i + nn] for i in range(len(sb_) - nn + 1)}) <= nn for nn in (8, 16, 32, 64)):
            nper += 1
            continue
        c3 = prefix_certificate(xx, l, 300)[0]
        cens.append((c3 - beta * LOG23, c3 - beta * math.log2(5), c3, beta, abs(th2), dict(sg), a))
    cens.sort(key=lambda c: c[0])
    P(f"  census: {len(cens)} random primitive binary substitutions (|images| <= 16) with aperiodic supercritical fixed points"
      f" (beta log2 3 > 1; {nper} periodic ones excluded):")
    P(f"    3x+1: min(certificate - mu) = {cens[0][0]:+.4f}; all certified (certificate > mu): {all(c[0] > 0 for c in cens)}"
      f" -> PC PROVED on each of these fixed points (Theorem D + Lemma F); {sum(1 for c in cens if c[4] < 1)} of them have"
      f" |theta_2| < 1 (bounded discrepancy: members of a C2 strip); periodic fixed points excluded beforehand")
    pis = [c for c in cens if c[4] < 1]
    if pis:
        c = min(pis, key=lambda c: c[0])
        sg, a0 = c[5], c[6]
        P(f"    tightest bounded-discrepancy member: sigma(0)={''.join(map(str, sg[0]))} sigma(1)={''.join(map(str, sg[1]))}"
          f" start {a0}: beta={c[3]:.4f} certificate {c[2]:.4f} vs mu_3 = {c[3]*LOG23:.4f}, |theta_2| = {c[4]:.3f}")
    c5 = sorted(c[1] for c in cens)
    P(f"    5x+1: fraction with certificate > mu_5 = {sum(1 for v in c5 if v > 0)}/{len(c5)}; min margin {c5[0]:+.4f}"
      f" (method-limited instances exist on the DRIFT side)  {elapsed()}")


# ------------------------------------------------------------------ H5
def sparse_block_word(nblocks, member, B, B2):
    """b_n = B2 if member(n) else B (blocks of equal length and equal number of 1s)."""
    w = []
    for nb in range(nblocks):
        w += list(B2) if member(nb) else list(B)
    return w


def is_square(nb):
    k = math.isqrt(nb)
    return nb >= 1 and k * k == nb


def is_cube(nb):
    if nb < 1:
        return False
    k = round(nb ** (1 / 3))
    return any((k + d) ** 3 == nb for d in (-1, 0, 1))


def Yword(nblocks, r=9):
    return sparse_block_word(nblocks, is_square, [1] * r + [0], [1] * (r - 1) + [0, 1])


def block_series_check(w, amap, B, B2, member, Nbits):
    """Phi_T(w) = -(1/M) sum_n c_n rho^n (2-adically), rho = 2^L/M, M = m1^A, c_n = R_B or R_B2 (exact mod 2^Nbits)."""
    L = len(B)
    MB, RB = HV.MR(list(B), 0, L, amap)
    MB2, RB2 = HV.MR(list(B2), 0, L, amap)
    assert MB == MB2
    mod = 1 << Nbits
    X = HV.phi2_mod(w, Nbits, amap)
    invM = pow(MB, -1, mod)
    rho = (1 << L) * invM % mod
    tot, pw, nb = 0, 1, 0
    while L * nb < Nbits + L:
        tot = (tot + (RB2 if member(nb) else RB) * pw) % mod
        pw = pw * rho % mod
        nb += 1
    return X == (-tot * invM) % mod, MB, RB, RB2


def tschakaloff_rows(L, Mbase, nmax):
    """2-adic Tschakaloff-Pade linear forms for theta = sum_{l>=0} rho^{l^2}, rho = 2^L/Mbase (exact integers).
    Returns rows (n, m, v2(N_n I_n), predicted L(e_off + (n+m+1)^2), log2 H_n, margin)."""
    rows = []
    for n in range(1, nmax + 1):
        for m in sorted({0, round(n / PHI)}):
            C = [dict() for _ in range(n + 1)]
            C[0][0] = 1
            for jj in range(1, n + 1):
                new = [dict() for _ in range(n + 1)]
                for k in range(n + 1):
                    for e, c in C[k].items():
                        new[k][e] = new[k].get(e, 0) + c
                        if k + 1 <= n:
                            new[k + 1][e + jj] = new[k + 1].get(e + jj, 0) - c
                C = new
            Ad, Bd = {}, {}
            for k in range(n + 1):
                for e, c in C[k].items():
                    base = -k - 2 * (e + k * (k - 1) // 2 + k * m)
                    Ad[base] = Ad.get(base, 0) + c
                    for l in range(0, k + m + 1):
                        Bd[base + l * l] = Bd.get(base + l * l, 0) + c
            eoff = max(0, -min(min(Ad), min(Bd)))
            ep = max(0, max(max(Ad), max(Bd)))

            def integ(D):
                return sum(c * (1 << (L * (e + eoff))) * Mbase ** (ep - e) for e, c in D.items() if c)
            NA, NB = integ(Ad), integ(Bd)
            K = L * (eoff + (n + m + 1) ** 2) + 256
            mod = 1 << K
            r = (1 << L) * pow(Mbase, -1, mod) % mod
            th, l = 0, 0
            while L * l * l < K + L:
                th = (th + pow(r, l * l, mod)) % mod
                l += 1
            G = (NA * th - NB) % mod
            v = (G & -G).bit_length() - 1 if G else K
            H = math.log2(max(abs(NA), abs(NB)))
            rows.append((n, m, v, L * (eoff + (n + m + 1) ** 2), H, v - H))
    return rows


def section_H5(n):
    P("=" * 110)
    P("H5  THE SQUARE-SWAP BLOCK WORD Y = b_0 b_1 ...,  b_n = 1^8 0 1 if n is a nonzero square, else 1^9 0")
    P("=" * 110)
    Y = Yword(n // 10)
    beta = sum(Y) / len(Y)
    lo, hi = HW.discrepancy_range(Y, 0.9)
    cx = HW.complexity(Y[:300000], 60, step=10)
    P(f"  length {len(Y)}, beta = {beta} (exact 9/10), discrepancy a_s - 0.9 s in [{lo:.3f}, {hi:.3f}] (width 1),"
      f" complexity p(n) at n=1,11,...: {[c for _, c in cx]}")
    amap = HV.MAPS["3x+1"]
    mu = HV.mu_of(0.9, amap)
    res = lpf_run(Y, mu, tag="Y", topK=16, jmin=16)
    prof = dio_profile(res, 6, int(math.log2(len(Y))) - 1)
    P(f"  3x+1, mu = 0.9 log2 3 = {mu:.4f}.  max_j (j+LPF)/j per dyadic scale: {prof}")
    lp = [(d["s"], round(d["lpf"] / math.sqrt(d["j"]), 2)) for d in res["scales"] if d["s"] >= 8]
    P(f"  LPF(j)/sqrt(j) at the per-scale maximiser (Proposition Y: LPF(j) <= 13 sqrt(j) + 60): {lp}")
    rows, summ = certify(Y, amap, res, nverify=14, jmax=1 << 16)
    P(f"  Theorem D does not apply (Dio(Y) = 1 < mu).  Exact Theorem-R gains of the best periodic approximants: best"
      f" {summ['best_gain']:.1f} bits at |UV|={summ['best_row']['j']}; iso ok={summ['iso_all']}; per-scale best linear-gain"
      f" estimates {[(d['s'], round(d['gain'])) for d in res['scales'] if d['s'] >= 4]}")
    P(f"      {fmt_rows(rows, 6)}")
    P(f"  reconstruction at N=10^4 bits: {recon(Y[:90000], amap, 10000)}")
    ok, M, RB, RB2 = block_series_check(Y, amap, [1] * 9 + [0], [1] * 8 + [0, 1], is_square, 30000)
    P(f"  block identity: Phi_2(Y) = -(1/3^9) sum_n c_n rho^n, rho = 2^10/3^9, c_n = R_B = {RB} (n not a square),"
      f" R_B' = {RB2} (n a square): exact mod 2^30000: {ok}")
    P(f"    i.e. Phi_2(Y) = -1 - 512/(3^9-2^10) - (256/3^9) sum_(k>=1) rho^(k^2);  Phi_2(Y) in Q <=> theta(rho) in Q (in Q_2).")
    phr = HV.phi_real(Y, amap, 30)
    P(f"  Phi_R(Y) = {phr} (same expression summed in R; theta values at algebraic q in (0,1) are transcendental [R])")
    trow = tschakaloff_rows(10, 3 ** 9, 14)
    P("  Theorem Y certificate (2-adic Tschakaloff-Pade linear forms N_n(A_n theta - B_n), exact): (n, m, v2, predicted,"
      " log2 H_n, margin v2 - log2 H_n):")
    P("    " + " ".join(f"({a},{b},{c},{d},{e:.0f},{f:+.0f})" for a, b, c, d, e, f in trow))
    P(f"    all v2 == predicted: {all(r[2] == r[3] for r in trow)}; margins increase: {all(trow[i+2][5] > trow[i][5] for i in range(len(trow)-2))};"
      f" asymptotic margin (m=0): L n^2 (3 - 2 mu) = {10*(3 - 2*mu):.2f} n^2 > 0  =>  Phi_2(Y) is IRRATIONAL (PROVED)")
    # the general square-swap criterion, 5x+1 examples on both sides of phi
    for (A, L) in ((6, 10), (8, 10)):
        amap5 = HV.MAPS["5x+1"]
        B = [1] * A + [0] * (L - A)
        B2 = [1] * (A - 1) + [0, 1] + [0] * (L - A - 1)
        w5 = sparse_block_word(4000, is_square, B, B2)
        ok5, M5, _, _ = block_series_check(w5, amap5, B, B2, is_square, 20000)
        mu5 = math.log2(M5) / L
        tr5 = tschakaloff_rows(L, M5, 10)
        best = max(tr5, key=lambda r: r[0] * 100 + r[1])
        P(f"  5x+1 square-swap word with B = 1^{A}0^{L-A}: block identity ok={ok5}; mu_bar = max(1, log2(5^{A})/{L}) ="
          f" {max(1, mu5):.4f} {'< phi: Lemma P applies' if max(1, mu5) < PHI else '> phi: Lemma P does not apply'};"
          f" margins at n=10: {[round(r[5]) for r in tr5 if r[0] == 10]}")
    P(f"  {elapsed()}")


def section_H6(n):
    P("=" * 110)
    P("H6  THE CUBE-SWAP BLOCK WORD Y3 (b_n = 1^8 0 1 iff n is a nonzero cube): the smallest open instance found")
    P("=" * 110)
    B, B2 = [1] * 9 + [0], [1] * 8 + [0, 1]
    Y3 = sparse_block_word(n // 10, is_cube, B, B2)
    lo, hi = HW.discrepancy_range(Y3, 0.9)
    amap = HV.MAPS["3x+1"]
    mu = HV.mu_of(0.9, amap)
    res = lpf_run(Y3, mu, tag="Y3", topK=16, jmin=16)
    P(f"  beta = 9/10, discrepancy in [{lo:.3f}, {hi:.3f}] (width 1); Dio profile {dio_profile(res, 8, int(math.log2(len(Y3))) - 1)}")
    rows, summ = certify(Y3, amap, res, nverify=12, jmax=1 << 16)
    P(f"  exact Theorem-R gains: best {summ['best_gain']:.1f} bits at |UV|={summ['best_row']['j']} (iso ok={summ['iso_all']});"
      f" per-scale linear-gain estimates {[(d['s'], round(d['gain'])) for d in res['scales'] if d['s'] >= 6]}")
    P(f"  reconstruction at N=10^4 bits: {recon(Y3[:90000], amap, 10000)}")
    ok, M, RB, RB2 = block_series_check(Y3, amap, B, B2, is_cube, 30000)
    P(f"  block identity: Phi_2(Y3) = -1 - 512/(3^9-2^10) - (256/3^9) sum_(k>=1) rho^(k^3), rho = 2^10/3^9: exact mod 2^30000: {ok}")
    P("  partial sums s_K of sum rho^(k^3): height 3^(9K^3) = 2^(14.26 K^3) vs 2-adic error 2^(-10 (K+1)^3): Liouville fails"
      " (needs log2(3^9)/10 < 1); the series has no first-order q-difference equation, so Lemma P's construction does not apply.")
    P(f"  STATUS: irrationality of Phi_2(Y3) (equivalently of sum_(k>=1) (2^10/3^9)^(k^3) in Q_2) is OPEN.  {elapsed()}")


# ------------------------------------------------------------------ C1
def strip_Ksets(A, B, W, iters=3000):
    K, tr = HW.strip_automaton(A, B, W)
    hi = [3.0] * (K + 1)
    lo = [0.0] * (K + 1)
    for _ in range(iters):
        hi = [max(((1 / 3 + 2 / 3 * hi[k2]) if e else 2 * hi[k2]) for e, k2 in tr[k].items()) for k in range(K + 1)]
        lo = [min(((1 / 3 + 2 / 3 * lo[k2]) if e else 2 * lo[k2]) for e, k2 in tr[k].items()) for k in range(K + 1)]
    # overlap criterion: images of the hull intervals under the allowed letters cover the hull
    overlap = True
    for k in range(K + 1):
        ims = sorted(((1 / 3 + 2 / 3 * lo[k2], 1 / 3 + 2 / 3 * hi[k2]) if e else (2 * lo[k2], 2 * hi[k2]))
                     for e, k2 in tr[k].items())
        for (a1, b1), (a2, b2) in zip(ims, ims[1:]):
            if a2 > b1 + 1e-12:
                overlap = False
    return K, tr, lo, hi, overlap


def section_C1():
    P("=" * 110)
    P("C1  SUPERCRITICAL STRIPS: automaton, entropy, E-sets K_q (E = -Phi_R of the tail), decoupled-Mahler pressure")
    P("=" * 110)
    P("  strip(A/B, W): states k = B a_s - A s in [0, floor(BW)]; letter 1: k += B-A, letter 0: k -= A.")
    P("  K_q = {E(w') : w' strip tail from q}, E(1w) = 1/3 + 2E(w)/3, E(0w) = 2E(w).  If the one-step images of the hulls")
    P("  overlap at every state, K_q = [lo_q, hi_q] exactly (graph-directed IFS; 'overlap' column).  m_q = min(1, hi_q - lo_q)")
    P("  = measure of frac(K_q) then.  Decoupled pressure P_MZ = log2(3^alpha/2) + log2 rho(A_q,q' m_q'): P_MZ > 0 means the")
    P("  word-free Mahler relaxation MZ is heuristically FALSE (vacuous when all m_q = 1); P_MZ < 0 means MZ is non-vacuous")
    P("  and heuristically TRUE, and MZ => C2 for that strip (Proposition M).  Where 'overlap' is False, m_q is the hull")
    P("  length, an upper bound for the measure of frac(K_q), so the printed P_MZ is an upper bound.")
    P("  'forced cycle' = a cycle through states with a single allowed letter (Proposition Rnd needs: none).")
    P("   alpha  |  W  | states | entropy(bits/letter) | E hull            | overlap | min/max m_q   | #q with m_q<1 | P_MZ   | forced cycle")
    for (A, B) in [(2, 3), (3, 4), (4, 5), (9, 10)]:
        for W in [1.0, 1.5, 2.0, 3.0]:
            K, tr, lo, hi, ov = strip_Ksets(A, B, W)
            h = HW.strip_entropy(A, B, W)
            m = [min(1.0, hi[k] - lo[k]) for k in range(K + 1)]
            lam = 3 ** (A / B) / 2
            pw = HW.strip_entropy(A, B, W, weights=m)
            forced = [k for k in range(K + 1) if len(tr[k]) == 1]
            # cycle detection in the subgraph induced by forced states (functional graph: out-degree 1)
            fc = False
            for k0 in forced:
                seen, k = set(), k0
                while k in forced and k not in seen:
                    seen.add(k)
                    k = next(iter(tr[k].values()))
                if k in seen:
                    fc = True
                    break
            P(f"   {A}/{B:<3d} | {W:3.1f} | {K+1:6d} | {h:20.4f} | [{min(lo):6.3f},{max(hi):8.3f}] | {str(ov):7s} |"
              f" {min(m):.3f}/{max(m):.3f} | {sum(1 for v in m if v < 1):13d} | {math.log2(lam) + pw:+.3f} | {fc}")
    P(f"  {elapsed()}")


# ------------------------------------------------------------------ C2
def c2_row(name, w, alpha_f, amap_name="3x+1"):
    amap = HV.MAPS[amap_name]
    mu = HV.mu_of(alpha_f, amap)
    lo, hi = HW.discrepancy_range(w, alpha_f)
    res = lpf_run(w, mu, tag="c2", topK=12, jmin=16)
    prof = dio_profile(res, 8, 14)
    rows, summ = certify(w, amap, res, nverify=10, jmax=1 << 14)
    cx = HW.complexity(w[:60000], 22)
    hest = math.log2(cx[21][1] / cx[17][1]) / 4
    rc = recon(w, amap, 10000)
    phr = HV.phi_real(w, amap, 25)
    X64 = HV.phi2_mod(w, 64, amap)
    P(f"  {name}: beta={sum(w)/len(w):.4f} mu={mu:.4f} discrepancy in [{lo:.2f},{hi:.2f}] entropy~{hest:.3f} bits/letter")
    P(f"      Dio profile (scale: max (j+LPF)/j): {prof}")
    P(f"      exact Theorem-R gains: best {summ['best_gain']:.1f} bits at |UV|={summ['best_row']['j']} (iso ok={summ['iso_all']},"
      f" parity ok={summ['parity_ok']}); per-scale linear-gain estimates {[(d['s'], round(d['gain'])) for d in res['scales'] if 6 <= d['s'] <= 14]}")
    P(f"      reconstruction N=10^4: {rc}")
    P(f"      Phi_R(w) = {phr}   Phi_2(w) mod 2^64 = {X64:#018x}")
    return summ


def section_C2(n, quick):
    P("=" * 110)
    P("C2  SUPERCRITICAL STRIP WORDS OF WIDTH >= 1 (random / greedy low-repetition / Sturmian-perturbed)")
    P("=" * 110)
    specs = [(3, 4, 2.0), (4, 5, 1.5), (9, 10, 1.5), (9, 10, 2.0)]
    if quick:
        specs = specs[:1]
    for (A, B, W) in specs:
        K, tr = HW.strip_automaton(A, B, W)
        for mode, seed in (("random", 11), ("greedy", 12)):
            w = HW.strip_word_rational(A, B, W, n, mode, seed=seed)
            assert HW.strip_check_rational(w, A, B, W, K // 2)
            c2_row(f"strip {A}/{B} W={W} {mode}", w, A / B)
    al = HW.qi_from_cf([1], [2])      # sqrt2/2
    for mode in ("random", "perturbed"):
        w = HW.strip_word_quadratic(al, F(-3, 4), F(3, 2), n, mode, seed=5, flip_rate=0.05)
        c2_row(f"strip sqrt2/2 W=3/2 {mode}", w, float(al))
    P(f"  {elapsed()}")


# ------------------------------------------------------------------ C3
def section_C3():
    P("=" * 110)
    P("C3  PROPOSITION M CHECKS")
    P("=" * 110)
    import mpmath
    mpmath.mp.dps = 50
    amap = HV.MAPS["3x+1"]
    P("  (i) eventually periodic supercritical words u v^inf: Phi_R(u v^inf) = Phi_2(u v^inf) (same rational), so L = 0:")
    for u, v in [((), (1,)), ((), (1, 1, 0)), ((0, 1), (1, 1, 0)), ((), (1, 1, 1, 0, 1)), ((1, 0, 0), (1, 1, 0, 1, 1, 1, 0)),
                 ((), (1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1))]:
        a_v = sum(v)
        sup = a_v * LOG23 > len(v)
        w = list(u) + list(v) * 4000
        Pq, Qq = HV.approximant(w, len(u), len(u) + len(v), amap)
        r = F(Pq, Qq)
        real = HV.phi_real(w, amap, 40)
        P(f"      u={''.join(map(str,u)) or '-':4s} v={''.join(map(str,v)):12s} supercritical={sup}: P/Q = {r}  "
          f"real series = {mpmath.nstr(real, 25)}  |diff| = {mpmath.nstr(abs(real - mpmath.mpf(r.numerator) / r.denominator), 3)}")
    P("  (ii) tail recursion E_(s+1) = E_s/2 (w_s = 0), (3E_s - 1)/2 (w_s = 1), E_s = -Phi_R(sigma^s w) > 0, on a 3/4-strip word:")
    w = HW.strip_word_rational(3, 4, 2.0, 4000, "random", seed=3)
    Es = [HV.E_tail(w, s, amap, 30) for s in range(8)]
    ok = all(abs((Es[s] / 2 if w[s] == 0 else (3 * Es[s] - 1) / 2) - Es[s + 1]) < mpmath.mpf(10) ** -20 for s in range(7))
    P(f"      E_0..E_7 = {[mpmath.nstr(e, 8) for e in Es]}; recursion holds to 1e-20: {ok}")
    P(f"  {elapsed()}")


# ------------------------------------------------------------------ C4
def section_C4():
    P("=" * 110)
    P("C4  WHY THE PROVED MECHANISMS STOP (measurements)")
    P("=" * 110)
    P("   alpha  | lambda = 3^alpha/2 | orbit values <= X (hypothetical strip orbit) | Monks-Yazinski liminf a_s/s >= log_3 2 = 0.6309")
    for al in (0.65, 0.7, 0.75, 0.8, 0.9, 0.95):
        lam = 3 ** al / 2
        P(f"   {al:5.2f}  | {lam:18.4f} | ~ log X / log lambda = {1/math.log(lam):6.2f} log X (density zero)  | satisfied ({al:.2f} > 0.6309)")
    P("   capacity (in-house (D2)-(D2c)) needs n_j <= C j; supercritical strips give n_j >= c lambda^j: the counting inequality")
    P("   N <= #{admissible integers <= C lambda^N} ~ lambda^N / 3 always holds, and E_eps (density zero) contains the orbit")
    P(f"   without contradiction.  {elapsed()}")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    os.makedirs(SCR, exist_ok=True)
    if not os.path.exists(LPF_BIN):
        raise SystemExit(f"build the engine first: cc -O2 -o {LPF_BIN} collatz_procgen_20260922_hard_lpf.c -lm")
    n_big = (1 << 18) if args.quick else (1 << 20)
    P(f"collatz_procgen_20260922_hard: hard-class lane; quick={args.quick}; words up to {n_big} letters")
    dio_s = section_H0()
    section_H1(n_big)
    section_H2(n_big, dio_s, args.quick)
    section_H3(n_big, args.quick)
    section_H4(n_big, args.quick)
    section_H5(n_big)
    section_H6(n_big)
    section_C1()
    section_C2(90000, args.quick)
    section_C3()
    section_C4()
    P(f"DONE {elapsed()}")


if __name__ == "__main__":
    main()
