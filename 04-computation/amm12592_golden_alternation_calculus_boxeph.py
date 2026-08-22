"""LANE D3 -- THE GOLDEN TRANSIENT BOUND, part 1: the alternation-transport
calculus (boxeph, 2026-08-03).  Exact int arithmetic throughout; no floats in
any decision; no numpy.

Certifies the algebra behind Lemmas G1/G1'/G2 of
05-knowledge/results/amm12592-golden-transient-bound-boxeph.md:

  G1  (transport identity)  For j_t = (-1)^(d-t) a_t and the T6 Pascal kernel
      K_delta ((1,1) for delta=0, (1,2,1) for delta=1, acting by
      (K*j)_c = sum_s K_s j_{c-s}):
          (K*j)_c = (-1)^(d-c) (D^{1+delta} a)_c ,
      D = backward difference (Da)_c = a_c - a_{c-1}.  In the NEW phase
      (-1)^(d'-c), d' = d+delta:  image = (-1)^delta * (-1)^(d'-c) (D^{1+delta}a)_c
      -- exactly alternating again iff D^{1+delta} a is single-signed.
  G1z (zero-sum + defect-mass law)  sum_c (D^{1+delta}a)_c = 0 for finitely
      supported a; hence  sum_c |(K*j)_c| = 2 sum_c ((D^{1+delta}a)_c)^+
      = 2 sum_c ((D^{1+delta}a)_c)^- :  the transported L1 EQUALS twice the
      one-signed (defect) mass of the difference profile.
  G1r (ratio law)  geometric magnitudes a_{t-1}/a_t = rho (inward growth):
      interior image magnitude = (rho-1)^{1+delta} * a_c exactly; rho = 2 is
      magnitude-preserving for the (1,2,1) kernel ((rho-1)^2 = 1).
  G1x (hazard record)  the naive claim "|j_t| <= M_t alternating implies
      |(K*j)| <= |D^2 M|" is FALSE: explicit counterexample recorded.
  G2  (initial data are exactly alternating)  the row-0 junk of the T4 closed
      form at d_0 = floor(g R) + D0 obeys sign(j_t) = (-1)^(d-t) on its full
      support [0, R-2-d_0] (T4b corollary, re-verified independently), and its
      magnitude profile has the exact ratio-2 crossover near t2 = 2d-R+3:
      cells t < t2 have a_{t-1} < 2 a_t ("sub-2" = kernel-contracting band),
      cells t > t2 have a_{t-1} > 2 a_t ("super-2" = kernel-expanding band).
      Additive convexity signature of a reported (D^2 a sign profile).
  G3pre (mechanism traces)  plain-rule replays at adjacent die/close pairs:
      first interior alternation defect (adjacent same-sign pair) row+cell,
      first unabsorbed edge transport (front advance), FDA margin profile.

Outputs: 05-knowledge/results/amm12592_golden_alternation_calculus_boxeph.out
         05-knowledge/results/amm12592_golden_alternation_calculus_boxeph.json
"""
import sys, os, json, time, io, contextlib, importlib.util, random
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
OUT = os.path.join(RESULTS, "amm12592_golden_alternation_calculus_boxeph.out")
JS = os.path.join(RESULTS, "amm12592_golden_alternation_calculus_boxeph.json")

_argv = sys.argv[:]
sys.argv = [_argv[0], "noop"]
spec = importlib.util.spec_from_file_location(
    "fastflow", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fastflow = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fastflow)
sys.argv = _argv
floorgs = fastflow.floor_gamma_star
two_G_coeffs = fastflow.two_G_coeffs
initial_junk = fastflow.initial_junk

RES = {}
if os.path.exists(JS):
    try:
        RES = json.load(open(JS))
    except Exception:
        RES = {}


def log(line):
    print(line, flush=True)
    with open(OUT, "a") as f:
        f.write(line + "\n")


def save():
    json.dump(RES, open(JS, "w"), indent=1)


def kernel_image(j, delta):
    """(K_delta * j)_c = sum_s K_s j_{c-s}, K = (1,1) or (1,2,1)."""
    K = (1, 1) if delta == 0 else (1, 2, 1)
    w = {}
    for t, v in j.items():
        for s, ks in enumerate(K):
            w[t + s] = w.get(t + s, 0) + ks * v
    return {c: v for c, v in w.items() if v}


def diff_profile(a, order, lo, hi):
    """(D^order a)_c = sum_s (-1)^s C(order,s) a_{c-s}, dict a, cells lo..hi."""
    out = {}
    for c in range(lo, hi + 1):
        v = sum((-1) ** s * comb(order, s) * a.get(c - s, 0)
                for s in range(order + 1))
        if v:
            out[c] = v
    return out


# ---------------------------------------------------------------- G1 identity
def cert_G1():
    rng = random.Random(12592)
    n_id = n_zero = n_l1 = 0
    for trial in range(400):
        d = rng.randrange(6, 60)
        delta = trial & 1
        N = rng.randrange(2, min(d - 2, 40))
        a = {t: rng.randrange(-10**9, 10**9) for t in range(N + 1)}
        j = {t: (-1) ** ((d - t) & 1) * a[t] for t in a if a[t]}
        w = kernel_image(j, delta)
        D = diff_profile(a, 1 + delta, 0, N + 2)
        # identity in the OLD phase
        ok_id = all(w.get(c, 0) == (-1) ** ((d - c) & 1) * D.get(c, 0)
                    for c in range(0, N + 3))
        # zero-sum + defect-mass law
        tot = sum(D.values())
        l1 = sum(abs(v) for v in w.values())
        pos = sum(v for v in D.values() if v > 0)
        neg = -sum(v for v in D.values() if v < 0)
        ok_zero = (tot == 0)
        ok_l1 = (l1 == 2 * pos == 2 * neg) if l1 else (pos == neg == 0)
        n_id += ok_id
        n_zero += ok_zero
        n_l1 += ok_l1
        assert ok_id and ok_zero and ok_l1, ("G1", trial, d, delta)
    # ratio law on exact integer geometric profiles
    ratio_checks = []
    for rho in (2, 3, 4, 5):
        for delta in (0, 1):
            N = 30
            d = 40
            a = {t: rho ** (N - t) for t in range(N + 1)}   # a_{t-1}/a_t = rho
            j = {t: (-1) ** ((d - t) & 1) * a[t] for t in a}
            w = kernel_image(j, delta)
            # interior cells 2+delta .. N: |image_c| == (rho-1)^{1+delta} a_c
            ok = all(abs(w.get(c, 0)) == (rho - 1) ** (1 + delta) * a[c]
                     for c in range(1 + delta, N + 1))
            ratio_checks.append({"rho": rho, "delta": delta,
                                 "interior_multiplier_exact": ok})
            assert ok, ("G1r", rho, delta)
    # rho = 2, delta = 1: magnitude preserving exactly ((rho-1)^2 = 1)
    # hazard record: naive envelope claim is FALSE
    d, N = 20, 8
    M = {t: 100 for t in range(N + 1)}
    a = {t: (100 if t % 2 == 0 else 0) for t in range(N + 1)}
    j = {t: (-1) ** ((d - t) & 1) * a[t] for t in a if a[t]}
    w = kernel_image(j, 1)
    D2M = diff_profile(M, 2, 0, N + 2)
    viol = [(c, abs(w.get(c, 0)), abs(D2M.get(c, 0)))
            for c in range(0, N + 3) if abs(w.get(c, 0)) > abs(D2M.get(c, 0))]
    assert viol, "expected counterexample to the naive envelope claim"
    RES["G1"] = {"identity_trials": n_id, "zero_sum_trials": n_zero,
                 "defect_mass_L1_trials": n_l1, "ratio_law": ratio_checks,
                 "naive_envelope_counterexample":
                     {"a": "M on even cells, 0 on odd; M=100, N=8",
                      "first_violations": viol[:4],
                      "verdict": "naive claim FALSE (hazard record)"}}
    log(f"G1: identity {n_id}/400, zero-sum {n_zero}/400, "
        f"defect-mass-L1 {n_l1}/400, ratio law "
        f"{sum(r['interior_multiplier_exact'] for r in ratio_checks)}/8, "
        f"naive-envelope counterexample: {len(viol)} violating cells "
        f"(e.g. {viol[0]})")
    save()


# ----------------------------------------------------- G2 initial-data structure
def cert_G2():
    grid = [(256, 0), (256, 1), (512, 0), (512, 5), (1024, 0), (1024, 15),
            (2048, 38), (4096, 89), (8192, 192)]
    recs = []
    for R, D0 in grid:
        d = floorgs(R) + D0
        j, junkL1, c0 = initial_junk(R, d)
        supp = sorted(j)
        F = supp[-1]
        # sign law: j_t = (-1)^(d-t) |j_t| on the FULL support
        sign_ok = all((j[t] > 0) == (((d - t) & 1) == 0) for t in supp)
        # support exactly [0, R-2-d] (T4b)
        supp_ok = (supp == list(range(0, R - 2 - d + 1)))
        a = {t: abs(j[t]) for t in supp}
        # ratio-2 crossover: exact comparisons a_{t-1} ? 2 a_t
        t2_pred = 2 * d - R + 3
        below = [t for t in range(1, F + 1) if a[t - 1] < 2 * a[t]]
        above = [t for t in range(1, F + 1) if a[t - 1] > 2 * a[t]]
        # transition window: last sub-2 cell and first super-2 cell
        last_sub2 = max(below) if below else None
        first_sup2 = min(above) if above else None
        # is the ratio profile monotone in the comparison sense: all cells
        # below some t* sub-2, all above super-2 (up to a small mixed window)?
        mixed = [t for t in range(1, F + 1)
                 if (t <= (t2_pred - 8) and a[t - 1] > 2 * a[t]) or
                    (t >= (t2_pred + 8) and a[t - 1] < 2 * a[t])]
        # additive convexity: sign of (D^2 a)_t on 2..F
        conc = [t for t in range(2, F + 1)
                if a[t] - 2 * a[t - 1] + a[t - 2] < 0]
        conv = [t for t in range(2, F + 1)
                if a[t] - 2 * a[t - 1] + a[t - 2] > 0]
        rec = {"R": R, "D0": D0, "d": d, "front": F,
               "sign_law_ok": sign_ok, "support_T4b_ok": supp_ok,
               "t2_pred_2d-R+3": t2_pred,
               "last_sub2_cell": last_sub2, "first_super2_cell": first_sup2,
               "cells_outside_pred_window8": len(mixed),
               "n_concave_cells": len(conc), "n_convex_cells": len(conv),
               "concave_range": ([min(conc), max(conc)] if conc else None),
               "convex_range": ([min(conv), max(conv)] if conv else None)}
        recs.append(rec)
        log(f"G2 R={R} D0={D0} d={d}: sign_law={sign_ok} supp=[0,{F}] "
            f"(T4b {supp_ok}); ratio-2 crossover pred {t2_pred}, "
            f"last_sub2={last_sub2} first_super2={first_sup2} "
            f"mixed_outside_w8={len(mixed)}; "
            f"D2a<0 on {len(conc)} cells {rec['concave_range']}, "
            f"D2a>0 on {len(conv)} cells {rec['convex_range']}")
        save()
    RES["G2"] = recs
    save()


# ------------------------------------------- G3pre: plain-run mechanism traces
def plain_replay(R, D0, maxrows=None, keep_stats_rows=None):
    """Plain rule replay with full junk; per-row alternation/defect stats."""
    g = two_G_coeffs(R)
    d = floorgs(R) + D0
    j, _, _ = initial_junk(R, d)
    stats = []
    first_defect = None
    first_front_adv = None
    outcome = None
    for i in range(1, R - 1):
        dn = floorgs(R + i) + D0
        delta = dn - d
        w = kernel_image(j, delta)
        if dn + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn + i]
        if delta == 1 and dn - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn - 1 + i]
            w[1] = w.get(1, 0) + g[dn - 1 + i]
        w = {t: v for t, v in w.items() if v}
        prev_front = max(j) if j else None
        jn = {}
        for t in sorted(w):
            v = w[t]
            lo = -2 * comb(dn - 1, t)
            hi = 2 * comb(dn - 1, t - 1) if t >= 1 else 0
            u = min(hi, max(lo, v))
            if u != v:
                jn[t] = v - u
        if dn in jn:
            outcome = {"DIE": i}
        # stats
        supp = sorted(jn)
        defects = [t for k, t in enumerate(supp[:-1])
                   if supp[k + 1] == t + 1 and
                   (jn[t] > 0) == (jn[t + 1] > 0)]
        if defects and first_defect is None:
            first_defect = {"row": i, "cells": defects[:6],
                            "front": supp[-1] if supp else None}
        if prev_front is not None and supp and supp[-1] == prev_front + 1 + delta \
                and first_front_adv is None:
            first_front_adv = {"row": i, "new_front": supp[-1]}
        if keep_stats_rows is None or i in keep_stats_rows:
            # FDA margin at the front: landing-cell end minus |edge transport|
            fda = None
            if prev_front is not None:
                land = prev_front + 1 + delta
                v_edge = j.get(prev_front, 0)
                end = 2 * comb(dn - 1, land) if v_edge < 0 else \
                    2 * comb(dn - 1, land - 1)
                fda = (end - abs(v_edge)).bit_length() if end >= abs(v_edge) \
                    else -((abs(v_edge) - end).bit_length())
            stats.append({"i": i, "d": dn, "front": (supp[-1] if supp else None),
                          "ndefects": len(defects),
                          "fda_sgnbits": fda})
        j = jn
        d = dn
        if outcome:
            break
        if not j and dn + i > R - 1:
            outcome = {"CLOSED_capture": i}
            break
        if maxrows and i >= maxrows:
            outcome = {"STOPPED": i}
            break
    return {"R": R, "D0": D0, "outcome": outcome,
            "first_interior_defect": first_defect,
            "first_front_advance": first_front_adv, "rows": stats}


def cert_G3pre():
    runs = [(512, 4), (512, 5), (1024, 14), (1024, 15)]
    recs = []
    for R, D0 in runs:
        rows = set(range(0, 61)) | set(range(60, 700, 10))
        r = plain_replay(R, D0, keep_stats_rows=rows)
        recs.append(r)
        fd = r["first_interior_defect"]
        fa = r["first_front_advance"]
        d0 = floorgs(R) + D0
        t2 = 2 * d0 - R + 3
        log(f"G3pre plain R={R} D0={D0}: outcome={r['outcome']} "
            f"first_defect={fd} first_front_advance={fa} "
            f"(row-0 crossover t2={t2}, row-0 front={R-2-d0})")
        save()
    RES["G3pre"] = recs
    save()


if __name__ == "__main__":
    t0 = time.time()
    log("=" * 72)
    log("LANE D3 script A: alternation-transport calculus certificates "
        "(exact ints)")
    stages = sys.argv[1:] or ["G1", "G2", "G3pre"]
    if "G1" in stages:
        cert_G1()
    if "G2" in stages:
        cert_G2()
    if "G3pre" in stages:
        cert_G3pre()
    log(f"done in {round(time.time()-t0,1)}s")
    save()
