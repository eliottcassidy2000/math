"""LANE F1 -- THE ENTRY SCAN: feed-phase-only exact runner with per-row
ENTRY diagnostics and early stop at the Theorem S-cone-fc certificate.

Session: boxeph multifront, 2026-08-04.  All decisions exact int/Fraction;
floats only in display fields; no numpy; no sympy.

Semantics: transport/clamp/feed loop VERBATIM from the certified
amm12592_S_invariant_cone_certificates_boxeph.py fcscan (itself certified
bit-identical to the fast engine); initial_junk / floor_gamma_star imported
from the certified fast engine.  The ONLY functional differences:
  * feed coefficients of 2G_R are computed for the used index range
    [d_0 - 2, R-1] only (memory: the low coefficients are never read --
    feed indices d_i + i and d_i - 1 + i are >= d_0 - 1 and increase), and
    freed once consumed;
  * the run STOPS after the certificate row i_fc plus a persistence window
    (default 64 rows): Theorem S-cone-fc (PROVED) covers the remaining
    evolution, so the post-feed tail need not be simulated;
  * per-row ENTRY diagnostics (sign collapse, support, positive part,
    F2/F3 first rows, exact marginal-surface ratios).

Validation mode `validate` re-runs R=1024,D0=32 and R=4096,D0=128 and
compares (i_pf, d_fe, full feed-end junk dict) against the stored
amm12592_S_cone_feedend_* snapshots and (i_fc, m, a0, drain, Tmax,
capture_ub, budget_margin) against the stored amm12592_S_cone_fcscan_*
records -- all must be EXACTLY equal.

Output -> 05-knowledge/results/amm12592_entry_scan_R{R}_D0{D0}_boxeph.json
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb, isqrt
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

_spec = importlib.util.spec_from_file_location(
    "fastflow", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fastflow = importlib.util.module_from_spec(_spec)
with contextlib.redirect_stdout(io.StringIO()):
    _spec.loader.exec_module(fastflow)
initial_junk = fastflow.initial_junk
floor_gamma_star = fastflow.floor_gamma_star


def two_G_coeffs_top(R, j_start):
    """{j: [x^j] 2G_R} for j in [j_start, R-1] only (exact incremental)."""
    j_start = max(1, j_start)
    b = comb(R - 1, j_start)
    g = {}
    for j in range(j_start, R):
        g[j] = (-b if j & 1 else b) - 1
        b = b * (R - 1 - j) // (j + 1)
    return g


def caps_row(d, tmax):
    out = [2]
    P = 1
    for t in range(tmax):
        P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        out.append(2 * P)
    return out


def scan(R, D0, persist_rows=64, stride=256, postfeed_diag_rows=12,
         max_scan_past_pf=None, flush_every=1024):
    """Feed-phase ENTRY scan.  Returns result dict; writes ledger JSON."""
    t0 = time.time()
    d_prev = floor_gamma_star(R) + D0
    d0 = d_prev
    g = two_G_coeffs_top(R, d0 - 2)
    glow = max(1, d0 - 2)          # all indices < glow already freed
    j, junkL1, c0 = initial_junk(R, d_prev)
    minus2 = 1 if c0 == -2 else 0
    if max_scan_past_pf is None:
        max_scan_past_pf = R // 8
    i_pf = None
    d_fe = None
    fe_rec = None
    lastpos_row = 0 if any(v > 0 for v in j.values()) else None
    first_F2 = None
    first_F3 = None
    i_fc = None
    fc_rec = None
    persist_ok = True
    persist_seen = 0
    outcome = None
    traj = []
    pf_diags = []
    ledger_path = os.path.join(
        RESULTS, f"amm12592_entry_scan_R{R}_D0{D0}_boxeph.json")

    def sample(i, d):
        m = max(j) if j else -1
        npos = sum(1 for v in j.values() if v > 0)
        L1 = 0
        mb = 0
        for v in j.values():
            L1 += v if v > 0 else -v
            bl = (v if v > 0 else -v).bit_length()
            if bl > mb:
                mb = bl
        traj.append({"i": i, "d": d, "m": m, "ncells": len(j),
                     "npos": npos, "L1_bits": L1.bit_length(),
                     "max_bits": mb})

    def flush(final=False):
        json.dump({"R": R, "D0": D0, "i_pf": i_pf, "d_feedend": d_fe,
                   "feedend": fe_rec, "lastpos_row": lastpos_row,
                   "first_F2_row": first_F2, "first_F3_row": first_F3,
                   "i_fc": i_fc, "fc": fc_rec,
                   "F3_persist_rows_checked": persist_seen,
                   "F3_persist_ok": persist_ok,
                   "outcome": outcome, "postfeed_diags": pf_diags,
                   "trajectory": traj,
                   "elapsed_s": round(time.time() - t0, 1),
                   "final": final},
                  open(ledger_path, "w"))

    sample(0, d_prev)
    for i in range(1, R - 1):
        d = floor_gamma_star(R + i) + D0
        delta = d - d_prev
        assert delta in (0, 1), (i, d, d_prev)
        postfeed = (d + i > R)
        if postfeed and i_pf is None:
            i_pf = i
            d_fe = d_prev
            g = None                      # feed exhausted: free coefficients
            m = max(j) if j else -1
            a0 = -j.get(0, 0)
            npos = sum(1 for v in j.values() if v > 0)
            assert a0 == (R - 2) - 2 * minus2, ("T8 identity", i, a0, minus2)
            fe_rec = {"i_pf": i_pf, "d_fe": d_fe, "m": m, "a0": a0,
                      "a0_minus_dfe": a0 - d_fe, "ncells": len(j),
                      "npos": npos, "allneg": npos == 0,
                      "contiguous": sorted(j) == list(range(m + 1)),
                      "minus2": minus2,
                      "m_x1000_over_sqrtR": (m * 1000) // isqrt(R * 1000000)
                      if m >= 0 else None,
                      "junk": {str(t): str(v) for t, v in sorted(j.items())}
                      if m <= 400 else None}
        if postfeed and j:
            a = {t: -v for t, v in j.items()}
            m = max(a)
            a0 = a.get(0, 0)
            cap = caps_row(d, m + 3)
            F2 = a0 <= d - 1
            F3 = all(2 * a.get(t - 1, 0) + a.get(t - 2, 0) <= cap[t]
                     for t in range(2, m + 3))
            if F2 and first_F2 is None:
                first_F2 = i
            if F3 and first_F3 is None:
                first_F3 = i
            if len(pf_diags) < postfeed_diag_rows or (i_fc is not None and
                                                     i <= i_fc + 2):
                # exact marginal-surface profile: r_t = spill_t / cap_t
                worst_num, worst_den, worst_t = 0, 1, None
                over = []
                min_marg_num, min_marg_den, min_t = None, None, None
                for t in range(2, m + 3):
                    sp = 2 * a.get(t - 1, 0) + a.get(t - 2, 0)
                    if sp == 0 or cap[t] == 0:
                        continue
                    if sp * worst_den > worst_num * cap[t]:
                        worst_num, worst_den, worst_t = sp, cap[t], t
                    if sp > cap[t]:
                        over.append(t)
                    marg_n, marg_d = cap[t] - sp, cap[t]
                    if (min_marg_num is None or
                            marg_n * min_marg_den < min_marg_num * marg_d):
                        min_marg_num, min_marg_den, min_t = marg_n, marg_d, t
                pf_diags.append(
                    {"i": i, "d": d, "delta": delta, "m": m, "a0": a0,
                     "a0_vs_dminus1": (d - 1) - a0,
                     "npos": sum(1 for v in j.values() if v > 0),
                     "F2": F2, "F3": F3,
                     "max_r_float": (worst_num / worst_den
                                     if worst_den else None),
                     "max_r_t": worst_t, "over_cells": over,
                     "min_F3_margin_float": (min_marg_num / min_marg_den
                                             if min_marg_den else None),
                     "min_F3_margin_t": min_t})
            if i_fc is None and F2 and F3 and m + 2 < d:
                # exact clocks (verbatim fcscan semantics)
                drain = -((-a0) // 2)
                alive = [t for t in range(1, m + 1) if a.get(t, 0)]
                need = {t: a[t] for t in alive}
                cum = {t: 0 for t in alive}
                Tt = {}
                dd = d
                capk = cap[:]
                Kmax = R - 2 - i
                k = 0
                while need and k < Kmax:
                    k += 1
                    dn = floor_gamma_star(R + i + k) + D0
                    dl = dn - dd
                    if dl:
                        for t in range(0, m + 3):
                            capk[t] = capk[t] * dd // (dd - t)
                    a0k = max(0, a0 - 2 * (k - 1))
                    for t in list(need):
                        if t == 1:
                            dec = 2 * (dn - 1) - (1 + dl) * a0k
                            if dec > 0:
                                cum[1] += dec
                        else:
                            dec = capk[t] - cap[t]
                            if dec > 0:
                                cum[t] += dec
                        if cum[t] >= need[t]:
                            Tt[t] = k
                            del need[t]
                    dd = dn
                if not need:
                    Tmax = max(Tt.values()) if Tt else 0
                    ub = i + max(drain, Tmax)
                    if ub <= R - 2:
                        i_fc = i
                        fc_rec = {"i_fc": i, "i_fc_minus_i_pf": i - i_pf,
                                  "d": d, "m": m, "a0": a0, "drain": drain,
                                  "Tmax": Tmax,
                                  "T_worst_t": max(Tt, key=Tt.get) if Tt
                                  else None,
                                  "capture_ub": ub,
                                  "budget_margin": (R - 2) - ub}
                        flush()
            elif i_fc is not None:
                persist_seen += 1
                if not F3:
                    persist_ok = False       # would falsify the theorem
                if persist_seen >= persist_rows:
                    outcome = {"status": "CERT_STOP",
                               "note": "Theorem S-cone-fc covers the rest"}
                    flush(final=True)
                    break
        if i_pf is not None and i_fc is None and i - i_pf > max_scan_past_pf:
            outcome = {"status": "ENTRY_FAIL_WINDOW",
                       "window": max_scan_past_pf}
            flush(final=True)
            break
        # ---- transport (verbatim engine semantics)
        w = {}
        if delta == 0:
            for t, v in j.items():
                w[t] = w.get(t, 0) + v
                w[t + 1] = w.get(t + 1, 0) + v
        else:
            for t, v in j.items():
                w[t] = w.get(t, 0) + v
                w[t + 1] = w.get(t + 1, 0) + 2 * v
                w[t + 2] = w.get(t + 2, 0) + v
        fed = False
        if d + i <= R - 1:
            w[0] = w.get(0, 0) + g[d + i]; fed = True
        if delta == 1 and d - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[d - 1 + i]
            w[1] = w.get(1, 0) + g[d - 1 + i]; fed = True
        if g is not None:
            while glow <= d + i - 2:
                g.pop(glow, None)
                glow += 1
        # ---- clamp (verbatim engine semantics)
        jn = {}
        c0 = 0
        died = False
        if w:
            ts = sorted(w)
            ta, tb = ts[0], ts[-1]
            assert tb <= d, ("support beyond bottom cell", i, tb, d)
            P = comb(d - 1, ta)
            Pprev = comb(d - 1, ta - 1) if ta >= 1 else 0
            for t in range(ta, tb + 1):
                v = w.get(t, 0)
                if v:
                    lo, hi = -2 * P, 2 * Pprev
                    u = min(hi, max(lo, v))
                    assert (v - u) % 2 == 0 and v % 2 == 0, ("parity", i, t)
                    if u != v:
                        jn[t] = v - u
                    if t == 0:
                        c0 = u
                Pprev = P
                P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        if d in jn:
            outcome = {"status": "DIE", "row": i,
                       "const_bits": abs(jn[d]).bit_length()}
            died = True
        elif c0 == -2:
            minus2 += 1
        j = jn
        d_prev = d
        if any(v > 0 for v in j.values()):
            lastpos_row = i
        if i % stride == 0:
            sample(i, d)
        if died:
            flush(final=True)
            break
        if not j and not fed and d + i > R - 1:
            outcome = {"status": "CLOSED_EARLY", "capture_row": i}
            flush(final=True)
            break
        if i % flush_every == 0:
            flush()
    if outcome is None:
        outcome = {"status": "SCAN_END"}
        flush(final=True)
    print(f"[entry] R={R} D0={D0}: i_pf={i_pf} lastpos={lastpos_row} "
          f"firstF2={first_F2} firstF3={first_F3} i_fc={i_fc} fc={fc_rec} "
          f"persist_ok={persist_ok}({persist_seen}) outcome={outcome} "
          f"({round(time.time()-t0,1)}s)", flush=True)
    return {"R": R, "D0": D0, "i_pf": i_pf, "fe": fe_rec, "fc": fc_rec,
            "lastpos_row": lastpos_row, "outcome": outcome}


def validate():
    """Bit-exact agreement with the stored S-cone feedend + fcscan records."""
    ok_all = True
    for R, D0 in ((1024, 32), (4096, 128)):
        res = scan(R, D0, persist_rows=32)
        fe = json.load(open(os.path.join(
            RESULTS, f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json")))
        fc = json.load(open(os.path.join(
            RESULTS, f"amm12592_S_cone_fcscan_R{R}_D0{D0}_boxeph.json")))
        mine_fe = json.load(open(os.path.join(
            RESULTS, f"amm12592_entry_scan_R{R}_D0{D0}_boxeph.json")))
        checks = {
            "i_pf": res["i_pf"] == fe["i_postfeed"],
            "d_fe": mine_fe["d_feedend"] == fe["d_feedend"],
            "junk_dict": mine_fe["feedend"]["junk"] ==
                         {t: str(v) for t, v in fe["junk"].items()},
            "i_fc": res["fc"]["i_fc"] == fc["i_fc"],
            "m": res["fc"]["m"] == fc["fc"]["m"],
            "a0": res["fc"]["a0"] == fc["fc"]["a0"],
            "drain": res["fc"]["drain"] == fc["fc"]["drain"],
            "Tmax": res["fc"]["Tmax"] == fc["fc"]["Tmax"],
            "capture_ub": res["fc"]["capture_ub"] == fc["fc"]["capture_ub"],
            "budget_margin": res["fc"]["budget_margin"] ==
                             fc["fc"]["budget_margin"],
        }
        ok = all(checks.values())
        ok_all = ok_all and ok
        print(f"[validate] R={R} D0={D0}: "
              f"{'ALL-PASS' if ok else ('FAIL ' + str(checks))}", flush=True)
    print(f"[validate] overall: {'ALL-PASS' if ok_all else 'FAIL'}",
          flush=True)
    return ok_all


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        validate()
    else:
        R, D0 = int(sys.argv[1]), int(sys.argv[2])
        scan(R, D0)
