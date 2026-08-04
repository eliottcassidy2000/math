"""LANE E1 -- Hypothesis S via the INVARIANT CONE: instrumented exact runner,
entry certificates, preservation margins, corner/robustness experiments.

Session: boxeph multifront, 2026-08-03.  All decisions exact int; no floats
in any decision (floats appear only in *display* fields of ledgers, never in
comparisons); no numpy; no sympy.

Semantics: identical to the certified fast engine
amm12592_transient_fast_junkflow_boxeph.py (T2+T6 conjugacy chain); the
initializer (initial_junk), 2G_R coefficients and the Lucas/Fibonacci floor
engine are IMPORTED from it, and the transport/clamp row loop is copied
verbatim (certified by mode `certify`: per-row trace equality + outcome
equality against run_fast).

New content (Lane E1):
  * post-feed detection: first row i with d_i + i > R (both feed branches
    dead forever after, since d_i + i is strictly increasing);
  * feed-end state snapshot (sparse junk, exact ints) + ENTRY certificate
    against the cone constants (LAM, MU, C0);
  * per-row POST-FEED cone diagnostics (exact int comparisons; margins as
    bit-length differences):
      - all-negative (assert), support max m, contiguity,
      - a_0 clock: a_0 vs d-1 (cell-1 spill criterion), a_0 vs 2(R-2-i),
      - Lambda-bound: a_t <= (LAM-1)*2C(d-1,t) for all t>=1,
      - coupling margins: for t in [2..m]: 2*spill_t <= capnext_t where
        spill_t = 2a_{t-1}+a_{t-2} (delta=1) or a_{t-1} (delta=0),
        capnext_t = 2C(d_i - 1, t) (the cap of the row being computed)
        -- the PROVED preservation lemma's hypothesis (coupling <= 1/2),
      - cell-1 coupling exact: (1+delta)*a_0 vs capnext_1,
      - freeze margins at the front (C-F),
      - extinction events (cell death rows, re-ignition detector);
  * `corner` mode: run the post-feed flow from the CONE CORNER state
    (a_t = (LAM-1)*2C(d-1,t) on t in [1, m0], a_0 = d + C0) -- the
    comparison-lemma certificate C(R) experiment;
  * `fromsnap` mode: run from a saved feed-end snapshot scaled up by 2^k
    (robustness radius experiment, comparison lemma).

Outputs -> 05-knowledge/results/amm12592_S_cone_{run,entry,corner}_*.json
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

_spec = importlib.util.spec_from_file_location(
    "fastflow", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fastflow = importlib.util.module_from_spec(_spec)
with contextlib.redirect_stdout(io.StringIO()):
    _spec.loader.exec_module(fastflow)
two_G_coeffs = fastflow.two_G_coeffs
initial_junk = fastflow.initial_junk
floor_gamma_star = fastflow.floor_gamma_star

# Cone constants (the S-entry hypothesis instantiation)
LAM = 2 ** 11          # uniform cap-ratio bound Lambda
MU = 2                 # front bound m <= MU * ceil(sqrt(R))
C0 = 64                # a_0 <= d + C0 at feed-end


def caps_row(d, tmax):
    """[2*C(d-1,t) for t in 0..tmax], exact incremental."""
    out = [2]
    P = 1
    for t in range(tmax):
        P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        out.append(2 * P)
    return out


def bits(n):
    return n.bit_length()


def entry_certificate(R, D0, i_pf, d_pf, j):
    """ENTRY check at the feed-end state (junk entering row i_pf, degree
    d_pf = d_{i_pf - 1}).  All comparisons exact ints."""
    m = max(j) if j else -1
    a0 = -j.get(0, 0)
    cap = caps_row(d_pf, m + 1)
    allneg = all(v < 0 for v in j.values())
    contig = (sorted(j) == list(range(0, m + 1))) if j else True
    sq = isqrt(R)
    m0 = MU * (sq if sq * sq == R else sq + 1)
    # E-A: a_t <= (LAM-1) * 2C(d-1,t) for t >= 1; record worst ratio bits
    ea_ok = True
    worst = None            # (t, bits(a_t) - bits(cap_t)) max
    maxA_bits = None
    for t, v in j.items():
        if t == 0:
            continue
        a = -v
        if a > (LAM - 1) * cap[t]:
            ea_ok = False
        rb = bits(a) - bits(cap[t])
        if maxA_bits is None or rb > maxA_bits:
            maxA_bits = rb
            worst = t
    rec = {"R": R, "D0": D0, "i_postfeed": i_pf, "d_feedend": d_pf,
           "m": m, "sqrtR_x1000": isqrt(R * 1000000),
           "a0": a0, "a0_minus_d": a0 - d_pf,
           "E_N_allneg": allneg, "E_C_contiguous": contig,
           "E_F_front_ok": m <= m0, "m0": m0,
           "E_A_lambda_ok": ea_ok, "maxA_bits_excess_over_cap": maxA_bits,
           "maxA_argmax_t": worst,
           "E_0_debt_ok": a0 <= d_pf + C0 and a0 <= R - 2,
           "deadline_slack": 2 * (R - 2 - i_pf) - a0,
           "prod_6Lam_m2_le_d": 6 * LAM * (m + 2) <= d_pf,
           "prod_6Lam_m2_bits": bits(6 * LAM * (m + 2)) - bits(d_pf)}
    rec["ENTRY_PASS"] = (allneg and m <= m0 and ea_ok and
                         a0 <= d_pf + C0 and a0 <= R - 2)
    return rec


def postfeed_diag(i, d, delta, j, R):
    """Cone diagnostics for the state j ENTERING row i (degree of the row
    being computed is d = d_i, so the caps below are 2C(d-1, t))."""
    if not j:
        return {"i": i, "empty": True}
    m = max(j)
    a0 = -j.get(0, 0)
    cap = caps_row(d, m + 3)
    a = {t: -v for t, v in j.items()}
    # coupling for cells t in [2..m]: 2*spill <= cap_t ?  (PROVED lemma hyp)
    worst_t, worst_num, worst_den, coup_ok = None, None, None, True
    for t in range(2, m + 1):
        if delta == 1:
            spill = 2 * a.get(t - 1, 0) + a.get(t - 2, 0)
        else:
            spill = a.get(t - 1, 0)
        if spill == 0:
            continue
        if 2 * spill > cap[t]:
            coup_ok = False
        rb = bits(2 * spill) - bits(cap[t])
        if worst_num is None or rb > worst_num:
            worst_num, worst_t = rb, t
    # cell-1 coupling: (1+delta)*a0 <= cap_1 = 2(d-1) ?
    c1_ok = (1 + delta) * a0 <= cap[1]
    # Lambda bound on live cells
    lam_ok = all(a[t] <= (LAM - 1) * cap[t] for t in a if t >= 1)
    # freeze (C-F) at the front
    if delta == 1:
        fr_ok = (a.get(m, 0) <= comb(d, m + 2) * 2 and
                 2 * a.get(m, 0) + a.get(m - 1, 0) <= 2 * comb(d, m + 1))
    else:
        fr_ok = a.get(m, 0) <= cap[m + 1]
    return {"i": i, "d": d, "delta": delta, "m": m, "nlive": len(j),
            "a0": a0, "a0_vs_dminus1": (d - 1) - a0,
            "deadline_slack": 2 * (R - 2 - i) - a0,
            "coupling_ok": coup_ok, "coupling_worst_bits": worst_num,
            "coupling_worst_t": worst_t, "cell1_spill_ok": c1_ok,
            "lambda_ok": lam_ok, "freeze_ok": fr_ok}


def flow(R, D0, j=None, i_start=None, d_start=None, tag="run",
         ledger_stride=8, flush_every=512, diag_from_postfeed=True,
         trace_keys=False):
    """Exact rule-A junk flow.  If j/i_start/d_start given, start there
    (post-feed restart for corner/robustness modes); else full run from
    row 0.  Returns result dict; writes partial ledgers.
    Transport/clamp semantics IDENTICAL to fastflow.run_fast."""
    t0 = time.time()
    g = two_G_coeffs(R)
    if j is None:
        d_prev = floor_gamma_star(R) + D0
        j, junkL1, c0 = initial_junk(R, d_prev)
        i0 = 1
    else:
        d_prev = d_start
        i0 = i_start
        c0 = 0
    front0 = max(j) if j else -1
    minus2 = 0
    if i0 == 1 and c0 == -2:
        minus2 = 1
    i_pf = None
    snap = None
    entry = None
    diags = []
    events = []            # (row, t, "die"/"born")
    alive_prev = set(t for t in j if t >= 1)
    trace = []
    outcome = None
    capture = None
    ledger_path = os.path.join(
        RESULTS, f"amm12592_S_cone_{tag}_R{R}_D0{D0}_boxeph.json")

    def flush(final=False):
        rec = {"R": R, "D0": D0, "tag": tag, "i_postfeed": i_pf,
               "outcome": outcome, "capture_row": capture,
               "minus2_rows": minus2, "front0": front0,
               "entry_certificate": entry, "events": events[:4000],
               "n_events": len(events), "diags": diags,
               "elapsed_s": round(time.time() - t0, 1),
               "final": final}
        if trace_keys:
            rec["trace"] = trace
        json.dump(rec, open(ledger_path, "w"))

    for i in range(i0, R - 1):
        d = floor_gamma_star(R + i) + D0
        delta = d - d_prev
        assert delta in (0, 1), (i, d, d_prev)
        postfeed = (d + i > R)
        if postfeed and i_pf is None:
            i_pf = i
            snap = dict(j)
            entry = entry_certificate(R, D0, i_pf, d_prev, j)
            if tag == "run":     # only genuine runs may write the snapshot
                json.dump({"R": R, "D0": D0, "i_postfeed": i_pf,
                           "d_feedend": d_prev, "entry": entry,
                           "junk": {str(t): str(v) for t, v in snap.items()}},
                          open(os.path.join(RESULTS,
                            f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json"),
                            "w"))
        if postfeed and diag_from_postfeed and (
                (i - i_pf) % ledger_stride == 0 or (i - i_pf) < 64):
            diags.append(postfeed_diag(i, d, delta, j, R))
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
        # ---- clamp (verbatim engine semantics)
        jn = {}
        junkL1 = 0
        c0 = 0
        w_bottom = w.get(d, 0)
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
                        junkL1 += abs(v - u)
                    if t == 0:
                        c0 = u
                Pprev = P
                P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        if d in jn:
            outcome = {"status": "DIE", "row": i,
                       "const_bits": abs(jn[d]).bit_length()}
            died = True
        else:
            if c0 == -2:
                minus2 += 1
        # ---- post-feed invariance asserts + extinction bookkeeping
        if postfeed:
            assert all(v < 0 for v in jn.values()), ("positivity broken", i)
            alive = set(t for t in jn if t >= 1)
            for t in sorted(alive_prev - alive):
                events.append((i, t, "die"))
            for t in sorted(alive - alive_prev):
                events.append((i, t, "born"))
            alive_prev = alive
        else:
            alive_prev = set(t for t in jn if t >= 1)
        if trace_keys and (i % 64 == 0 or died or not jn):
            ts2 = sorted(jn)
            trace.append({"i": i, "d": d, "nclamped": len(jn),
                          "tmin": ts2[0] if ts2 else None,
                          "tmax": ts2[-1] if ts2 else None,
                          "junkL1_bits": junkL1.bit_length(), "c0": c0,
                          "e_in0": w_bottom})
        j = jn
        d_prev = d
        if died:
            break
        if not j and not fed and d + i > R - 1:
            capture = i
            outcome = {"status": "CLOSED", "capture_row": i}
            break
        if i % flush_every == 0:
            flush()
    if outcome is None:
        outcome = {"status": "CLOSED", "capture_row": R - 2} if not j \
            else {"status": "OPEN_RESIDUAL"}
        if not j:
            capture = R - 2
    # final summary diagnostics
    flush(final=True)
    return {"R": R, "D0": D0, "tag": tag, "outcome": outcome,
            "capture_row": capture, "minus2_rows": minus2, "i_pf": i_pf,
            "entry": entry, "n_events": len(events),
            "elapsed_s": round(time.time() - t0, 1)}


# ------------------------------------------------------------- certification
def certify():
    """Per-row equality vs the certified engine at 4 configs incl. a death."""
    out = {}
    for R, D0 in ((128, 8), (256, 16), (512, 32), (512, 4)):
        fast = fastflow.run_fast(R, D0)
        mine = flow(R, D0, tag="cert", diag_from_postfeed=True,
                    trace_keys=True)
        ok = (fast["outcome"] == mine["outcome"] and
              fast["capture_row"] == mine["capture_row"] and
              fast["minus2_rows"] == mine["minus2_rows"])
        # spot per-row: compare my sparse trace rows against engine trace
        ftr = {r["i"]: r for r in fast["trace"]}
        led = json.load(open(os.path.join(
            RESULTS, f"amm12592_S_cone_cert_R{R}_D0{D0}_boxeph.json")))
        bad = 0
        for r in led.get("trace", []):
            fr = ftr.get(r["i"])
            if fr is None:
                continue
            for k in ("d", "nclamped", "tmin", "tmax", "junkL1_bits", "c0",
                      "e_in0"):
                if fr.get(k) != r.get(k):
                    bad += 1
        out[f"R{R}_D0{D0}"] = {"outcome_equal": ok, "row_mismatches": bad,
                               "fast": fast["outcome"],
                               "mine": mine["outcome"]}
        print(f"[cert] R={R} D0={D0}: outcome_equal={ok} "
              f"row_mismatches={bad} ({mine['outcome']})", flush=True)
    json.dump(out, open(os.path.join(
        RESULTS, "amm12592_S_cone_certify_vs_engine_boxeph.json"), "w"),
        indent=1)
    return out


def corner(R, D0, lam=None, mu=None, c0v=None, scale_bits=0):
    """Comparison-lemma corner certificate C(R): start the post-feed flow
    from the cone corner state at the feed-end row (no feed run needed:
    i_pf and d are computed from the floor profile directly)."""
    lam = lam if lam is not None else LAM
    mu = mu if mu is not None else MU
    c0v = c0v if c0v is not None else C0
    # find i_pf = first i with d_i + i > R
    lo, hi = 0, R - 2
    while lo < hi:
        mid = (lo + hi) // 2
        if floor_gamma_star(R + mid) + D0 + mid > R:
            hi = mid
        else:
            lo = mid + 1
    i_pf = lo
    d_pf = floor_gamma_star(R + i_pf - 1) + D0     # degree of prior row
    sq = isqrt(R)
    m0 = mu * (sq if sq * sq == R else sq + 1)
    cap = caps_row(d_pf, m0)
    j = {0: -2 * ((d_pf + c0v) // 2)}
    for t in range(1, m0 + 1):
        j[t] = -(lam - 1) * cap[t] * (1 << scale_bits)
    res = flow(R, D0, j=j, i_start=i_pf, d_start=d_pf,
               tag=f"corner_l{lam}_m{mu}_s{scale_bits}")
    print(f"[corner] R={R} D0={D0} lam={lam} mu={mu} scale=2^{scale_bits}: "
          f"{res['outcome']} capture={res['capture_row']} "
          f"({res['elapsed_s']}s)", flush=True)
    return res


def fromsnap(R, D0, scale_bits):
    """Robustness: rerun post-feed from the saved feed-end snapshot with all
    cells t>=1 multiplied by 2^scale_bits (a_0 kept: its clock is exact)."""
    p = os.path.join(RESULTS,
                     f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json")
    S = json.load(open(p))
    j = {int(t): int(v) * ((1 << scale_bits) if int(t) >= 1 else 1)
         for t, v in S["junk"].items()}
    res = flow(R, D0, j=j, i_start=S["i_postfeed"], d_start=S["d_feedend"],
               tag=f"snapx{scale_bits}")
    print(f"[fromsnap] R={R} D0={D0} x2^{scale_bits}: {res['outcome']} "
          f"capture={res['capture_row']} ({res['elapsed_s']}s)", flush=True)
    return res


def iconescan(R, D0, from_snap=False):
    """Mid-flight cone-entry scan: find the first post-feed row i_cone at
    which the CURRENT state satisfies the Theorem S-cone hypotheses
    (post-layer + half-cap H4* + BUD with current-state clocks), then
    verify persistence of the half-cap condition and capture.
    All decisions exact."""
    from fractions import Fraction
    g_hi = Fraction(156759, 262144)
    t0 = time.time()
    g = two_G_coeffs(R)
    if from_snap:
        S = json.load(open(os.path.join(
            RESULTS, f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json")))
        j = {int(t): int(v) for t, v in S["junk"].items()}
        i0, d_prev = S["i_postfeed"], S["d_feedend"]
    else:
        d_prev = floor_gamma_star(R) + D0
        j, _, _ = initial_junk(R, d_prev)
        i0 = 1
    i_cone = None
    cone_rec = None
    persist_ok = True
    outcome = None
    capture = None
    for i in range(i0, R - 1):
        d = floor_gamma_star(R + i) + D0
        delta = d - d_prev
        postfeed = (d + i > R)
        if postfeed and j:
            a = {t: -v for t, v in j.items()}
            m = max(a)
            a0 = a.get(0, 0)
            cap = caps_row(d, m + 3)
            postlayer = a0 <= d - 1
            h4 = all(2 * (2 * a.get(t - 1, 0) + a.get(t - 2, 0)) <= cap[t]
                     for t in range(2, m + 3))
            if i_cone is None:
                if postlayer and h4:
                    drain = -((-a0) // 2)
                    K2 = 0
                    for t in range(2, m + 1):
                        v = a.get(t, 0)
                        if v:
                            K2 = max(K2, -((-v) // (cap[t] // 2)))
                    a1 = a.get(1, 0)
                    if a1:
                        D0c = 2 * (d - 1) - a0
                        # n*D0c + n(n-1) >= a1  (quadratic, exact ceil)
                        n = 0
                        lo, hi = 0, a1 // max(1, D0c) + 2
                        while lo < hi:
                            mid = (lo + hi) // 2
                            if mid * D0c + mid * (mid - 1) >= a1:
                                hi = mid
                            else:
                                lo = mid + 1
                        n = lo
                        K1 = int(-((-(n + 1)) // (1 - g_hi)))
                    else:
                        K1 = 0
                    ub = i + max(drain, max(K1, K2))
                    if ub <= R - 2:
                        i_cone = i
                        cone_rec = {"i_cone": i, "d": d, "m": m, "a0": a0,
                                    "drain": drain, "K1": K1, "K2": K2,
                                    "capture_ub": ub,
                                    "budget_margin": (R - 2) - ub}
            else:
                if not h4:
                    persist_ok = False   # would falsify the theorem
        # transport + clamp (engine semantics)
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
        jn = {}
        if w:
            ts = sorted(w)
            ta, tb = ts[0], ts[-1]
            P = comb(d - 1, ta)
            Pprev = comb(d - 1, ta - 1) if ta >= 1 else 0
            for t in range(ta, tb + 1):
                v = w.get(t, 0)
                if v:
                    lo, hi = -2 * P, 2 * Pprev
                    u = min(hi, max(lo, v))
                    if u != v:
                        jn[t] = v - u
                Pprev = P
                P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        if d in jn:
            outcome = {"status": "DIE", "row": i}
            break
        j = jn
        d_prev = d
        if not j and not fed and d + i > R - 1:
            capture = i
            outcome = {"status": "CLOSED", "capture_row": i}
            break
    res = {"R": R, "D0": D0, "outcome": outcome, "capture_row": capture,
           "i_cone": i_cone, "cone": cone_rec, "h4_persist_ok": persist_ok,
           "elapsed_s": round(time.time() - t0, 1)}
    json.dump(res, open(os.path.join(
        RESULTS, f"amm12592_S_cone_iconescan_R{R}_D0{D0}_boxeph.json"), "w"),
        indent=1)
    print(f"[icone] R={R} D0={D0}: i_cone={i_cone} rec={cone_rec} "
          f"persist={persist_ok} outcome={outcome} ({res['elapsed_s']}s)",
          flush=True)
    return res


def fcscan(R, D0, from_snap=False, i_stop_scan=None):
    """FULL-CAP cone scan (Theorem S-cone-fc): find the first post-feed row
    i_fc whose state satisfies
      F1 j <= 0, supp in [0,m], m+2 < d;
      F2 a_0 <= d - 1;
      F3 2a_{t-1} + a_{t-2} <= 2C(d_row - 1, t) =: capref_t, t in [2, m+2]
         (d_row = degree of the row about to be computed at i_fc);
      F4 i_fc + max(ceil(a0/2), max_t T_t) <= R - 2, with EXACT clocks:
         t >= 2: T_t = first K with sum_{k=1..K} (2C(d(i_fc+k)-1,t)
                  - capref_t) >= a_t          (cap-drift staircase);
         t = 1:  T_1 = first K with sum_{k=1..K} max(0, 2(d(i_fc+k)-1)
                  - (1+delta_k) * max(0, a0 - 2(k-1))) >= a_1
                  (exact delta word; drain-assisted).
    Then Theorem S-cone-fc forces capture by the F4 bound; the scan also
    verifies F3-persistence and the actual capture.  All exact."""
    t0 = time.time()
    g = two_G_coeffs(R)
    if from_snap:
        S = json.load(open(os.path.join(
            RESULTS, f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json")))
        j = {int(t): int(v) for t, v in S["junk"].items()}
        i0, d_prev = S["i_postfeed"], S["d_feedend"]
    else:
        d_prev = floor_gamma_star(R) + D0
        j, _, _ = initial_junk(R, d_prev)
        i0 = 1
    i_fc = None
    fc_rec = None
    persist_ok = True
    capref = None
    outcome = None
    capture = None
    for i in range(i0, R - 1):
        d = floor_gamma_star(R + i) + D0
        delta = d - d_prev
        postfeed = (d + i > R)
        if postfeed and j:
            a = {t: -v for t, v in j.items()}
            m = max(a)
            a0 = a.get(0, 0)
            cap = caps_row(d, m + 3)
            F2 = a0 <= d - 1
            F3 = all(2 * a.get(t - 1, 0) + a.get(t - 2, 0) <= cap[t]
                     for t in range(2, m + 3))
            if i_fc is None and F2 and F3 and m + 2 < d and \
                    (i_stop_scan is None or i <= i_stop_scan):
                # exact clocks
                drain = -((-a0) // 2)
                alive = [t for t in range(1, m + 1) if a.get(t, 0)]
                need = {t: a[t] for t in alive}
                cum = {t: 0 for t in alive}
                Tt = {}
                dd = d
                capk = cap[:]      # 2C(dd-1, t)
                Kmax = R - 2 - i
                k = 0
                while need and k < Kmax:
                    k += 1
                    dn = floor_gamma_star(R + i + k) + D0
                    dl = dn - dd
                    if dl:
                        # C(dn-1,t) = C(dd-1,t) * dd/(dd-t)  (dn = dd+1)
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
                        capref = cap[:]
                        fc_rec = {"i_fc": i, "i_pf_offset": None, "d": d,
                                  "m": m, "a0": a0, "drain": drain,
                                  "Tmax": Tmax,
                                  "T_worst_t": max(Tt, key=Tt.get) if Tt
                                  else None,
                                  "capture_ub": ub,
                                  "budget_margin": (R - 2) - ub}
            elif i_fc is not None:
                if not F3:
                    persist_ok = False
        # transport + clamp (engine semantics)
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
        jn = {}
        if w:
            ts = sorted(w)
            ta, tb = ts[0], ts[-1]
            P = comb(d - 1, ta)
            Pprev = comb(d - 1, ta - 1) if ta >= 1 else 0
            for t in range(ta, tb + 1):
                v = w.get(t, 0)
                if v:
                    lo, hi = -2 * P, 2 * Pprev
                    u = min(hi, max(lo, v))
                    if u != v:
                        jn[t] = v - u
                Pprev = P
                P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        if d in jn:
            outcome = {"status": "DIE", "row": i}
            break
        j = jn
        d_prev = d
        if not j and not fed and d + i > R - 1:
            capture = i
            outcome = {"status": "CLOSED", "capture_row": i}
            break
    res = {"R": R, "D0": D0, "outcome": outcome, "capture_row": capture,
           "i_fc": i_fc, "fc": fc_rec, "F3_persist_ok": persist_ok,
           "elapsed_s": round(time.time() - t0, 1)}
    json.dump(res, open(os.path.join(
        RESULTS, f"amm12592_S_cone_fcscan_R{R}_D0{D0}_boxeph.json"), "w"),
        indent=1)
    print(f"[fcscan] R={R} D0={D0}: i_fc={i_fc} rec={fc_rec} "
          f"persist={persist_ok} outcome={outcome} ({res['elapsed_s']}s)",
          flush=True)
    return res


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "certify":
        certify()
    elif mode == "run":
        R, D0 = int(sys.argv[2]), int(sys.argv[3])
        res = flow(R, D0)
        print(f"[run] R={R} D0={D0}: {res['outcome']} "
              f"capture={res['capture_row']} minus2={res['minus2_rows']} "
              f"i_pf={res['i_pf']} entryPASS="
              f"{res['entry'] and res['entry']['ENTRY_PASS']} "
              f"({res['elapsed_s']}s)", flush=True)
    elif mode == "corner":
        R, D0 = int(sys.argv[2]), int(sys.argv[3])
        sb = int(sys.argv[4]) if len(sys.argv) > 4 else 0
        corner(R, D0, scale_bits=sb)
    elif mode == "fromsnap":
        R, D0, sb = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        fromsnap(R, D0, sb)
    elif mode == "iconescan":
        R, D0 = int(sys.argv[2]), int(sys.argv[3])
        iconescan(R, D0, from_snap=("--fromsnap" in sys.argv))
    elif mode == "fcscan":
        R, D0 = int(sys.argv[2]), int(sys.argv[3])
        fcscan(R, D0, from_snap=("--fromsnap" in sys.argv))
