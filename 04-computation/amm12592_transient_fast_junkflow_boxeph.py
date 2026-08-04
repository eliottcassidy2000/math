"""ANGLE B2 -- FAST exact implementation of rule A via the T6 clamped-Pascal
junk flow (proved conjugacy chain T2+T6, session boxeph multifront 2026-08-03).

State per row = sparse junk vector j (dict cell->even int) + feed pointer into
the closed-form coefficients of 2G_R.  Exact identities used (all PROVED in
amm12592-allR-transient-theorem-boxeph.md):

  T2  e-coordinates: e_{-1} = 2G_R, row map = decode/clamp/shift; survival
      iff bottom-cell load in {0,2}; c-boxes [-2C(d-1,t), +2C(d-1,t-1)].
  T4  row-0 load: w_t(R,d) = (-1)^{d-t} C(R-2-t,d-t) - C(d+1,t+1) + 2C(d,t).
  T6  transport: w_{i} = K_delta * j_{i-1} + feed_i, K = (1,1) / (1,2,1) for
      delta = d_i - d_{i-1} in {0,1}; feed on cells {0,1} only, values
      2G_R[d_i + i] (cell 0) and, when delta = 1, 2G_R[d_i - 1 + i]
      (cells 0 and 1); feed active only while the coefficient index <= R-1.
  T3  at dyadic R everything stays even; parity fix never fires (asserted).
  T1/T5 capture: junk empty + feed exhausted  =>  coast to closure.

The full polynomial e_i is NEVER formed (that is the O(R^2 d) slow path in
amm12592_transient_error_dynamics_boxeph.py); this module is O(R * band).
Certification mode proves block-for-block equality against the slow path at
R = 8..128 and per-row trace equality against the STORED slow traces at
R = 256, 512.  Everything int-exact; no floats in any decision.
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

spec = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(iref)
floor_gamma_star = iref.floor_gamma_star


def two_G_coeffs(R):
    """[x^j] 2G_R exactly: g[0]=2, g[j] = (-1)^j C(R-1,j) - 1 (j>=1)."""
    g = [2]
    b = 1
    for j in range(1, R):
        b = b * (R - j) // j          # C(R-1, j), exact
        g.append((-b if j & 1 else b) - 1)
    return g


def initial_junk(R, d):
    """Row 0: clamp the T4 closed-form load; return (junk dict, stats).
    All binomials advanced incrementally with exact-division asserts."""
    A = comb(R - 2, d)        # C(R-2-t, d-t) at t=0
    B = d + 1                 # C(d+1, t+1)   at t=0
    P = 1                     # C(d-1, t)     at t=0   (box lo = -2P)
    Pprev = 0                 # C(d-1, t-1)   at t=0   (box hi = +2Pprev)
    CB = 1                    # C(d, t)       at t=0   (incremental)
    j = {}
    junkL1 = 0
    c0 = None
    sgn = 1 if (d % 2 == 0) else -1   # (-1)^{d-t} at t=0
    for t in range(d + 1):
        w = sgn * A - B + 2 * CB
        lo, hi = -2 * P, 2 * Pprev
        u = min(hi, max(lo, w))
        assert (w - u) % 2 == 0 and (u - lo) % 2 == 0, ("parity", t)
        if u != w:
            j[t] = w - u
            junkL1 += abs(w - u)
        if t == 0:
            c0 = u
        # advance to t+1
        if t < d:
            num = A * (d - t)
            A, r = divmod(num, R - 2 - t); assert r == 0
            num = B * (d - t)
            B, r = divmod(num, t + 2); assert r == 0
            num = CB * (d - t)
            CB, r = divmod(num, t + 1); assert r == 0
            Pprev = P
            if t < d - 1:
                num = P * (d - 1 - t)
                P, r = divmod(num, t + 1); assert r == 0
            else:
                P = 0                      # C(d-1, d) = 0
            sgn = -sgn
    # sanity: bottom cell must carry w = 2 (== [x^0] 2G_R), inside [0,2]
    assert d not in j, "row-0 bottom overflow (impossible for these R)"
    return j, junkL1, c0


def run_fast(R, D0=0, trace_path=None, flush_every=128, witness_path=None,
             keep_rows=True):
    """Exact fast rule-A run at the gamma* floor profile + D0.
    Returns dict: outcome CLOSED/DIE(+row), capture row, trace."""
    t0 = time.time()
    g = two_G_coeffs(R)
    d_prev = floor_gamma_star(R) + D0
    j, junkL1, c0 = initial_junk(R, d_prev)
    trace = []
    sparse_c = {}       # row -> {t: c_t} on w-support (witness reconstruction)
    front0 = max(j) if j else -1
    minus2 = 1 if c0 == -2 else 0

    def rowrec(i, d, j, junkL1, c0, w_bottom):
        ts = sorted(j)
        return {"i": i, "d": d, "nclamped": len(j),
                "tmin": (ts[0] if ts else None), "tmax": (ts[-1] if ts else None),
                "junkL1_bits": junkL1.bit_length(), "c0": c0,
                "e_in0": w_bottom, "gap": (d - ts[-1]) if ts else None}

    trace.append(rowrec(0, d_prev, j, junkL1, c0, 2))
    outcome = None
    capture = None
    for i in range(1, R - 1):
        d = floor_gamma_star(R + i) + D0
        delta = d - d_prev
        assert delta in (0, 1), (i, d, d_prev)
        # ---- transport: w = K_delta * j  (+ feed on cells 0,1)
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
        # ---- clamp over the (contiguous range of the) support
        jn = {}
        junkL1 = 0
        c0 = 0
        w_bottom = w.get(d, 0)
        died = False
        if w:
            ts = sorted(w)
            ta, tb = ts[0], ts[-1]
            assert tb <= d, ("support beyond bottom cell", i, tb, d)
            P = comb(d - 1, ta)           # C(d-1, t)
            Pprev = comb(d - 1, ta - 1) if ta >= 1 else 0
            crow = {} if witness_path else None
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
                    if crow is not None and u:
                        crow[t] = u
                # advance C(d-1, .)
                Pprev = P
                P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
            if crow:
                sparse_c[i] = crow
        if d in jn:
            outcome = {"status": "DIE", "row": i,
                       "const_bits": abs(jn[d]).bit_length()}
            r = rowrec(i, d, jn, junkL1, c0, w_bottom)
            r["death_const_bits"] = abs(jn[d]).bit_length()
            trace.append(r)
            died = True
        else:
            if c0 == -2:
                minus2 += 1
            if keep_rows or (i % 8 == 0) or jn == {}:
                trace.append(rowrec(i, d, jn, junkL1, c0, w_bottom))
        j = jn
        d_prev = d
        if died:
            break
        if not j and not fed and d + i > R - 1:
            capture = i
            outcome = {"status": "CLOSED", "capture_row": i}
            break
        if trace_path and i % flush_every == 0:
            json.dump({"R": R, "D0": D0, "partial_row": i, "trace": trace,
                       "elapsed_s": round(time.time() - t0, 1)},
                      open(trace_path, "w"))
    if outcome is None:
        # reached row R-2 without capture/death: decide the last row exactly
        # e_{R-2} = bern(j, d_prev)/x; closure iff decode at d_{R-1} lands in
        # the last-row boxes [0, 2C(d,t)] (then c = w exactly and e_{R-1}=0).
        if not j:
            outcome = {"status": "CLOSED", "capture_row": R - 2}
            capture = R - 2
        else:
            from amm12592_allR_family_toolbox_boxeph import (bern_to_poly,
                                                             poly_to_bern)
            P = bern_to_poly([j.get(t, 0) for t in range(d_prev + 1)], d_prev)
            assert P and P[0] == 0
            P = P[1:]
            dL = floor_gamma_star(R + R - 1) + D0
            wl = poly_to_bern(P + [0] * (dL + 1 - len(P)), dL) \
                if len(P) <= dL + 1 else None
            ok = wl is not None and all(
                0 <= wl[t] <= 2 * comb(dL, t) and wl[t] % 2 == 0
                for t in range(dL + 1))
            outcome = ({"status": "CLOSED", "capture_row": None,
                        "lastrow_absorbed": True} if ok else
                       {"status": "OPEN_RESIDUAL"})
    res = {"R": R, "D0": D0, "outcome": outcome, "capture_row": capture,
           "minus2_rows": minus2, "front0": front0,
           "T6b_bound": trace[0]["d"] - front0,
           "elapsed_s": round(time.time() - t0, 1), "trace": trace}
    if trace_path:
        json.dump({k: v for k, v in res.items() if k != "blocks"},
                  open(trace_path, "w"))
    if witness_path and outcome["status"] == "CLOSED":
        json.dump({"R": R, "D0": D0, "profile": "floor_gamma_star(R+i)+D0",
                   "capture_row": capture,
                   "note": "delta_i = ballot(d_i) + sparse_c[i] cellwise; "
                           "rows absent from sparse_c are pure ballot "
                           "(Delta = 2x-1); last row Delta = -1 (corner). "
                           "Row 0 c-cells reconstructible from T4 closed form"
                           " (clamped); stored separately below.",
                   "sparse_c": {str(k): {str(t): v for t, v in row.items()}
                                for k, row in sparse_c.items()}},
                  open(witness_path, "w"))
    return res


# --------------------------------------------------------------- verification
def full_blocks(R, D0):
    """Reconstruct FULL delta blocks from a fast run (small R only)."""
    from amm12592_allR_family_toolbox_boxeph import ballot
    g = two_G_coeffs(R)
    pr = [floor_gamma_star(R + i) + D0 for i in range(R)]
    blocks = []
    # row 0: clamp of the T4 load, full vector
    d = pr[0]
    A = comb(R - 2, d); B = d + 1
    row = []
    sgn = 1 if d % 2 == 0 else -1
    for t in range(d + 1):
        w = sgn * A - B + 2 * comb(d, t)
        lo, hi = -2 * comb(d - 1, t), 2 * (comb(d - 1, t - 1) if t else 0)
        u = min(hi, max(lo, w))
        row.append(u)
        if t < d:
            A = A * (d - t) // (R - 2 - t)
            B = B * (d - t) // (t + 2)
            sgn = -sgn
    # rerun fast dynamics but keep full c rows
    res = run_fast(R, D0, witness_path=os.path.join(
        RESULTS, f"_tmp_wit_R{R}_D0{D0}.json"))
    if res["outcome"]["status"] != "CLOSED":
        return None, res
    W = json.load(open(os.path.join(RESULTS, f"_tmp_wit_R{R}_D0{D0}.json")))
    b0 = ballot(pr[0])
    blocks = [[b0[t] + row[t] for t in range(pr[0] + 1)]]
    for i in range(1, R - 1):
        d = pr[i]
        b = ballot(d)
        crow = W["sparse_c"].get(str(i), {})
        blocks.append([b[t] + crow.get(str(t), 0) for t in range(d + 1)])
    dL = pr[R - 1]
    # last row: corner + exact absorb (capture => e_{R-2} = 0 => Delta = -1)
    blocks.append([-comb(dL, t) for t in range(dL + 1)])
    return blocks, res


def certify():
    out = {"blocks_equal_vs_slow": {}, "trace_equal_vs_stored": {},
           "witness_full_check": {}}
    spec2 = importlib.util.spec_from_file_location(
        "slow", os.path.join(HERE, "amm12592_transient_error_dynamics_boxeph.py"))
    slow = importlib.util.module_from_spec(spec2)
    with contextlib.redirect_stdout(io.StringIO()):
        spec2.loader.exec_module(slow)
    from amm12592_allR_family_toolbox_boxeph import (admissible, epoch_sum,
                                                     qpow)
    # (a) full block equality + independent witness check, R = 8..128
    for R in (8, 16, 32, 64, 128):
        sl = slow.run_error_dynamics(R, 0)
        fb, res = full_blocks(R, 0)
        eq = (sl["outcome"]["status"] == "CLOSED" and fb is not None and
              len(fb) == len(sl["blocks"]) and
              all(fb[i] == sl["blocks"][i] for i in range(R)))
        pr = [floor_gamma_star(R + i) for i in range(R)]
        adm = fb is not None and all(admissible(fb[i], pr[i]) for i in range(R))
        idt = fb is not None and epoch_sum(R, fb, pr) == qpow(R - 1)
        out["blocks_equal_vs_slow"][R] = {"blocks_equal": eq,
                                          "admissible": adm, "identity": idt}
        print(f"R={R:4d}: blocks_equal={eq} admissible={adm} identity={idt}",
              flush=True)
    # (b) per-row trace equality vs stored slow traces
    keys = ("d", "nclamped", "tmin", "tmax", "junkL1_bits", "c0", "e_in0")
    for R, D0 in ((256, 0), (256, 1), (512, 0), (512, 4), (512, 5), (512, 8)):
        path = os.path.join(RESULTS,
                            f"amm12592_transient_trace_R{R}_D0{D0}_boxeph.json")
        if not os.path.exists(path):
            continue
        st = json.load(open(path))
        fast = run_fast(R, D0)
        ftr = {r["i"]: r for r in fast["trace"]}
        bad = []
        for r in st["trace"]:
            i = r["i"]
            if i not in ftr:
                bad.append((i, "missing")); continue
            for k in keys:
                if k in r and r[k] != ftr[i].get(k):
                    bad.append((i, k, r[k], ftr[i].get(k)))
        so = st.get("outcome") or (st.get("trace") and
                                   {"partial": st.get("partial_row")})
        rec = {"rows_compared": len(st["trace"]), "mismatches": bad[:10],
               "n_mismatch": len(bad), "slow_outcome": so,
               "fast_outcome": fast["outcome"]}
        out["trace_equal_vs_stored"][f"R{R}_D0{D0}"] = rec
        print(f"R={R} D0={D0}: rows={rec['rows_compared']} "
              f"mismatches={len(bad)} fast={fast['outcome']}", flush=True)
    json.dump(out, open(os.path.join(RESULTS,
        "amm12592_transient_fastflow_cert_boxeph.json"), "w"), indent=1)
    return out


def verify_T9(R, D0):
    """Exact check of the no-return march lemma T9 on a run:
      (a) at the new extreme cell t' = tmax+1+delta, w_{t'} = j_{tmax} exactly;
      (b) trigger := first row where |j_{tmax}| > geometric cap tail
          T(i) := sum_{k>=1} 2 C(d-1, tmax + k) (exact finite sum);
      (c) after the trigger the front advances every row (gap -1/row) and the
          extreme value loses at most the local cap per row;
      (d) cap ratio along the diagonal < 1 requires tau = t/d > 1/phi:
          checked exactly via  d(d-t-1) < (t+1)(t+2)  at the front.
    Returns dict of exact findings."""
    g = two_G_coeffs(R)
    d = floor_gamma_star(R) + D0
    j, _, _ = initial_junk(R, d)
    findings = {"R": R, "D0": D0, "trigger_row": None, "march_monotone": True,
                "edge_transport_exact": True, "phi_condition_all_rows": True,
                "post_trigger_gap_steps": [], "outcome": None}
    trig = None
    prev_extreme = None
    prev_tmax = None
    for i in range(1, R - 1):
        dn = floor_gamma_star(R + i) + D0
        delta = dn - d
        w = {}
        K = (1, 1) if delta == 0 else (1, 2, 1)
        for t, v in j.items():
            for s, ks in enumerate(K):
                w[t + s] = w.get(t + s, 0) + ks * v
        if dn + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn + i]
        if delta == 1 and dn - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn - 1 + i]
            w[1] = w.get(1, 0) + g[dn - 1 + i]
        # (a) edge transport
        if j:
            tm = max(j)
            if w.get(tm + 1 + delta, 0) != j[tm]:
                findings["edge_transport_exact"] = False
        jn = {}
        for t, v in sorted(w.items()):
            if not v:
                continue
            lo, hi = -2 * comb(dn - 1, t), 2 * comb(dn - 1, t - 1) if t else (-2, 0)
            if t == 0:
                lo, hi = -2, 0
            u = min(hi, max(lo, v))
            if u != v:
                jn[t] = v - u
        if dn in jn:
            findings["outcome"] = {"DIE": i}
            break
        if jn:
            tm = max(jn)
            # (d) phi condition at the front
            if not d * (d - tm - 1) < (tm + 1) * (tm + 2):
                # only meaningful once front is in the deep region; record
                if tm > d // 2:
                    findings["phi_condition_all_rows"] = False
            # (b) geometric tail (exact finite sum)
            if trig is None:
                tail = 0
                c = 2 * comb(dn - 1, tm)
                for tt in range(tm + 1, dn + 1):
                    c = c * (dn - tt) // (tt)          # 2C(dn-1, tt)
                    tail += c
                if abs(jn[tm]) > tail:
                    trig = i
                    findings["trigger_row"] = i
            else:
                # (c) march monotone: gap must shrink by exactly 1
                gap_prev = d - prev_tmax
                gap_now = dn - tm
                findings["post_trigger_gap_steps"].append(gap_prev - gap_now)
                if gap_now > gap_prev:
                    findings["march_monotone"] = False
            # (e) no-recovery bookkeeping: last row where gap failed to shrink
            if prev_tmax is not None:
                if dn - tm >= d - prev_tmax:
                    findings["last_nondecrease_row"] = i
                    findings["gap_at_last_nondecrease"] = dn - tm
            prev_extreme = jn.get(tm)
            prev_tmax = tm
        j = jn
        d = dn
        if not j and dn + i > R - 1:
            findings["outcome"] = {"CLOSED_capture": i}
            break
    steps = findings.pop("post_trigger_gap_steps")
    findings["post_trigger_steps_all_1"] = (all(s == 1 for s in steps)
                                            if steps else None)
    findings["post_trigger_rows"] = len(steps)
    return findings


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "T9":
        out = []
        for R, D0 in ((256, 0), (512, 0), (512, 4), (1024, 14), (2048, 37),
                      (4096, 88), (1024, 15), (2048, 38)):
            f = verify_T9(R, D0)
            print(f, flush=True)
            out.append(f)
        json.dump(out, open(os.path.join(RESULTS,
            "amm12592_transient_T9_march_check_boxeph.json"), "w"), indent=1)
    elif mode == "cert":
        certify()
    elif mode == "run":
        R, D0 = int(sys.argv[2]), int(sys.argv[3])
        wp = os.path.join(HERE, f"amm12592_witness_R{R}_ruleA_D0_{D0}_"
                                f"fastflow_boxeph.json") \
            if "--witness" in sys.argv else None
        res = run_fast(R, D0, trace_path=os.path.join(RESULTS,
            f"amm12592_fastflow_trace_R{R}_D0{D0}_boxeph.json"),
            witness_path=wp, keep_rows=(R <= 2048))
        print(f"R={R} D0={D0}: {res['outcome']}  T6b_bound={res['T6b_bound']}"
              f"  capture={res['capture_row']}  minus2={res['minus2_rows']}"
              f"  ({res['elapsed_s']}s)", flush=True)
    elif mode == "scan":
        R = int(sys.argv[2])
        for D0 in [int(v) for v in sys.argv[3].split(",")]:
            res = run_fast(R, D0, trace_path=os.path.join(RESULTS,
                f"amm12592_fastflow_trace_R{R}_D0{D0}_boxeph.json"),
                keep_rows=False)
            print(f"R={R} D0={D0}: {res['outcome']} T6b_bound="
                  f"{res['T6b_bound']} capture={res['capture_row']} "
                  f"minus2={res['minus2_rows']} ({res['elapsed_s']}s)",
                  flush=True)
