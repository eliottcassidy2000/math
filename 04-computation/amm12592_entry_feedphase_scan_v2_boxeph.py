"""LANE F1 -- ENTRY SCAN v2: same exact semantics as
amm12592_entry_feedphase_scan_boxeph.py (v1, itself validated bit-exactly
against the certified E1 fcscan/feedend records), re-implemented for speed:

  * junk/load/caps as contiguous LISTS (support of the T4b junk block is
    [0, m]; interior zeros are carried harmlessly);
  * caps updated by PASCAL'S RULE: 2C(d,t) = 2C(d-1,t) + 2C(d-1,t-1),
    one big-int ADD per cell on delta = 1 rows, nothing on delta = 0 rows
    (v1 recomputed C(d-1,t) by mul//div per cell per row);
  * the exact clock staircase also advances its cap row by Pascal;
  * parity asserts on every 128th row (full row) instead of every cell
    (T3 is PROVED; the strided assert remains as a tripwire);
  * identical outputs, identical ledger schema, same file names.

Validation mode compares v2 against the v1 ledgers at (1024, 32),
(2048, 64), (4096, 128): every reported field must be equal (i_pf, d_fe,
feed-end junk dict, lastpos/firstF2/firstF3 rows, full fc record).
All decisions exact int; no floats; no numpy; no gmpy2.
"""
import sys, os, json, time
from math import comb, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

import importlib.util as _ilu
_spec = _ilu.spec_from_file_location(
    "escan1", os.path.join(HERE, "amm12592_entry_feedphase_scan_boxeph.py"))
_v1 = _ilu.module_from_spec(_spec)
import io as _io, contextlib as _ctx
with _ctx.redirect_stdout(_io.StringIO()):
    _spec.loader.exec_module(_v1)
initial_junk = _v1.initial_junk
floor_gamma_star = _v1.floor_gamma_star
two_G_coeffs_top = _v1.two_G_coeffs_top


def scan(R, D0, persist_rows=64, stride=256, postfeed_diag_rows=12,
         max_scan_past_pf=None, flush_every=2048):
    t0 = time.time()
    d_prev = floor_gamma_star(R) + D0
    d0 = d_prev
    g = two_G_coeffs_top(R, d0 - 2)
    glow = max(1, d0 - 2)
    jd, junkL1, c0 = initial_junk(R, d_prev)
    mj = max(jd)
    j = [jd.get(t, 0) for t in range(mj + 1)]
    del jd
    minus2 = 1 if c0 == -2 else 0
    if max_scan_past_pf is None:
        max_scan_past_pf = R // 8
    # caps[t] = 2*C(dm1, t), dm1 tracks d_i - 1 of the row being clamped
    dm1 = d_prev - 1
    caps = [2]
    P = 1
    need_len = mj + 8
    for t in range(need_len):
        P = P * (dm1 - t) // (t + 1) if t < dm1 else 0
        caps.append(2 * P)
    i_pf = None
    d_fe = None
    fe_rec = None
    lastpos_row = 0 if any(v > 0 for v in j) else None
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

    def cur_m():
        for t in range(len(j) - 1, -1, -1):
            if j[t]:
                return t
        return -1

    def sample(i, d):
        m = cur_m()
        npos = sum(1 for v in j if v > 0)
        L1 = 0
        mb = 0
        for v in j:
            if v:
                L1 += v if v > 0 else -v
                bl = (v if v > 0 else -v).bit_length()
                if bl > mb:
                    mb = bl
        traj.append({"i": i, "d": d, "m": m, "ncells": sum(1 for v in j if v),
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
        # ---- advance caps to degree d-1 (Pascal, descending, in place)
        if delta:
            dm1 += 1
            for t in range(len(caps) - 1, 0, -1):
                caps[t] += caps[t - 1]
        postfeed = (d + i > R)
        if postfeed and i_pf is None:
            i_pf = i
            d_fe = d_prev
            g = None
            m = cur_m()
            a0 = -j[0] if j and j[0] < 0 else 0
            npos = sum(1 for v in j if v > 0)
            assert a0 == (R - 2) - 2 * minus2, ("T8 identity", i, a0, minus2)
            fe_rec = {"i_pf": i_pf, "d_fe": d_fe, "m": m, "a0": a0,
                      "a0_minus_dfe": a0 - d_fe,
                      "ncells": sum(1 for v in j if v),
                      "npos": npos, "allneg": npos == 0,
                      "contiguous": all(j[t] for t in range(m + 1)),
                      "minus2": minus2,
                      "m_x1000_over_sqrtR": (m * 1000) // isqrt(R * 1000000)
                      if m >= 0 else None,
                      "junk": {str(t): str(j[t]) for t in range(m + 1)
                               if j[t]} if 0 <= m <= 400 else None}
        if postfeed and any(j):
            m = cur_m()
            a0 = -j[0] if j[0] < 0 else 0
            while len(caps) < m + 4:      # extend at degree dm1
                tt = len(caps) - 1        # caps[tt] known; add caps[tt+1]
                caps.append(caps[tt] * (dm1 - tt) // (tt + 1)
                            if tt < dm1 else 0)
            def A(t):
                return -j[t] if 0 <= t <= m and j[t] < 0 else 0
            F2 = a0 <= d - 1
            F3 = all(2 * A(t - 1) + A(t - 2) <= caps[t]
                     for t in range(2, m + 3))
            if F2 and first_F2 is None:
                first_F2 = i
            if F3 and first_F3 is None:
                first_F3 = i
            if len(pf_diags) < postfeed_diag_rows or (i_fc is not None and
                                                     i <= i_fc + 2):
                worst_num, worst_den, worst_t = 0, 1, None
                over = []
                min_marg_num, min_marg_den, min_t = None, None, None
                for t in range(2, m + 3):
                    sp = 2 * A(t - 1) + A(t - 2)
                    if sp == 0 or caps[t] == 0:
                        continue
                    if sp * worst_den > worst_num * caps[t]:
                        worst_num, worst_den, worst_t = sp, caps[t], t
                    if sp > caps[t]:
                        over.append(t)
                    marg_n, marg_d = caps[t] - sp, caps[t]
                    if (min_marg_num is None or
                            marg_n * min_marg_den < min_marg_num * marg_d):
                        min_marg_num, min_marg_den, min_t = marg_n, marg_d, t
                pf_diags.append(
                    {"i": i, "d": d, "delta": delta, "m": m, "a0": a0,
                     "a0_vs_dminus1": (d - 1) - a0,
                     "npos": sum(1 for v in j if v > 0),
                     "F2": F2, "F3": F3,
                     "max_r_float": (worst_num / worst_den
                                     if worst_den else None),
                     "max_r_t": worst_t, "over_cells": over,
                     "min_F3_margin_float": (min_marg_num / min_marg_den
                                             if min_marg_den else None),
                     "min_F3_margin_t": min_t})
            if i_fc is None and F2 and F3 and m + 2 < d:
                drain = -((-a0) // 2)
                alive = [t for t in range(1, m + 1) if A(t)]
                need = {t: A(t) for t in alive}
                cum = {t: 0 for t in alive}
                Tt = {}
                dd = d
                capk = caps[:m + 3]
                capref = caps[:m + 3]
                Kmax = R - 2 - i
                k = 0
                while need and k < Kmax:
                    k += 1
                    dn = floor_gamma_star(R + i + k) + D0
                    dl = dn - dd
                    if dl:
                        capk = capk[:]
                        for t in range(len(capk) - 1, 0, -1):
                            capk[t] += capk[t - 1]
                    a0k = max(0, a0 - 2 * (k - 1))
                    for t in list(need):
                        if t == 1:
                            dec = 2 * (dn - 1) - (1 + dl) * a0k
                            if dec > 0:
                                cum[1] += dec
                        else:
                            dec = capk[t] - capref[t]
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
                    persist_ok = False
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
        # ---- transport: w = K_delta * j (+ feed), list form
        mj = len(j) - 1
        if delta == 0:
            w = j + [0]
            for t in range(mj, -1, -1):
                v = j[t]
                if v:
                    w[t + 1] += v
        else:
            w = j + [0, 0]
            for t in range(mj, -1, -1):
                v = j[t]
                if v:
                    w[t + 1] += 2 * v
                    w[t + 2] += v
        fed = False
        if d + i <= R - 1:
            w[0] += g[d + i]; fed = True
        if delta == 1 and d - 1 + i <= R - 1:
            w[0] += g[d - 1 + i]
            w[1] += g[d - 1 + i]; fed = True
        if g is not None:
            while glow <= d + i - 2:
                g.pop(glow, None)
                glow += 1
        # ---- clamp, list form
        while len(caps) < len(w) + 1:
            tt = len(caps) - 1
            caps.append(caps[tt] * (dm1 - tt) // (tt + 1) if tt < dm1 else 0)
        assert len(w) - 1 <= d, ("support beyond bottom cell", i, len(w), d)
        rowcheck = (i & 127) == 0
        jn = [0] * len(w)
        haspos = False
        newm = -1
        c0 = 0
        died = False
        for t in range(len(w)):
            v = w[t]
            if v:
                if t == 0:
                    lo, hi = -2, 0
                else:
                    lo, hi = -caps[t], caps[t - 1]
                if v > hi:
                    q = v - hi
                elif v < lo:
                    q = v - lo
                else:
                    q = 0
                if rowcheck:
                    assert v % 2 == 0 and q % 2 == 0, ("parity", i, t)
                if q:
                    jn[t] = q
                    if q > 0:
                        haspos = True
                    if t > newm:
                        newm = t
                if t == 0:
                    c0 = v - q
        if newm == d:
            outcome = {"status": "DIE", "row": i,
                       "const_bits": abs(jn[d]).bit_length()}
            died = True
        elif c0 == -2:
            minus2 += 1
        del j[:]
        j = jn[:newm + 1] if newm >= 0 else []
        keep = max(newm + 8, 16)
        if len(caps) > keep + 64:      # trim receded front (re-extend on
            del caps[keep:]            # demand by exact tail recurrence)
        d_prev = d
        if haspos:
            lastpos_row = i
        if i % stride == 0:
            sample(i, d)
        if died:
            flush(final=True)
            break
        if newm < 0 and not fed and d + i > R - 1:
            outcome = {"status": "CLOSED_EARLY", "capture_row": i}
            flush(final=True)
            break
        if i % flush_every == 0:
            flush()
    if outcome is None:
        outcome = {"status": "SCAN_END"}
        flush(final=True)
    print(f"[entry-v2] R={R} D0={D0}: i_pf={i_pf} lastpos={lastpos_row} "
          f"firstF2={first_F2} firstF3={first_F3} i_fc={i_fc} fc={fc_rec} "
          f"persist_ok={persist_ok}({persist_seen}) outcome={outcome} "
          f"({round(time.time()-t0,1)}s)", flush=True)
    return {"R": R, "D0": D0, "i_pf": i_pf, "fe": fe_rec, "fc": fc_rec,
            "lastpos_row": lastpos_row, "first_F2": first_F2,
            "first_F3": first_F3, "outcome": outcome}


def validate():
    """v2 vs the v1 ledgers (v1 itself validated vs the E1 records)."""
    ok_all = True
    for R, D0 in ((1024, 32), (2048, 64), (4096, 128)):
        p = os.path.join(RESULTS,
                         f"amm12592_entry_scan_R{R}_D0{D0}_boxeph.json")
        v1led = json.load(open(p))
        res = scan(R, D0, persist_rows=v1led["F3_persist_rows_checked"])
        v2led = json.load(open(p))
        keys = ("i_pf", "d_feedend", "feedend", "lastpos_row",
                "first_F2_row", "first_F3_row", "i_fc", "fc",
                "F3_persist_ok")
        checks = {k: v1led[k] == v2led[k] for k in keys}
        ok = all(checks.values())
        ok_all = ok_all and ok
        print(f"[validate-v2] R={R} D0={D0}: "
              f"{'ALL-PASS' if ok else 'FAIL ' + str([k for k, v in checks.items() if not v])}",
              flush=True)
    print(f"[validate-v2] overall: {'ALL-PASS' if ok_all else 'FAIL'}",
          flush=True)
    return ok_all


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        validate()
    else:
        R, D0 = int(sys.argv[1]), int(sys.argv[2])
        scan(R, D0)
