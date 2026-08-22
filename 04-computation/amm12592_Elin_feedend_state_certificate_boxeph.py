"""LANE D2 (E-lin) -- exact certificate of the FEED-END STATE of rule A at
linear slack D0 = ceil(eps R), eps in {1/32, 1/16}.

For each (R, D0) the engine's flow is replayed exactly to the last fed row
i_feed (checked against the proved closed form), and the state at the first
autonomous row is certified:

  (N)  all-negative junk:      j_t <= 0 for every cell (and j != {});
  (C)  contiguous low block:   support(j) = [0, m] with no holes;
  (W)  window:                 m <= d/6  (front deep inside the low-cell
                               regime where the A-step lemma propagates);
  (F)  front-freeze row-1:     |j_m| <= 2C(d, m+2)  and
                               2|j_m| + |j_{m-1}| <= 2C(d, m+1)
                               (the exact conditions under which the next
                               row cannot advance the front; re-established
                               at later rows by the run itself);
  (D)  debt identity:          j_0 = (2 - R) + 2 * #{c0 = -2 rows so far}
                               (T8 + the x=1 Bernstein evaluation: post-feed
                               s_i = e_i(1) = j_0);
  (DL) drain deadline:         |j_0| <= 2 (R - 2 - i_feed)
                               (capture consistent with 2/row drain);
  (A)  cap-ratio profile:      bits of A_t = |j_t| / (2 C(d-1,t)) recorded
                               per cell (the envelope any future proof of
                               Hypothesis S must carry).

All decisions exact; ratios reported as bit-lengths only.
"""
import os, sys, json, importlib.util, contextlib, io
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
OUT = os.path.join(RESULTS, "amm12592_Elin_feedend_state_boxeph.json")

spec = importlib.util.spec_from_file_location(
    "fast", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fast = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fast)
fg = fast.floor_gamma_star


def ceil_div(a, b):
    return -((-a) // b)


def feedend_state(R, D0):
    g = fast.two_G_coeffs(R)
    d = fg(R) + D0
    j, _, c0 = fast.initial_junk(R, d)
    minus2 = 1 if c0 == -2 else 0
    i_feed_obs = None
    for i in range(1, R - 1):
        dn = fg(R + i) + D0
        delta = dn - d
        w = {}
        K = (1, 1) if delta == 0 else (1, 2, 1)
        for t, v in j.items():
            for s, ks in enumerate(K):
                w[t + s] = w.get(t + s, 0) + ks * v
        fed = False
        if dn + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn + i]; fed = True
        if delta == 1 and dn - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn - 1 + i]
            w[1] = w.get(1, 0) + g[dn - 1 + i]; fed = True
        jn = {}
        c0 = 0
        if w:
            ts_w = sorted(w)
            ta, tb = ts_w[0], ts_w[-1]
            P = comb(dn - 1, ta)
            Pprev = comb(dn - 1, ta - 1) if ta >= 1 else 0
            for t in range(ta, tb + 1):
                v = w.get(t, 0)
                if v:
                    lo, hi = -2 * P, 2 * Pprev
                    u = min(hi, max(lo, v))
                    if t == 0:
                        c0 = u
                    if u != v:
                        jn[t] = v - u
                Pprev = P
                P = P * (dn - 1 - t) // (t + 1) if t < dn - 1 else 0
        assert dn not in jn, ("death during certified feed phase!", R, D0, i)
        if c0 == -2:
            minus2 += 1
        j, d = jn, dn
        if not fed:
            i_feed_obs = i - 1
            break
    # state at first autonomous row (row i = i_feed + 1 already computed:
    # the loop broke at the first row with fed = False, i.e. j is the junk
    # AFTER row i_feed + 1 ... re-read: fed=False means THIS row had no
    # feed, so i_feed = i - 1 and j = state after the first autonomous row.
    ts = sorted(j)
    m = ts[-1] if ts else None
    rec = {"R": R, "D0": D0, "d": d, "i_feed_obs": i_feed_obs,
           "first_autonomous_row": i_feed_obs + 1,
           "minus2_rows_so_far": minus2}
    rec["N_all_negative"] = bool(ts) and all(j[t] < 0 for t in ts)
    rec["C_contiguous_0_to_m"] = (ts == list(range(m + 1))) if ts else False
    rec["m"] = m
    rec["W_m_le_d_over_6"] = (6 * m <= d) if m is not None else None
    if m is not None and m >= 1:
        f1 = abs(j[m]) <= 2 * comb(d, m + 2)
        f2 = 2 * abs(j[m]) + abs(j[m - 1]) <= 2 * comb(d, m + 1)
        rec["F_front_freeze_row1"] = bool(f1 and f2)
    debt_id = (j.get(0, 0) == (2 - R) + 2 * minus2)
    rec["D_debt_identity"] = bool(debt_id)
    rec["j0"] = str(j.get(0, 0))
    rec["DL_drain_deadline"] = bool(abs(j.get(0, 0))
                                    <= 2 * (R - 2 - i_feed_obs))
    rec["A_bits_profile"] = [
        (t, abs(j[t]).bit_length() - (2 * comb(d - 1, t)).bit_length())
        for t in ts]        # ~log2 A_t (within 1); front cells last
    rec["junkL1_bits"] = sum(abs(v) for v in j.values()).bit_length()
    return rec


if __name__ == "__main__":
    Rs = [128, 256, 512, 1024, 2048, 4096, 8192]
    if "--big" in sys.argv:
        Rs.append(16384)
    out = []
    if os.path.exists(OUT):
        out = json.load(open(OUT))
    done = {(r["R"], r["D0"]) for r in out}
    for R in Rs:
        for D0 in (ceil_div(R, 32), ceil_div(R, 16)):
            if (R, D0) in done:
                continue
            rec = feedend_state(R, D0)
            out.append(rec)
            json.dump(out, open(OUT, "w"), indent=1)
            ok = all(rec[k] for k in
                     ("N_all_negative", "C_contiguous_0_to_m",
                      "W_m_le_d_over_6", "F_front_freeze_row1",
                      "D_debt_identity", "DL_drain_deadline"))
            print(f"R={R:6d} D0={D0:5d}: i_feed={rec['i_feed_obs']} "
                  f"m={rec['m']} j0={rec['j0']} minus2={rec['minus2_rows_so_far']} "
                  f"ALL_OK={ok}", flush=True)
    allok = all(all(r[k] for k in ("N_all_negative", "C_contiguous_0_to_m",
                                   "W_m_le_d_over_6", "F_front_freeze_row1",
                                   "D_debt_identity", "DL_drain_deadline"))
                for r in out)
    print("FEEDEND_ALL_OK =", allok, "->", OUT, flush=True)
