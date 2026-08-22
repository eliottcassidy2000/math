"""Lane E2 -- Hypothesis S via super-block impossibility: the post-feed
capture CERTIFICATE.

Input: an exact feed-end state (output of amm12592_S_feedend_extract_boxeph.py)
= the junk vector j_{i0} of plain rule A at the first autonomous row i0, at
degree d_{i0}, for a dyadic epoch R at slack D0.

The certificate machine-checks the hypotheses and the inequality chain of the
theorems in 05-knowledge/results/amm12592-S-superblock-proof-boxeph.md:

  T1 (march budget / late-death, PROVED unconditionally):
      death at row i requires  i >= d_{i0} + i0 - m  >=  R - m,
      so death can only occur in the last ~m rows and requires the junk front
      to advance at the maximal T6a rate on all but at most m-1 rows of the
      entire post-feed phase (advance deficit <= m-1).
  LM (one-step majorant, PROVED): the extended C-A cap-ratio inequalities,
      including the exact cell-0 debt drain (j0' = j0 + 2 per row, from
      w_0 = j_0, cap 2) and the debt pump into cells 1 and 2
      (|j0|/d resp. |j0|/(d(d-1)) in cap-ratio units) -- the "T7 block
      version" of the drain.
  LF (front freeze in cap-ratio form, PROVED = C-F rewritten):
      delta=1: A_m <= d(d-m-1)/((m+1)(m+2))  and
               2A_m + A_{m-1} m/(d-m) <= d/(m+1);
      delta=0: A_m <= (d-1-m)/(m+1).
  T2 (certificate soundness, PROVED): if the fixed-point majorant iteration
      B (dyadic upper bounds, ceilings only -- every rounding is upward, so
      B dominates every true cap-ratio A at every row) keeps the freeze
      inequalities valid at every row until (cells 1..m extinct AND debt
      drained), and this happens within the row budget R-2-i0, then the
      autonomous flow of rule A captures (junk -> empty) by row R-2 with no
      death: Hypothesis S holds at this (R, D0).
  SB (super-block impossibility, PROVED corollary): under the certificate the
      junk values are <= sup_k B_t(k) * 2C(d_max-1, t) while the cap tail
      beyond the front always contains the middle binomial 2C(d-1,(d-1)/2);
      the certified value/tail gap shows no super block (D1 sec 3b sense)
      can EVER form post-feed.

All decisions exact integer arithmetic (fixed-point ceilings at P = 96 bits,
always rounded toward the safe side); floats appear only in *reporting*.
Self-test mode replays the TRUE autonomous flow at small R and checks
domination of every cell at every row plus the capture-row bound.

Usage:
  python3 amm12592_S_superblock_certificate_boxeph.py cert R D0
  python3 amm12592_S_superblock_certificate_boxeph.py selftest R D0
  python3 amm12592_S_superblock_certificate_boxeph.py negctrl R D0
  python3 amm12592_S_superblock_certificate_boxeph.py all       (ledger run)
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

spec = importlib.util.spec_from_file_location(
    "fastflow", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fastflow = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fastflow)
floor_gamma_star = fastflow.floor_gamma_star

P = 96                    # fixed-point bits for the majorant
S = 1 << P

# D2 sweep capture rows (exact, verified): consistency oracle only.
KNOWN_CAPTURE = {(128, 4): 79, (128, 8): 78, (256, 8): 153, (256, 16): 150,
                 (512, 16): 317, (512, 32): 313, (1024, 32): 624,
                 (1024, 64): 616, (2048, 64): 1261, (2048, 128): 1237,
                 (4096, 128): 2519, (4096, 256): 2486, (8192, 256): 5040,
                 (8192, 512): 4964, (16384, 512): 10090, (16384, 1024): 9937}


def ceildiv(a, b):
    assert a >= 0 and b > 0
    return -((-a) // b)


def flog2(a, b):
    """floor(log2(a/b)) for positive ints (exact)."""
    assert a > 0 and b > 0
    e = a.bit_length() - b.bit_length()
    if e >= 0:
        return e if a >= (b << e) else e - 1
    return e if (a << (-e)) >= b else e - 1


def slack_bits(lhs, rhs):
    """log2(rhs/lhs) rounded down; None-safe for lhs = 0 (infinite slack)."""
    if lhs == 0:
        return None
    return flog2(rhs, lhs)


def load_state(R, D0):
    p = os.path.join(RESULTS, "amm12592_S_feedend_R%d_D0%d_boxeph.json" % (R, D0))
    st = json.load(open(p))
    st["junk"] = {int(t): int(v) for t, v in st["junk"].items()}
    st["j0"] = int(st["j0"])
    return st


def published_record(R, D0):
    p = os.path.join(RESULTS, "amm12592_Elin_feedend_state_boxeph.json")
    if not os.path.exists(p):
        return None
    for r in json.load(open(p)):
        if r["R"] == R and r["D0"] == D0:
            return r
    return None


def certify(R, D0, verbose=True):
    t0 = time.time()
    st = load_state(R, D0)
    i0, d0c, m = st["first_autonomous_row"], st["d"], st["m"]
    j = st["junk"]
    rep = {"R": R, "D0": D0, "i0": i0, "d_i0": d0c, "m": m, "checks": {}}
    ck = rep["checks"]

    # ---------------- hypotheses on the state (exact) --------------------
    ck["H_all_negative"] = all(v < 0 for v in j.values())
    ck["H_all_even"] = all(v % 2 == 0 for v in j.values())
    ck["H_support"] = (min(j) >= 0 and max(j) == m)
    ck["H_W_6m_le_d"] = (6 * m <= d0c)
    ck["H_m_ge_3"] = (m >= 3)
    ck["H_debt_negative"] = (j.get(0, 0) < 0)
    j0abs = abs(j.get(0, 0))
    rep["j0abs"] = j0abs
    rep["debt_edge_c0"] = (j0abs - d0c)          # observed edge law |j0|-d
    # regression vs the published D2 feed-end record (if present)
    pub = published_record(R, D0)
    if pub is not None:
        ok = (pub["first_autonomous_row"] == i0 and pub["d"] == d0c
              and pub["m"] == m and int(pub["j0"]) == j.get(0, 0)
              and pub["junkL1_bits"] == st["junkL1_bits"]
              and pub["i_feed_obs"] == st["i_feed_obs"])
        # published A-bits = bit_length(|j_t|) - bit_length(2C(d-1,t))
        prof_ok = True
        caps0 = [2 * comb(d0c - 1, t) for t in range(m + 1)]
        Aprof = {t: (-j[t]).bit_length() - caps0[t].bit_length()
                 for t in j if -j[t] > 0}
        for t, bits in pub.get("A_bits_profile", []):
            if Aprof.get(t) != bits:
                prof_ok = False
        ck["REGRESSION_vs_published_record"] = ok and prof_ok
    else:
        ck["REGRESSION_vs_published_record"] = "no published record (new scale)"

    # ---------------- Theorem 1: march budget / late death ----------------
    # autonomy for good: d_{i0} + i0 >= R.  Then for every i > i0 (strict
    # increase of d_i + i) both feed branches (d+i <= R-1; d+i <= R with
    # delta = 1) are dead forever.
    ck["T1_autonomy"] = (d0c + i0 >= R)
    i_death_min = d0c + i0 - m
    rep["i_death_min"] = i_death_min
    rep["death_window_rows"] = max(0, (R - 1) - i_death_min + 1)
    rep["deficit_budget"] = (R - 1) - i_death_min   # frozen-row budget, death
    ck["T1_death_only_last_rows"] = (i_death_min > R - m - 3)

    # ---------------- majorant iteration --------------------------------
    K = (R - 2) - i0                     # available transitions
    caps0 = [2 * comb(d0c - 1, t) for t in range(m + 1)]
    N = [0] * (m + 1)                    # fixed point: A_t <= N_t / S, t>=1
    for t in range(1, m + 1):
        a = abs(j.get(t, 0))             # |j_t|; sign hypothesis gated by H_
        if a:
            N[t] = ceildiv(a * S, caps0[t])
    supN = list(N)
    k_debt = j0abs // 2                  # exact rows to drain the debt
    ck["DL_debt_within_budget"] = (k_debt <= K)

    # degree sequence, computed lazily in chunks
    D = [d0c]
    def deg(idx):
        while len(D) <= idx:
            i = i0 + len(D)
            D.append(floor_gamma_star(R + i) + D0)
            assert D[-1] - D[-2] in (0, 1)
        return D[idx]

    freeze_ok = True
    first_freeze_fail = None
    min_fs = [None, None, None]          # freeze slacks (bits): d1i, d1ii, d0
    k = 0
    k_last_junk = 0                      # last k with some N_t > 0
    while True:
        j0k = max(0, j0abs - 2 * k)
        if all(n == 0 for n in N) and j0k == 0:
            break                        # captured (majorant): row i0 + k
        if k >= K:
            break                        # out of budget -> FAIL below
        d = deg(k); dn = deg(k + 1); delta = dn - d
        # ---- freeze checks (support cannot pass m on this transition)
        if delta == 1:
            l1, r1 = N[m] * (m + 1) * (m + 2), S * d * (d - m - 1)
            l2, r2 = (2 * N[m] * (d - m) * (m + 1) + N[m - 1] * m * (m + 1),
                      S * d * (d - m))
            okf = (l1 <= r1) and (l2 <= r2)
            for idx, (l, r) in enumerate(((l1, r1), (l2, r2))):
                sb = slack_bits(l, r)
                if sb is not None and (min_fs[idx] is None or sb < min_fs[idx]):
                    min_fs[idx] = sb
        else:
            l3, r3 = N[m] * (m + 1), S * (d - 1 - m)
            okf = (l3 <= r3)
            sb = slack_bits(l3, r3)
            if sb is not None and (min_fs[2] is None or sb < min_fs[2]):
                min_fs[2] = sb
        if not okf and freeze_ok:
            freeze_ok = False
            first_freeze_fail = k
        # ---- one majorant step (every rounding upward)
        Nn = [0] * (m + 1)
        if delta == 1:
            for t in range(1, m + 1):
                v = ceildiv(N[t] * (d - t), d)
                if t == 1:
                    v += ceildiv(S * j0k, d)
                else:
                    v += ceildiv(2 * N[t - 1] * t, d)
                    if t == 2:
                        v += ceildiv(S * j0k, d * (d - 1))
                    else:
                        v += ceildiv(N[t - 2] * t * (t - 1),
                                     d * (d - t + 1))
                Nn[t] = max(0, v - S)
        else:
            for t in range(1, m + 1):
                v = N[t]
                if t == 1:
                    v += ceildiv(S * j0k, 2 * (d - 1))
                else:
                    v += ceildiv(N[t - 1] * t, d - t)
                Nn[t] = max(0, v - S)
        N = Nn
        k += 1
        for t in range(1, m + 1):
            if N[t] > supN[t]:
                supN[t] = N[t]
            if N[t]:
                k_last_junk = k
        # runaway guard (majorant divergence => FAIL, don't burn cycles)
        if max(N) > (S << 40):
            break

    captured = all(n == 0 for n in N) and (j0abs - 2 * k <= 0)
    k_end = k
    ck["CAPTURE_within_budget"] = captured and (k_end <= K)
    ck["FREEZE_all_rows"] = freeze_ok
    rep["k_end"] = k_end
    rep["k_debt"] = k_debt
    rep["k_junk_stable"] = k_last_junk   # cells 1..m extinct (and stay) after
    rep["budget_K"] = K
    rep["capture_row_bound"] = i0 + k_end
    rep["freeze_first_fail_k"] = first_freeze_fail
    rep["freeze_min_slack_bits"] = {"d1_extreme": min_fs[0],
                                    "d1_mate": min_fs[1], "d0": min_fs[2]}
    # exact capture prediction: debt is exactly 2/row (LM), so if the interior
    # is (majorant-)extinct by k_debt the true capture row IS i0 + k_debt.
    if captured and k_last_junk <= k_debt:
        rep["capture_row_exact_pred"] = i0 + k_debt
    else:
        rep["capture_row_exact_pred"] = None
    act = KNOWN_CAPTURE.get((R, D0))
    if act is not None:
        rep["capture_row_actual"] = act
        ck["CONSISTENCY_bound_ge_actual"] = (i0 + k_end >= act)
        if rep["capture_row_exact_pred"] is not None:
            ck["CONSISTENCY_exact_prediction"] = \
                (rep["capture_row_exact_pred"] == act)

    # ---------------- super-block impossibility margin -------------------
    d_max = deg(k_end)
    vbits = j0abs.bit_length()           # cell-0 value (debt)
    for t in range(1, m + 1):
        if supN[t]:
            vb = flog2(supN[t] * 2 * comb(d_max - 1, t), S) + 1
            vbits = max(vbits, vb)
    mid = (d0c - 1) // 2
    tail_bits = flog2(2 * comb(d0c - 1, mid), 1)
    rep["superblock_value_bits_max"] = vbits
    rep["superblock_tailmin_bits"] = tail_bits
    ck["SB_no_superblock_ever"] = (vbits < tail_bits)
    rep["superblock_margin_bits"] = tail_bits - vbits

    hyp = all(v is True for kk, v in ck.items() if kk.startswith("H_"))
    rep["PASS"] = (hyp and ck["T1_autonomy"] and ck["FREEZE_all_rows"]
                   and ck["CAPTURE_within_budget"]
                   and ck["DL_debt_within_budget"]
                   and ck["SB_no_superblock_ever"])
    rep["elapsed_s"] = round(time.time() - t0, 1)
    if verbose:
        tag = "PASS" if rep["PASS"] else "FAIL"
        print("[%s] R=%d D0=%d  i0=%d d=%d m=%d |j0|=%d (=d%+d)  "
              "k_end=%d (junk %d, debt %d, budget %d)  capture<=%d%s  "
              "freeze slacks %s  SB margin %d bits  (%.1fs)"
              % (tag, R, D0, i0, d0c, m, j0abs, rep["debt_edge_c0"],
                 k_end, k_last_junk, k_debt, K, i0 + k_end,
                 ("" if act is None else " (actual %d, pred %s)"
                  % (act, rep["capture_row_exact_pred"])),
                 rep["freeze_min_slack_bits"], rep["superblock_margin_bits"],
                 rep["elapsed_s"]))
        for kk, v in ck.items():
            if v is not True and not isinstance(v, str):
                print("    CHECK FAILED: %s = %s" % (kk, v))
    return rep


# ------------------------------------------------------------ self-test
def selftest(R, D0, verbose=True):
    """Hostile check of the majorant: replay the TRUE autonomous flow from the
    state and verify domination of every cell at every row, plus support/sign
    invariants and the capture-row bound.  Exact integer comparisons only."""
    st = load_state(R, D0)
    i0, d0c, m = st["first_autonomous_row"], st["d"], st["m"]
    j = dict(st["junk"])
    j0abs = -j.get(0, 0)
    # rebuild the majorant trajectory, recording N per row
    rows = []
    rep = certify(R, D0, verbose=False)
    assert rep["PASS"], "certificate itself failed; selftest moot"
    caps0 = [2 * comb(d0c - 1, t) for t in range(m + 1)]
    N = [0] * (m + 1)
    for t in range(1, m + 1):
        a = -j.get(t, 0)
        if a:
            N[t] = ceildiv(a * S, caps0[t])
    d_prev = d0c
    k = 0
    ok = True
    while j:
        # domination check at row i0+k
        d = d_prev
        for t in range(0, m + 1):
            a = -j.get(t, 0)
            if t == 0:
                if a != max(0, j0abs - 2 * k):
                    ok = False; print("  drain mismatch at k=%d" % k)
            elif a * S > N[t] * 2 * comb(d - 1, t):
                ok = False
                print("  DOMINATION FAIL k=%d t=%d" % (k, t))
        if max(j) > m or any(v > 0 for v in j.values()):
            ok = False
            print("  support/sign FAIL at k=%d" % k)
        # true one step
        i = i0 + k + 1
        dn = floor_gamma_star(R + i) + D0
        delta = dn - d_prev
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
        assert dn + i > R - 1, "feed would fire post-feed?!"
        assert not (delta == 1 and dn - 1 + i <= R - 1)
        jn = {}
        for t, v in sorted(w.items()):
            lo, hi = -2 * comb(dn - 1, t), 2 * comb(dn - 1, t - 1) if t else 0
            if t == 0:
                lo, hi = -2, 0
            u = min(hi, max(lo, v))
            if u != v:
                jn[t] = v - u
        assert dn not in jn, "DEATH in true flow (certificate unsound!)"
        # majorant one step (same code path as certify)
        j0k = max(0, j0abs - 2 * k)
        d = d_prev
        Nn = [0] * (m + 1)
        if delta == 1:
            for t in range(1, m + 1):
                v = ceildiv(N[t] * (d - t), d)
                if t == 1:
                    v += ceildiv(S * j0k, d)
                else:
                    v += ceildiv(2 * N[t - 1] * t, d)
                    if t == 2:
                        v += ceildiv(S * j0k, d * (d - 1))
                    else:
                        v += ceildiv(N[t - 2] * t * (t - 1), d * (d - t + 1))
                Nn[t] = max(0, v - S)
        else:
            for t in range(1, m + 1):
                v = N[t]
                if t == 1:
                    v += ceildiv(S * j0k, 2 * (d - 1))
                else:
                    v += ceildiv(N[t - 1] * t, d - t)
                Nn[t] = max(0, v - S)
        N = Nn
        j = jn
        d_prev = dn
        k += 1
        if k > rep["budget_K"] + 2:
            ok = False; print("  true flow exceeded budget?!"); break
    cap_row = i0 + k
    if cap_row > rep["capture_row_bound"]:
        ok = False
        print("  capture bound %d < true capture %d"
              % (rep["capture_row_bound"], cap_row))
    if verbose:
        print("[%s] selftest R=%d D0=%d: true capture at row %d, "
              "bound %d, exact-pred %s, domination %s over %d rows"
              % ("PASS" if ok else "FAIL", R, D0, cap_row,
                 rep["capture_row_bound"], rep["capture_row_exact_pred"],
                 "held" if ok else "VIOLATED", k))
    return ok


def negctrl(R, D0):
    """Falsifiability controls: corrupt the state two ways; the certificate
    must FAIL both times.  (Corruptions in memory only; files untouched.)"""
    import copy
    base = os.path.join(RESULTS,
                        "amm12592_S_feedend_R%d_D0%d_boxeph.json" % (R, D0))
    st = json.load(open(base))
    m = st["m"]
    # control 1: blow the front cell up to super scale (freeze must fail)
    st1 = copy.deepcopy(st)
    v = int(st1["junk"][str(m)])
    st1["junk"][str(m)] = str(v * (1 << 60))
    # control 2: flip one interior cell positive (negativity must fail)
    st2 = copy.deepcopy(st)
    st2["junk"][str(1)] = str(-int(st2["junk"][str(1)]))
    tmpdir = os.environ.get("TMPDIR", "/tmp")
    ok = True
    for tag, s in (("front*2^60", st1), ("sign-flip", st2)):
        tp = os.path.join(tmpdir, "S_negctrl_%s.json" % tag.replace("*", "_"))
        json.dump(s, open(tp, "w"))
        real_load = globals()["load_state"]
        def fake_load(RR, DD, _tp=tp):
            stx = json.load(open(_tp))
            stx["junk"] = {int(t): int(v) for t, v in stx["junk"].items()}
            stx["j0"] = int(stx["j0"])
            return stx
        globals()["load_state"] = fake_load
        try:
            rep = certify(R, D0, verbose=False)
        finally:
            globals()["load_state"] = real_load
        caught = not rep["PASS"]
        ok = ok and caught
        print("[%s] negative control %s: certificate %s"
              % ("PASS" if caught else "FAIL", tag,
                 "correctly FAILS" if caught else "WRONGLY PASSES"))
    return ok


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "cert":
        certify(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == "selftest":
        selftest(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == "negctrl":
        negctrl(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == "all":
        ledger = []
        for R in (128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768):
            for D0 in (-(-R // 32), -(-R // 16)):
                p = os.path.join(RESULTS,
                                 "amm12592_S_feedend_R%d_D0%d_boxeph.json"
                                 % (R, D0))
                if os.path.exists(p):
                    ledger.append(certify(R, D0))
        outp = os.path.join(RESULTS,
                            "amm12592_S_superblock_certificate_boxeph.json")
        json.dump(ledger, open(outp, "w"), indent=1)
        npass = sum(1 for r in ledger if r["PASS"])
        print("\n%d/%d certificates PASS -> %s" % (npass, len(ledger), outp))
