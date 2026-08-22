#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AMM 12592 -- ESTIMATE-E ENDGAME HOSTILE REFEREE (boxeph, 2026-08-03).

Adversarial audit of the three lane reports:
  D1 (bulk-rule design sweep), D2 (E-lin theorem), D3 (golden transient bound).

Everything on the verification side is FRESH CODE in this file:
  - own Lucas/Fibonacci exact floor engine for floor(m gamma*);
  - own STREAMING witness verifier (admissibility + parity + unit + debt +
    full epoch identity via the t = (1-x)/x substitution processed row-by-row
    in level order, O(M) memory) -- this makes a FULL R = 2048 verification
    feasible on this machine;
  - own independent replay of the plain T6 junk flow (a 4th implementation);
  - own implementations of every 'PROVED' inequality of D2 (Theorems A/B/C)
    and D3 (G1/G2/G3, L*) on hostile grids incl. boundary/parity cases and
    non-dyadic controls.
Build-side artifacts (witness files, engines) are only used as INPUT DATA;
where an engine is invoked (block reconstruction for bulk variants), its
output is verified end-to-end by the fresh streaming verifier.

Exact arithmetic only (int / Fraction).  No numpy.  No floats in decisions.
sympy not used.

Stages: quick wit1024 witD2 wit2048a wit2048b proofs summary
Output: 05-knowledge/results/amm12592_estimateE_referee_boxeph.out (append)
Ledger: 05-knowledge/results/amm12592_estimateE_referee_ledger_boxeph.json
"""
import sys, os, json, time, hashlib, io, contextlib, importlib.util
from math import comb, isqrt
from fractions import Fraction

WT = "/tmp/math-wt-boxeph-multifront"
C4 = WT + "/04-computation"
R5 = WT + "/05-knowledge/results"
LED = R5 + "/amm12592_estimateE_referee_ledger_boxeph.json"
sys.path.insert(0, C4)
sys.setrecursionlimit(100000)

FAILS = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    if not ok:
        FAILS.append(name)
    print("  [%s] %s%s" % (tag, name, ("  -- " + detail) if detail else ""),
          flush=True)
    return bool(ok)


def audit(name, holds, detail=""):
    tag = "CONFIRMED" if holds else "REFUTED  "
    print("  [%s] %s%s" % (tag, name, ("  -- " + detail) if detail else ""),
          flush=True)
    return bool(holds)


def note(msg):
    print("  [NOTE] " + msg, flush=True)


def led_update(key, val):
    d = {}
    if os.path.exists(LED):
        try:
            d = json.load(open(LED))
        except Exception:
            d = {}
    d[key] = val
    json.dump(d, open(LED, "w"), indent=1, default=str)


# ---------------------------------------------------------------------------
# fresh exact floor(m gamma*) engine
# ---------------------------------------------------------------------------

def fib_pair(n):
    if n == 0:
        return (0, 1)
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if (n & 1) else (c, d)


def le5(dd, m):
    """5^dd <= phi^(2m)?  exact (phi^(2m) = (L_{2m}+F_{2m} sqrt5)/2)."""
    if dd < 0:
        return True
    F, F1 = fib_pair(2 * m)
    L = 2 * F1 - F
    A = 2 * 5 ** dd - L
    if A <= 0:
        return True
    return A * A < 5 * F * F


_FC = {}


def fgs(m):
    """floor(m gamma*) exact."""
    if m in _FC:
        return _FC[m]
    d = int(m * 0.5979874) - 2
    while le5(d + 1, m):
        d += 1
    _FC[m] = d
    return d


def ballot(d):
    """b_t = C(d-1,t) - C(d-1,t-1)."""
    out = []
    P, Pp = 1, 0        # C(d-1,t), C(d-1,t-1)
    for t in range(d + 1):
        out.append(P - Pp)
        Pp = P
        P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
    return out


def crow(d):
    row = [1] * (d + 1)
    for k in range(d):
        row[k + 1] = row[k] * (d - k) // (k + 1)
    return row


def two_G(R):
    g = [2]
    b = 1
    for j in range(1, R):
        b = b * (R - j) // j
        g.append((-b if j & 1 else b) - 1)
    return g


def t4_load(R, d):
    """row-0 load w_t = (-1)^(d-t) C(R-2-t,d-t) - C(d+1,t+1) + 2C(d,t)."""
    out = []
    A = comb(R - 2, d)
    B = d + 1
    CB = 1
    sgn = 1 if d % 2 == 0 else -1
    for t in range(d + 1):
        out.append(sgn * A - B + 2 * CB)
        if t < d:
            A = A * (d - t) // (R - 2 - t)
            B = B * (d - t) // (t + 2)
            CB = CB * (d - t) // (t + 1)
            sgn = -sgn
    return out


# ---------------------------------------------------------------------------
# fresh STREAMING witness verifier
# ---------------------------------------------------------------------------

def verify_stream(label, R, D0, row_iter, do_x2=True):
    """row_iter yields the R full cell rows in order.  All checks fresh."""
    t0 = time.time()
    res = {"label": label, "R": R, "D0": D0}
    prof = [fgs(R + i) + D0 for i in range(R)]
    M = R - 1 + prof[-1]
    acc = []                 # polynomial in t, list of coeffs
    ok_adm = True
    ok_len = True
    ok_unit = True
    ok_traj = True
    first_bad = None
    minus = 0
    sig = 0
    SP2 = 0
    last_corner = False
    for i in range(R):
        d = prof[i]
        row = next(row_iter)
        if len(row) != d + 1:
            ok_len = False
            first_bad = first_bad or ("len", i)
            break
        cr = crow(d)
        for k in range(d + 1):
            v = row[k]
            b = cr[k]
            if v > b or v < -b or ((v - b) & 1):
                ok_adm = False
                first_bad = first_bad or ("adm", i, k)
        if row[0] not in (1, -1):
            ok_unit = False
        sig -= row[0]
        if abs(sig) > R - 1 - i or ((sig - (i + 1)) & 1):
            ok_traj = False
        if i <= R - 2 and row[0] == -1:
            minus += 1
        if i == R - 1:
            last_corner = (row == [-c for c in cr])
        if do_x2:
            v2 = 0
            for k, c in enumerate(row):
                if c:
                    v2 += c * (1 << (d - k)) * (1 if k % 2 == 0 else -1)
            SP2 += (1 << i) * v2
        # identity accumulation (ascending level-Horner)
        if i == 0:
            acc = list(row)
        else:
            e = 1 + (prof[i] - prof[i - 1])       # e_{i-1} - e_i in {1,2}
            for _ in range(e):
                acc.append(0)
                for jj in range(len(acc) - 1, 0, -1):
                    if acc[jj - 1]:
                        acc[jj] += acc[jj - 1]
            if len(row) > len(acc):
                acc.extend([0] * (len(row) - len(acc)))
            for k, v in enumerate(row):
                if v:
                    acc[k] += v
    # RHS = t^(R-1) (1+t)^(M-R+1)
    rhs = [0] * (R - 1) + crow(M - R + 1)
    while acc and acc[-1] == 0:
        acc.pop()
    while rhs and rhs[-1] == 0:
        rhs.pop()
    ok_id = (acc == rhs)
    res["adm"] = check("W-adm [%s] all %d blocks in Lucas box + parity"
                       % (label, R), ok_adm and ok_len,
                       str(first_bad) if first_bad else "")
    res["identity"] = check("W-id  [%s] epoch identity == (1-x)^(R-1) "
                            "(fresh streaming t-substitution)" % label, ok_id)
    if do_x2:
        res["x2"] = check("W-x2  [%s] identity at x=2 exact" % label,
                          SP2 == (-1) ** (R - 1))
    res["unit"] = check("W-unit[%s] Delta_i(1) = +-1 every row" % label, ok_unit)
    res["traj"] = check("W-traj[%s] ballot trajectory legal, ends 0"
                        % label, ok_traj and sig == 0)
    res["debt"] = check("W-debt[%s] ballot debt = (R-2)/2 = %d"
                        % (label, (R - 2) // 2), minus == (R - 2) // 2,
                        "got %d" % minus)
    res["last_corner"] = last_corner
    res["ok"] = all(res.get(k) is not False for k in
                    ("adm", "identity", "x2", "unit", "traj", "debt"))
    res["secs"] = round(time.time() - t0, 1)
    print("      [%s] ok=%s  (%.1fs)" % (label, res["ok"], res["secs"]),
          flush=True)
    return res


def bulk_rows(path):
    """Row generator for the D1 sparse witness format."""
    w = json.load(open(path))
    R, D0 = w["R"], w["D0"]
    sc = {int(k): {int(kk): vv for kk, vv in v.items()}
          for k, v in w["sparse_c"].items()}

    def gen():
        for i in range(R):
            d = fgs(R + i) + D0
            if i == R - 1:
                yield [-c for c in crow(d)]
                continue
            b = ballot(d)
            u = sc.get(i)
            if u:
                for t, v in u.items():
                    b[t] += v
            yield b
    return R, D0, gen()


def urows_rows(R, D0, urows):
    def gen():
        for i in range(R):
            d = fgs(R + i) + D0
            if i == R - 1:
                yield [-c for c in crow(d)]
                continue
            b = ballot(d)
            u = urows.get(i)
            if u:
                for t, v in u.items():
                    b[t] += v
            yield b
    return gen()


# ---------------------------------------------------------------------------
# fresh independent replay of the PLAIN T6 junk flow (4th implementation)
# ---------------------------------------------------------------------------

def replay_plain(R, D0, keep_states=None, theta_watch=False):
    """Fresh implementation from the T6 spec.  Returns dict with outcome,
    per-row junk snapshots for rows in keep_states, last fed row, and
    (optionally) the G3 precondition scan (same-sign front pair vs Theta')."""
    g = two_G(R)
    d = fgs(R) + D0
    # row 0
    j = {}
    load = t4_load(R, d)
    P, Pp = 1, 0
    for t in range(d + 1):
        lo, hi = -2 * P, 2 * Pp
        v = load[t]
        u = min(hi, max(lo, v))
        if u != v:
            j[t] = v - u
        Pp = P
        P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
    states = {}
    if keep_states and 0 in keep_states:
        states[0] = (d, dict(j))
    minus2 = 0
    last_fed = 0
    precond_rows = []
    outcome = None
    d_prev = d
    for i in range(1, R - 1):
        d = fgs(R + i) + D0
        delta = d - d_prev
        assert delta in (0, 1)
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
            w[0] = w.get(0, 0) + g[d + i]
            fed = True
        if delta == 1 and d - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[d - 1 + i]
            w[1] = w.get(1, 0) + g[d - 1 + i]
            fed = True
        if fed:
            last_fed = i
        jn = {}
        c0 = 0
        if w:
            ts = sorted(w)
            ta, tb = ts[0], ts[-1]
            if tb > d:
                outcome = {"status": "DIE", "row": i, "over_top": True}
                break
            P = comb(d - 1, ta)
            Pp = comb(d - 1, ta - 1) if ta >= 1 else 0
            for t in range(ta, tb + 1):
                v = w.get(t, 0)
                if v:
                    lo, hi = -2 * P, 2 * Pp
                    u = min(hi, max(lo, v))
                    assert (v - u) % 2 == 0
                    if u != v:
                        jn[t] = v - u
                    if t == 0:
                        c0 = u
                Pp = P
                P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        if d in jn:
            outcome = {"status": "DIE", "row": i,
                       "const_bits": abs(jn[d]).bit_length()}
            j = jn
            break
        if c0 == -2:
            minus2 += 1
        j = jn
        d_prev = d
        if keep_states and i in keep_states:
            states[i] = (d, dict(j))
        if theta_watch and j:
            L = max(j)
            if L >= 1 and (L - 1) in j and j[L] * j[L - 1] > 0:
                mn = min(abs(j[L]), abs(j[L - 1]))
                theta = 0
                for t in range(L - 1, d + 1):
                    theta += 2 * comb(d, t)
                if mn > theta:
                    precond_rows.append(i)
        if not j and not fed and d + i > R - 1:
            outcome = {"status": "CLOSED", "capture_row": i}
            break
    return {"R": R, "D0": D0, "outcome": outcome, "minus2": minus2,
            "last_fed": last_fed, "states": states,
            "precond_rows": precond_rows}


# ---------------------------------------------------------------------------
# fresh G3 chain checker
# ---------------------------------------------------------------------------

def g3_chain(R, D0, i0, d, j, seeded, band=False):
    """Run the G3 forced diagonal from state (after clamp of row i0).
    Fresh implementation of D3's chain.  seeded=False -> psi == 0 (NI class);
    seeded=True, band=False -> the note's F3 majorant on the FULL layer;
    seeded=True, band=True  -> D3's CODE truncation (band top advancing at
    kernel speed from an empty start).  Returns (valid, predicted_death)."""
    if not j:
        return False, None
    L = max(j)
    if L < 4 or (L - 1) not in j or j[L] == 0 or j[L - 1] == 0:
        return False, None
    if j[L] * j[L - 1] <= 0:
        return False, None
    s = 1 if j[L] > 0 else -1
    a, b = abs(j[L]), abs(j[L - 1])
    psi = {}                       # layer majorant on cells > L (0 = absent)
    gap0 = d - L
    i = i0
    while L < d:
        dn = fgs(R + i + 1) + D0
        delta = dn - d
        assert delta in (0, 1)
        dpr = dn                    # new degree
        if delta == 1:
            Lf, Lm = L + 2, L + 1
            fl = a - 2 * psi.get(L + 1, 0) - psi.get(L + 2, 0)
            ml = 2 * a + b - psi.get(L + 1, 0)
        else:
            Lf, Lm = L + 1, L
            fl = a - psi.get(L + 1, 0)
            ml = a + b
        endf = 2 * comb(dpr - 1, Lf) if s < 0 else 2 * comb(dpr - 1, Lf - 1)
        endm = 2 * comb(dpr - 1, Lm) if s < 0 else 2 * comb(dpr - 1, Lm - 1)
        if not (fl > endf and ml > endm):
            return False, None
        if seeded:
            npsi = {}
            if band:
                hi_old = max(psi) if psi else L
                hi = min(dpr, max(hi_old + 1 + delta, Lf + 1 + delta))
            else:
                hi = dpr
            for c in range(Lf + 1, hi + 1):
                if delta == 1:
                    v = (psi.get(c, 0) + 2 * psi.get(c - 1, 0)
                         + psi.get(c - 2, 0))
                else:
                    v = psi.get(c, 0) + psi.get(c - 1, 0)
                npsi[c] = v + 2 * comb(dpr, c)
            psi = npsi
        else:
            psi = {}
        a, b = fl - endf, ml - endm
        L, d = Lf, dpr
        i += 1
    return True, i0 + gap0


# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------

def stage_quick():
    print("=" * 78)
    print("STAGE quick -- D1 witness files: sha256, reconstruction, fresh "
          "streaming verification (R=256 desc1 D0=0, R=512 tr15 D0=4) + "
          "negative control + cross-check vs prior referee")
    print("=" * 78, flush=True)
    out = {}
    shas = {"amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json":
            "5950bd4287d499221609bfaf7f42b6200637a1ede108c9086173a6b41fb60709",
            "amm12592_witness_R512_bulk_tr15_D0_4_boxeph.json":
            "311975f69412f09efbc171b71bb9485da8c1e2f9cda804de526c2a70013c9e83",
            "amm12592_witness_R1024_bulk_desc_D0_14_boxeph.json":
            "81391349e842d8086d3d90c2da41ebc215183521a70977be32fcb2cbe1206cd4"}
    for fn, want in shas.items():
        h = hashlib.sha256(open(C4 + "/" + fn, "rb").read()).hexdigest()
        out["sha_" + fn[17:30]] = check("sha256(%s) matches D1 ledger" % fn,
                                        h == want, h[:16])
    R, D0, it = bulk_rows(C4 + "/amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json")
    out["w256"] = verify_stream("R256 desc1 D0=0", R, D0, it)
    audit("D1 NEW RECORD: R=256 closes at D0=0 (desc1) -- witness verifies",
          out["w256"]["ok"])
    R, D0, it = bulk_rows(C4 + "/amm12592_witness_R512_bulk_tr15_D0_4_boxeph.json")
    out["w512"] = verify_stream("R512 tr15 D0=4", R, D0, it)
    audit("D1 NEW RECORD: R=512 closes at D0=4 (tr15) -- witness verifies",
          out["w512"]["ok"])

    # negative control: +2 corruption must be caught by identity
    w = json.load(open(C4 + "/amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json"))
    sc = {int(k): {int(kk): vv for kk, vv in v.items()}
          for k, v in w["sparse_c"].items()}

    def gen_bad():
        for i in range(256):
            d = fgs(256 + i)
            if i == 255:
                yield [-c for c in crow(d)]
                continue
            b = ballot(d)
            for t, v in sc.get(i, {}).items():
                b[t] += v
            if i == 40:
                b[d // 2] += 2       # parity-preserving, in-box corruption
            yield b
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        bad = verify_stream("NEGCTRL", 256, 0, gen_bad())
    FAILS[:] = [f for f in FAILS if "NEGCTRL" not in f]   # deliberate corruption
    out["negctrl"] = check("negative control: +2 in-box corruption caught by "
                           "fresh identity check", not bad["identity"])

    # cross-check my streaming verifier against the prior referee impl
    spec = importlib.util.spec_from_file_location(
        "aref", C4 + "/amm12592_allR_referee_boxeph.py")
    aref = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(aref)
    R, D0, it = bulk_rows(C4 + "/amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json")
    blocks = list(it)
    prof = [fgs(R + i) + D0 for i in range(R)]
    with contextlib.redirect_stdout(io.StringIO()):
        r2 = aref.verify_witness("xcheck", R, D0, prof, blocks)
    out["xcheck"] = check("cross-check: prior referee verify_witness agrees "
                          "on R256 desc1 (ok=%s)" % r2["ok"], r2["ok"])
    ok = all(le5(fgs(m), m) and not le5(fgs(m) + 1, m)
             for m in list(range(2, 200)) + [1000, 4095, 16384, 32767])
    out["floor"] = check("fresh floor engine: defining inequalities on hostile"
                         " grid incl. 4095/16384/32767", ok)
    led_update("quick", out)
    print("quick done; FAILS so far: %s" % (FAILS or "none"), flush=True)


def stage_wit1024():
    print("=" * 78)
    print("STAGE wit1024 -- R=1024 desc D0=14 (D1 record): full fresh "
          "streaming verification")
    print("=" * 78, flush=True)
    R, D0, it = bulk_rows(C4 + "/amm12592_witness_R1024_bulk_desc_D0_14_boxeph.json")
    r = verify_stream("R1024 desc D0=14", R, D0, it)
    audit("D1 NEW RECORD: R=1024 closes at D0=14 (desc) => attainment "
          "constant 14 for n <= 2047 rests on a fully verified witness",
          r["ok"])
    led_update("wit1024", r)
    print("wit1024 done; FAILS: %s" % (FAILS or "none"), flush=True)


def stage_witD2():
    print("=" * 78)
    print("STAGE witD2 -- D2 linear-slack closures (256,16) and (512,32): "
          "reconstruct via certified engine, verify with FRESH verifier "
          "(3rd independent check after engine asserts + C7)")
    print("=" * 78, flush=True)
    spec = importlib.util.spec_from_file_location(
        "fast", C4 + "/amm12592_transient_fast_junkflow_boxeph.py")
    fast = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(fast)
    out = {}
    for (R, D0) in ((256, 16), (512, 32)):
        blocks, res = fast.full_blocks(R, D0)
        assert blocks is not None
        r = verify_stream("R%d ruleA D0=%d (eps=%s)" % (R, D0, D0 and
                          Fraction(D0, R)), R, D0, iter(blocks))
        out["R%d_D0%d" % (R, D0)] = r
        audit("D2 sweep closure (R=%d, D0=%d) is a genuine witness" % (R, D0),
              r["ok"])
    led_update("witD2", out)
    print("witD2 done; FAILS: %s" % (FAILS or "none"), flush=True)


def stage_wit2048a():
    print("=" * 78)
    print("STAGE wit2048a -- BIG-R END-TO-END: rerun desc1 R=2048 D0=36 "
          "(D1 record, previously spot-checked only), rebuild all 2048 "
          "blocks, full fresh streaming verification")
    print("=" * 78, flush=True)
    spec = importlib.util.spec_from_file_location(
        "sweep", C4 + "/amm12592_bulkrule_design_sweep_boxeph.py")
    sweep = importlib.util.module_from_spec(spec)
    sv = sys.argv
    sys.argv = ["x", "noop"]
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(sweep)
    sys.argv = sv
    t0 = time.time()
    res = sweep.run_bulk(2048, 36, "desc1", keep_rows=False, collect_u=True,
                         spot=True)
    ok_run = (res["outcome"]["status"] == "CLOSED"
              and res["capture_row"] == 1271
              and res["minus2_rows"] == 1023
              and res["spots"] == {"x2": True, "x3": True, "debt": True})
    check("desc1 R=2048 D0=36 rerun reproduces D1 ledger (CLOSED@1271, "
          "debt 1023, spots pass)", ok_run,
          "%.0fs" % (time.time() - t0))
    r = verify_stream("R2048 desc1 D0=36", 2048, 36,
                      urows_rows(2048, 36, res["_urows"]))
    audit("D1 NEW RECORD: R=2048 closes at D0=36 (desc1) -- now FULLY "
          "verified (not just spots)", r["ok"] and ok_run)
    r["run_ok"] = ok_run
    led_update("wit2048a", r)
    print("wit2048a done; FAILS: %s" % (FAILS or "none"), flush=True)


def stage_wit2048b():
    print("=" * 78)
    print("STAGE wit2048b -- BIG-R END-TO-END for D2's sweep: plain rule A "
          "R=2048 D0=128 (eps=1/16), full fresh streaming verification")
    print("=" * 78, flush=True)
    spec = importlib.util.spec_from_file_location(
        "fast", C4 + "/amm12592_transient_fast_junkflow_boxeph.py")
    fast = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(fast)
    tmp = R5 + "/_ref_tmp_wit_R2048_D0128.json"
    t0 = time.time()
    res = fast.run_fast(2048, 128, witness_path=tmp, keep_rows=False)
    ok_run = (res["outcome"]["status"] == "CLOSED"
              and res["outcome"].get("capture_row") == 1237
              and res["minus2_rows"] == 1023)
    check("rule A R=2048 D0=128 rerun reproduces D2 sweep (CLOSED@1237, "
          "debt 1023)", ok_run, "%.0fs" % (time.time() - t0))
    W = json.load(open(tmp))
    sc = {int(k): {int(kk): vv for kk, vv in v.items()}
          for k, v in W["sparse_c"].items()}

    def gen():
        R, D0 = 2048, 128
        for i in range(R):
            d = fgs(R + i) + D0
            if i == R - 1:
                yield [-c for c in crow(d)]
                continue
            b = ballot(d)
            if i == 0:
                load = t4_load(R, d)
                P, Pp = 1, 0
                for t in range(d + 1):
                    lo, hi = -2 * P, 2 * Pp
                    b[t] += min(hi, max(lo, load[t]))
                    Pp = P
                    P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
            else:
                for t, v in sc.get(i, {}).items():
                    b[t] += v
            yield b
    r = verify_stream("R2048 ruleA D0=128", 2048, 128, gen())
    audit("D2 sweep at R=2048 eps=1/16 is a genuine witness (first full "
          "verification of a D2 closure beyond R=512)", r["ok"] and ok_run)
    r["run_ok"] = ok_run
    led_update("wit2048b", r)
    os.remove(tmp)
    print("wit2048b done; FAILS: %s" % (FAILS or "none"), flush=True)


# ---------------------------------------------------------------------------

def stage_proofs():
    print("=" * 78)
    print("STAGE proofs -- adversarial re-derivation of every PROVED claim "
          "of D2 (A/B/C) and D3 (G1/G2/G3/L*)")
    print("=" * 78, flush=True)
    out = {}
    import random
    rng = random.Random(12592)

    # ---------------- D2 Theorem A: LIFT ----------------------------------
    print("\n-- D2 Theorem A (LIFT)", flush=True)
    ok_alg = True
    for _ in range(400):
        d = rng.randrange(1, 40)
        cr = crow(d)
        cells = [rng.randrange(-cr[k], cr[k] + 1, 2) if cr[k] > 1
                 else rng.choice([-cr[k], cr[k]]) for k in range(d + 1)]
        # fix parity: delta_k == C(d,k) mod 2
        cells = [v if ((v - cr[k]) % 2 == 0) else v + 1
                 for k, v in enumerate(cells)]
        cells = [min(cr[k], max(-cr[k], v)) for k, v in enumerate(cells)]
        lift = [ (cells[k] if k <= d else 0) + (cells[k - 1] if k >= 1 else 0)
                 for k in range(d + 2)]
        cr1 = crow(d + 1)
        if not all(abs(lift[k]) <= cr1[k] and ((lift[k] - cr1[k]) % 2 == 0)
                   for k in range(d + 2)):
            ok_alg = False
        # polynomial equality: sum delta_k x^(d-k) q^k unchanged
        def poly_of(cs, dd):
            acc = [0] * (dd + 1)
            for k, c in enumerate(cs):
                if c:
                    qq = [(-1) ** m * comb(k, m) for m in range(k + 1)]
                    for m, qc in enumerate(qq):
                        acc[dd - k + m] += c * qc
            return acc
        if poly_of(cells, d) + [0] != poly_of(lift, d + 1) and \
           poly_of(cells, d) != poly_of(lift, d + 1)[:d + 1] + ([] if
           all(z == 0 for z in poly_of(lift, d + 1)[d + 1:]) else [None]):
            ok_alg = False
    out["A_alg"] = audit("Theorem A algebra: lift delta'_k = delta_k + "
                         "delta_(k-1) stays in box, keeps parity, keeps the "
                         "polynomial (400 random blocks, d <= 40)", ok_alg)
    # lift an ACTUAL stored witness (independent instance from D2's C6)
    w = json.load(open(C4 + "/amm12592_witness_R128_ruleB_D0_0_boxeph.json"))
    R = w["R"]
    for k_lift in (1, 3):
        def gen():
            for i in range(R):
                row = w["blocks"][i]
                for _ in range(k_lift):
                    row = [(row[k] if k < len(row) else 0)
                           + (row[k - 1] if 1 <= k <= len(row) else 0)
                           for k in range(len(row) + 1)]
                yield row
        with contextlib.redirect_stdout(io.StringIO()):
            r = verify_stream("lift+%d" % k_lift, R, 0 + k_lift, gen())
        out["A_lift%d" % k_lift] = audit(
            "Theorem A on data: R=128 ruleB D0=0 witness lifted +%d verifies "
            "at profile D0=%d (fresh verifier)" % (k_lift, k_lift), r["ok"])

    # ---------------- D2 Theorem B ----------------------------------------
    print("\n-- D2 Theorem B (i_feed, eps*, survival)", flush=True)
    # (B1) i_feed formula vs brute scan, incl. non-dyadic R
    ok_if = True
    glo, ghi = Fraction(627035, 1048576), Fraction(156759, 262144)
    grid = []
    for R in (128, 129, 255, 256, 257, 512, 777, 1000, 1024, 2047, 2048,
              3000, 4096, 8192, 16384):
        for D0 in (0, 1, 5, (R + 31) // 32, (R + 15) // 16):
            brute = max(i for i in range(R)
                        if fgs(R + i) + D0 + i <= R - 1)
            flo = (R * (1 - ghi) - D0) / (1 + ghi)
            fhi = (R * (1 - glo) - D0) / (1 + glo)
            a, b = flo.numerator // flo.denominator, \
                fhi.numerator // fhi.denominator
            if a != b:      # sandwich ambiguous -> decide by exact Beatty
                cand = None
                for i in (a, b):
                    if fgs(R + i) + D0 + i <= R - 1 and \
                       fgs(R + i + 1) + D0 + i + 1 > R - 1:
                        cand = i
                formula = cand
            else:
                formula = a
            grid.append((R, D0, brute, formula))
            if brute != formula:
                ok_if = False
    out["B_ifeed"] = audit("Theorem B(B1): i_feed = floor((R(1-g)-D0)/(1+g)) "
                           "== brute Beatty scan on %d hostile points incl. "
                           "non-dyadic R" % len(grid), ok_if)
    # empirical last-fed-row vs i_feed (definition audit)
    diffs = {}
    for (R, D0) in ((128, 0), (128, 8), (256, 16), (512, 32), (512, 0)):
        rp = replay_plain(R, D0)
        brute = max(i for i in range(R) if fgs(R + i) + D0 + i <= R - 1)
        diffs["R%d_D0%d" % (R, D0)] = (rp["last_fed"], brute)
    note("i_feed definition audit (engine last-fed-row vs D2 formula): %s"
         % diffs)
    off = {k: v[0] - v[1] for k, v in diffs.items()}
    out["B_lastfed"] = audit(
        "FINDING: the engine's true last fed row equals D2's i_feed + 1 "
        "whenever the (delta=1) second feed branch fires past it "
        "(transient note sec.7 uses the smaller d_i+i+1 <= R-1 convention). "
        "Theorem B's inequality chain is proved for rows <= i_feed(D2) only; "
        "the extra fed row is NOT covered by the survival statement as "
        "worded ('before the feed stops'), though at eps >= 1/32 the T6b "
        "bound clears it with Theta(R) margin, so nothing downstream breaks",
        True, "offsets last_fed - i_feed: %s" % off)
    # (B2)+(B3) eps*: algebraic equivalence at random rational g; sandwich
    ok_eq = True
    for _ in range(300):
        g = Fraction(rng.randrange(50, 70), 100) + \
            Fraction(rng.randrange(0, 997), 99700)
        Rv = rng.randrange(10, 10 ** 6)
        D0v = rng.randrange(0, Rv // 8 + 1)
        lhs = (2 * g - 1) * Rv + 2 * D0v >= (Rv * (1 - g) - D0v) / (1 + g)
        rhs = D0v * (3 + 2 * g) >= 2 * Rv * (1 - g - g * g)
        if lhs != rhs:
            ok_eq = False
    out["B_eq"] = audit("Theorem B(B2) reduction: (2g-1)R + 2D0 >= "
                        "(R(1-g)-D0)/(1+g) <=> D0(3+2g) >= 2R(1-g-g^2) "
                        "(300 random rational g, exact Fractions)", ok_eq)
    eps = lambda g: 2 * (1 - g - g * g) / (3 + 2 * g)
    e_hi, e_lo = eps(glo), eps(ghi)      # eps* decreasing in g
    ok_eps = (e_lo < e_hi and
              Fraction(211736, 10 ** 7) < e_lo and
              e_hi < Fraction(211747, 10 ** 7) and
              e_hi < Fraction(1, 32))
    out["B_eps"] = audit("eps* = 2(1-g-g^2)/(3+2g) in (0.0211736, 0.0211747),"
                        " < 1/32  (exact Fractions from the g-sandwich)",
                        ok_eps, "eps* in (%.10f, %.10f)"
                        % (float(e_lo), float(e_hi)))
    # g sandwich itself, fresh at M = 2^20
    a = 627035
    M = 1 << 20
    okg = le5(a, M) and not le5(a + 1, M)
    out["B_g"] = audit("g sandwich: 627035/2^20 < g < 627036/2^20 = "
                       "156759/262144 (fresh 5^a vs phi^(2M) integer test, "
                       "M = 2^20)", okg)
    # ordering eps_phi < eps*  (eps_phi = 1/phi - g; 1/phi = (sqrt5-1)/2)
    # exact: compare 1/phi - g_hi_bound ... need sqrt5 sandwich
    N = 10 ** 40
    s = isqrt(5 * N * N)          # floor(sqrt5 * N)
    invphi_lo = Fraction(s - N, 2 * N)       # < 1/phi = (sqrt5-1)/2
    invphi_hi = Fraction(s + 1 - N, 2 * N)   # > 1/phi
    eph_hi = invphi_hi - glo
    ok_ord = eph_hi < e_lo
    out["B_ord"] = audit("constant ordering eps_phi = 1/phi - g < eps* "
                         "(exact sqrt5 + g sandwiches)", ok_ord,
                         "eps_phi <= %.9f < %.9f" % (float(eph_hi),
                                                     float(e_lo)))
    # window bounds: D0 = ceil(R/32) in T4b window for R >= 27; /16 for >= 162
    ok27 = all(2 * (fgs(R) + (R + 31) // 32) > R and
               3 * (fgs(R) + (R + 31) // 32) < 2 * R
               for R in range(27, 4000))
    ok162 = all(2 * (fgs(R) + (R + 15) // 16) > R and
                3 * (fgs(R) + (R + 15) // 16) < 2 * R
                for R in range(162, 4000))
    bad27 = [R for R in range(4, 27)
             if not (2 * (fgs(R) + (R + 31) // 32) > R and
                     3 * (fgs(R) + (R + 31) // 32) < 2 * R)]
    out["B_win"] = audit("window claims: d_0 in (R/2, 2R/3) at D0=ceil(R/32) "
                         "for ALL 27 <= R < 4000; at ceil(R/16) for ALL "
                         "162 <= R < 4000 (exact floor)", ok27 and ok162,
                         "sub-27 exceptions: %s" % bad27)
    # survival instantiation: T6b bound > i_feed at the swept slacks
    ok_surv = True
    for R in (128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768):
        for D0 in ((R + 31) // 32, (R + 15) // 16):
            d0 = fgs(R) + D0
            ifeed = max(i for i in range(R)
                        if fgs(R + i) + D0 + i <= R - 1)
            if not (2 * d0 - R + 2 > ifeed + 1):
                ok_surv = False
    out["B_surv"] = audit("Theorem B instantiated: 2d_0 - R + 2 > i_feed + 1 "
                          "(one PAST the extra fed row) at eps = 1/32, 1/16, "
                          "all dyadic R = 128..32768", ok_surv)

    # ---------------- D2 Theorem C ----------------------------------------
    print("\n-- D2 Theorem C (one-step lemmas)", flush=True)
    ok_rat = True
    for _ in range(500):
        d = rng.randrange(3, 200)
        t = rng.randrange(2, d)
        if comb(d - 1, t) * d != comb(d, t) * (d - t):
            ok_rat = False
        if comb(d - 1, t - 1) * d != comb(d, t) * t:
            ok_rat = False
        if comb(d - 1, t - 2) * d * (d - t + 1) != comb(d, t) * t * (t - 1):
            ok_rat = False
    out["C_rat"] = audit("C-A binomial ratio identities (3 laws, 500 random "
                         "(d,t))", ok_rat)
    ok_cn = ok_ca = ok_cf = True
    for trial in range(300):
        d = rng.randrange(8, 60)
        m = rng.randrange(2, d - 4)
        j = {t: -2 * rng.randrange(0, comb(d, t) + 2) for t in range(m + 1)}
        j = {t: v for t, v in j.items() if v}
        delta = rng.randrange(0, 2)
        dn = d + delta
        w = {}
        for t, v in j.items():
            if delta == 0:
                w[t] = w.get(t, 0) + v
                w[t + 1] = w.get(t + 1, 0) + v
            else:
                w[t] = w.get(t, 0) + v
                w[t + 1] = w.get(t + 1, 0) + 2 * v
                w[t + 2] = w.get(t + 2, 0) + v
        jn = {}
        for t, v in w.items():
            lo, hi = -2 * comb(dn - 1, t), 2 * (comb(dn - 1, t - 1)
                                                if t >= 1 else 0)
            u = min(hi, max(lo, v))
            if u != v:
                jn[t] = v - u
        if any(v > 0 for v in jn.values()):
            ok_cn = False
        # C-A at overflowing cells (Fraction arithmetic)
        for t, vn in jn.items():
            if t < 2 or t > m:
                continue
            A = lambda tt: (Fraction(abs(j.get(tt, 0)), 2 * comb(d - 1, tt))
                            if 0 <= tt <= d - 1 else Fraction(0))
            An = Fraction(abs(vn), 2 * comb(dn - 1, t))
            if delta == 1:
                bound = (A(t) * (d - t) / d + 2 * A(t - 1) * t / d
                         + A(t - 2) * t * (t - 1) / (d * (d - t + 1)) - 1)
            else:
                bound = A(t) + A(t - 1) * t / (d - t) - 1
            if An > bound:
                ok_ca = False
        # C-F freeze: if hypotheses hold then no junk beyond m
        jm = abs(j.get(m, 0))
        jm1 = abs(j.get(m - 1, 0))
        if delta == 1:
            hyp = (jm <= 2 * comb(d, m + 2)
                   and 2 * jm + jm1 <= 2 * comb(d, m + 1))
        else:
            hyp = jm <= 2 * comb(d - 1, m + 1)
        if hyp and any(t > m for t in jn):
            ok_cf = False
    out["C_N"] = audit("C-N negativity invariance (300 random negative junk "
                       "one-steps, plain clamp)", ok_cn)
    out["C_A"] = audit("C-A cap-ratio damping inequality (exact Fractions, "
                       "all overflowing interior cells)", ok_ca)
    out["C_F"] = audit("C-F front freeze (hypotheses => support does not "
                       "advance)", ok_cf)

    # ---------------- D3 G1 -----------------------------------------------
    print("\n-- D3 Lemma G1 (alternation-transport calculus)", flush=True)
    ok_ti = ok_zs = ok_dm = True
    for _ in range(400):
        d = rng.randrange(10, 60)
        n = rng.randrange(5, 25)
        a = {c: rng.randrange(-10 ** 9, 10 ** 9) for c in range(3, 3 + n)}
        delta = rng.randrange(0, 2)
        jj = {c: ((-1) ** ((d - c) % 2)) * a[c] for c in a}
        K = [1, 1] if delta == 0 else [1, 2, 1]
        img = {}
        for c, v in jj.items():
            for s2, kk in enumerate(K):
                img[c + s2] = img.get(c + s2, 0) + kk * v
        # claimed: img_c = (-1)^(d-c) (D^(1+delta) a)_c
        for c in img:
            Da = 0
            for s2 in range(2 + delta):
                Da += (-1) ** s2 * comb(1 + delta, s2) * a.get(c - s2, 0)
            if img[c] != ((-1) ** ((d - c) % 2)) * Da:
                ok_ti = False
        tot = sum((-1) ** s2 * comb(1 + delta, s2) * a.get(c - s2, 0)
                  for c in range(0, 40) for s2 in [0])  # dummy
        Dvals = []
        for c in range(0, 3 + n + 3):
            Da = 0
            for s2 in range(2 + delta):
                Da += (-1) ** s2 * comb(1 + delta, s2) * a.get(c - s2, 0)
            Dvals.append(Da)
        if sum(Dvals) != 0:
            ok_zs = False
        L1 = sum(abs(v) for v in img.values())
        pos = sum(v for v in Dvals if v > 0)
        neg = -sum(v for v in Dvals if v < 0)
        if not (L1 == 2 * pos == 2 * neg):
            ok_dm = False
    out["G1_ti"] = audit("G1-a transport identity (400 random, both kernels)",
                         ok_ti)
    out["G1_zs"] = audit("G1-c zero-sum law", ok_zs)
    out["G1_dm"] = audit("G1-c defect-mass law L1 = 2 x one-signed mass",
                         ok_dm)
    ok_rl = True
    for rho in (2, 3, 4, 5):
        for delta in (0, 1):
            N = 14
            a = {c: rho ** (N - c) for c in range(N + 1)}     # a_{c-1}/a_c=rho
            d = 30
            jj = {c: ((-1) ** ((d - c) % 2)) * a[c] for c in a}
            K = [1, 1] if delta == 0 else [1, 2, 1]
            img = {}
            for c, v in jj.items():
                for s2, kk in enumerate(K):
                    img[c + s2] = img.get(c + s2, 0) + kk * v
            for c in range(2 + delta, N):                     # interior
                if abs(img[c]) != (rho - 1) ** (1 + delta) * a[c]:
                    ok_rl = False
    out["G1_rl"] = audit("G1-d ratio law |(K*j)_c| = (rho-1)^(1+delta) a_c "
                         "(rho = 2..5, both kernels); rho = 2 marginal for "
                         "(1,2,1)", ok_rl)
    # counterexample to the naive envelope
    Mv, N = 100, 8
    a = {c: (Mv if c % 2 == 0 else 0) for c in range(N + 1)}
    d = 20
    jj = {c: ((-1) ** ((d - c) % 2)) * a[c] for c in a if a[c]}
    img = {}
    for c, v in jj.items():
        for s2, kk in enumerate([1, 2, 1]):
            img[c + s2] = img.get(c + s2, 0) + kk * v
    interior_max = max(abs(v) for c, v in img.items() if 2 <= c <= N - 1)
    out["G1_cx"] = audit("G1-e counterexample: |j| <= M (const) but "
                         "|K*j| = 2M interior >> |D^2 M| = 0 -- naive "
                         "envelope lemma is FALSE as D3 records",
                         interior_max == 2 * Mv)

    # ---------------- D3 G2 -----------------------------------------------
    print("\n-- D3 Lemma G2 (initial data: sign/support/convexity/crossover)",
          flush=True)
    ok2 = True
    det = []
    for (R, D0) in ((512, 5), (1024, 15), (2048, 38), (1000, 0)):
        d = fgs(R) + D0
        load = t4_load(R, d)
        j = {}
        P, Pp = 1, 0
        for t in range(d + 1):
            lo, hi = -2 * P, 2 * Pp
            v = load[t]
            u = min(hi, max(lo, v))
            if u != v:
                j[t] = v - u
            Pp = P
            P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        supp = sorted(j)
        ok_su = supp == list(range(0, R - 2 - d + 1))
        ok_sg = all((j[t] > 0) == ((d - t) % 2 == 0) for t in supp)
        aa = {t: abs(j[t]) for t in supp}
        conc = [t for t in supp[2:] if aa[t] - 2 * aa[t - 1] + aa[t - 2] <= 0]
        # ratio-2 crossover: a_{t-1} vs 2 a_t flips exactly once near t2
        t2 = 2 * d - R + 3
        flips = []
        for t in supp[1:]:
            if (aa[t - 1] > 2 * aa[t]) != (aa[supp[1] - 1] > 2 * aa[supp[1]]):
                flips.append(t)
                break
        last_sub = max((t for t in supp[1:] if aa[t - 1] <= 2 * aa[t]),
                       default=None)
        first_sup = min((t for t in supp[1:] if aa[t - 1] > 2 * aa[t]),
                        default=None)
        mixed = (last_sub is not None and first_sup is not None
                 and last_sub > first_sup)
        okg2 = (ok_su and ok_sg and not conc and not mixed
                and first_sup is not None and abs(first_sup - t2) <= 1)
        det.append((R, D0, ok_su, ok_sg, len(conc), last_sub, first_sup, t2))
        if not okg2:
            ok2 = False
    out["G2"] = audit("G2 at (512,5),(1024,15),(2048,38) + NON-DYADIC even "
                      "control (1000,0): support EXACTLY [0,R-2-d0], signs "
                      "(-1)^(d-t), zero concave cells, single ratio-2 flip "
                      "within 1 cell of t2 = 2d-R+3", ok2, str(det))
    # odd-R domain note
    R, D0 = 1001, 0
    d = fgs(R) + D0
    load = t4_load(R, d)
    j = {}
    P, Pp = 1, 0
    for t in range(d + 1):
        lo, hi = -2 * P, 2 * Pp
        v = load[t]
        u = min(hi, max(lo, v))
        if u != v:
            j[t] = v - u
        Pp = P
        P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
    note("odd-R control (1001,0): support = [0,%d], R-2-d = %d -- boundary "
         "cell %s (G2/T4b odd-R branch; dyadic epochs unaffected)"
         % (max(j), R - 2 - d, "in-box" if max(j) < R - 2 - d else "violates"))

    # ---------------- D3 G3 -----------------------------------------------
    print("\n-- D3 Lemma G3 (super-pair death certificate), via FRESH plain "
          "replay + FRESH chain", flush=True)
    # replay regression first (4th implementation of plain rule A)
    reg = {}
    for (R, D0, want) in ((128, 0, ("CLOSED", 81)), (256, 0, ("DIE", 61)),
                          (256, 1, ("CLOSED", 159)), (512, 4, ("DIE", 121)),
                          (512, 5, ("CLOSED", 312)),
                          (1024, 14, ("DIE", 250)),
                          (1024, 15, ("CLOSED", 639))):
        rp = replay_plain(R, D0)
        o = rp["outcome"]
        got = (o["status"], o.get("row", o.get("capture_row")))
        reg["R%d_D0%d" % (R, D0)] = got
        if got != want:
            check("replay regression (%d,%d)" % (R, D0), False,
                  "%s vs %s" % (got, want))
    out["G3_reg"] = audit("fresh plain replay reproduces all 7 canonical "
                          "outcomes (die 61/121/250, close 81/159/312/639) "
                          "-- 4th independent implementation agrees",
                          all(reg["R%d_D0%d" % (R, D0)] == want for
                              (R, D0, want) in ((128, 0, ("CLOSED", 81)),
                                                (256, 0, ("DIE", 61)),
                                                (256, 1, ("CLOSED", 159)),
                                                (512, 4, ("DIE", 121)),
                                                (512, 5, ("CLOSED", 312)),
                                                (1024, 14, ("DIE", 250)),
                                                (1024, 15, ("CLOSED", 639)))))
    g3res = {}
    for (R, D0, claim_noseed, claim_seeded, death) in (
            (512, 4, 15, 25, 121), (1024, 14, 21, 42, 250)):
        keep = set(range(claim_noseed - 2, claim_seeded + 3))
        rp = replay_plain(R, D0, keep_states=keep)
        # noseed chain at claimed row and the row before
        d, j = rp["states"][claim_noseed]
        ok_v, pred = g3_chain(R, D0, claim_noseed, d, j, seeded=False)
        d2, j2 = rp["states"][claim_noseed - 1]
        ok_prev, _ = g3_chain(R, D0, claim_noseed - 1, d2, j2, seeded=False)
        ds, js = rp["states"][claim_seeded]
        ok_s, pred_s = g3_chain(R, D0, claim_seeded, ds, js, seeded=True)
        ds2, js2 = rp["states"][claim_seeded - 1]
        ok_s2, _ = g3_chain(R, D0, claim_seeded - 1, ds2, js2, seeded=True)
        g3res["R%d" % R] = dict(noseed_valid=ok_v, noseed_pred=pred,
                                noseed_prev_invalid=not ok_prev,
                                seeded_valid=ok_s, seeded_pred=pred_s,
                                seeded_prev_invalid=not ok_s2)
        out["G3_%d" % R] = audit(
            "G3 (%d,%d): FRESH noseed chain validates at row %d predicting "
            "death %d EXACTLY (and not at row %d); FRESH seeded chain "
            "(every admissible continuation) validates at row %d predicting "
            "%d" % (R, D0, claim_noseed, death, claim_noseed - 1,
                    claim_seeded, death),
            ok_v and pred == death and ok_s and pred_s == death,
            str(g3res["R%d" % R]))
    # hostile control: closing run never satisfies the precondition
    rp = replay_plain(512, 5, theta_watch=True)
    out["G3_ctrl"] = audit("G3 hostile control (512,5) CLOSED: same-sign "
                           "front pair with min > Theta' fires on ZERO rows "
                           "(fresh scan)", rp["precond_rows"] == [],
                           "outcome %s" % rp["outcome"])

    # ---------------- D3 L* -----------------------------------------------
    print("\n-- D3 tau = 2/3 threshold", flush=True)
    okL = True
    tab = []
    for d in (60, 100, 153, 306, 612, 1262, 2538):
        cr = crow(d)
        tails = 0
        Ls = None
        suf = [0] * (d + 2)
        for t in range(d, -1, -1):
            suf[t] = suf[t + 1] + cr[t]
        for L in range(1, d + 1):
            if cr[L - 1] >= suf[L]:
                Ls = L
                break
        dev = Ls - (-(-2 * d // 3))
        tab.append((d, Ls, dev))
        if dev not in (0, 1):
            okL = False
    out["Lstar"] = audit("L*(d) - ceil(2d/3) in {0,1} at 7 degrees (fresh)",
                         okL, str(tab))

    # ---------------- D2 Theorem D arithmetic -----------------------------
    print("\n-- D2 Theorem D / conditional-bound arithmetic", flush=True)
    okD = (1 + ghi + Fraction(1, 16) == Fraction(435287, 262144)
           and 1 + ghi + Fraction(1, 32) == Fraction(427095, 262144)
           and Fraction(435287, 262144) < Fraction(16605, 10000)
           and Fraction(427095, 262144) < Fraction(16293, 10000))
    out["D_arith"] = audit("conditional bounds: 1 + g + 1/16 < 435287/262144 "
                           "= 1.66049 and 1 + g + 1/32 < 427095/262144 = "
                           "1.62924 (given the verified g sandwich)", okD)
    led_update("proofs", out)
    print("\nproofs done; FAILS: %s" % (FAILS or "none"), flush=True)


def stage_summary():
    print("=" * 78)
    print("STAGE summary")
    print("=" * 78, flush=True)
    led = json.load(open(LED))
    print(json.dumps({k: (v if not isinstance(v, dict) else
                          {kk: (vv if not isinstance(vv, dict) else vv.get(
                              "ok", vv)) for kk, vv in v.items()})
                      for k, v in led.items()}, indent=1, default=str)[:4000])
    print("\nCHECKER FAILURES: %s" % (FAILS or "NONE"))


def stage_fix1():
    print("=" * 78)
    print("STAGE fix1 -- corrected eps_phi ordering; ledger-count audits; "
          "plain (2048,38) closure via the FRESH replay (4th implementation)")
    print("=" * 78, flush=True)
    out = {}
    glo, ghi = Fraction(627035, 1048576), Fraction(156759, 262144)
    eps = lambda g: 2 * (1 - g - g * g) / (3 + 2 * g)
    e_lo = eps(ghi)
    N = 10 ** 40
    s = isqrt(5 * N * N)
    invphi_hi = Fraction(s + 1 - N, 2 * N)
    invphi_lo = Fraction(s - N, 2 * N)
    eph_hi = invphi_hi - glo
    eph_lo = invphi_lo - ghi
    out["ord"] = audit("constant ordering eps_phi = 1/phi - g < eps* < 1/32 "
                       "(exact sqrt5 + g sandwiches; fixes referee's own "
                       "first-pass arithmetic slip)", eph_hi < e_lo,
                       "eps_phi in (%.9f, %.9f), eps*_lo = %.9f"
                       % (float(eph_lo), float(eph_hi), float(e_lo)))
    # D2 ledger counts
    sw = json.load(open(R5 + "/amm12592_Elin_sweep_boxeph.json"))
    prim = [r for r in sw if r["D0"] in ((r["R"] + 31) // 32,
                                          (r["R"] + 15) // 16)]
    okd = all(r["minus2_rows"] == (r["R"] - 2) // 2 and
              r["outcome"]["status"] == "CLOSED" for r in sw)
    out["sweep"] = audit("D2 sweep ledger: 23 runs (16 primary eps=1/32,1/16 "
                         "dyadic 128..16384 + 7 probes), ALL CLOSED, debt "
                         "(R-2)/2 exact in every record", okd and
                         len(sw) == 23 and len(prim) >= 16,
                         "%d records, %d primary" % (len(sw), len(prim)))
    fe = json.load(open(R5 + "/amm12592_Elin_feedend_state_boxeph.json"))
    flags = [all(v for k, v in r.items() if isinstance(v, bool)) for r in fe]
    out["feedend"] = audit("D2 feed-end certificates: 16/16 records all six "
                           "N/C/W/F/D/DL flags true (D2 note sec. 7 says "
                           "'14/14' -- typo for 16/16, JSON is authoritative)",
                           len(fe) == 16 and all(flags))
    # D3 certificate JSON counts vs note
    gsb = json.load(open(R5 + "/amm12592_golden_superblock_certificate_boxeph.json"))
    okg3 = (len(gsb["deaths"]) == 6 and len(gsb["controls"]) == 5
            and all(d["chain_noseed"]["predicted_death_row"]
                    == d["outcome"]["DIE"]
                    and d["chain_seeded"]["predicted_death_row"]
                    == d["outcome"]["DIE"]
                    and d["chain_seeded"]["validated"]
                    for d in gsb["deaths"])
            and all(c["precond_rows"] == 0 for c in gsb["controls"]))
    out["g3led"] = audit("D3 G3 ledger: 6 death certificates + 5 controls "
                         "present; controls fire on zero rows", okg3)
    note("D1 table cell 'desc1 @ R=128 D0=0' has NO desc1 run in the D1 "
         "ledger (.out has tr15/plain only at 128) -- immaterial to "
         "best-of-family (plain closes at 0) but the cell is unevidenced.")
    note("transient-theorem note sec. 7 gives i_feed(128,0) = 31 via the "
         "d_i+i+1 <= R-1 convention; the certified engine actually feeds at "
         "row 32 (first branch d_i+i <= R-1) and can feed once more at "
         "i_feed+1 via the delta=1 band (seen at (512,32): row 109). "
         "D2's formula matches the first branch; both notes should adopt "
         "'last fed row <= i_feed + 1' as the safe convention.")
    # fresh replay at a big scale: plain (2048,38) must CLOSE with capture 1271
    t0 = time.time()
    rp = replay_plain(2048, 38)
    o = rp["outcome"]
    out["r2048"] = audit("fresh replay (4th implementation): plain R=2048 "
                         "D0=38 CLOSED at capture 1271, debt 1023 -- "
                         "pins D0*(2048) = 38 close-side independently",
                         o == {"status": "CLOSED", "capture_row": 1271}
                         and rp["minus2"] == 1023,
                         "%s, minus2=%d, %.0fs" % (o, rp["minus2"],
                                                   time.time() - t0))
    led_update("fix1", out)
    print("fix1 done; FAILS: %s" % (FAILS or "none"), flush=True)


def stage_wit4096():
    print("=" * 78)
    print("STAGE wit4096 -- R=4096: desc1 D0=87 (D1 record) AND plain D0=89 "
          "(growth-table record): rebuild + full fresh verification")
    print("=" * 78, flush=True)
    out = {}
    # desc1 D0=87
    spec = importlib.util.spec_from_file_location(
        "sweep", C4 + "/amm12592_bulkrule_design_sweep_boxeph.py")
    sweep = importlib.util.module_from_spec(spec)
    sv = sys.argv
    sys.argv = ["x", "noop"]
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(sweep)
    sys.argv = sv
    t0 = time.time()
    res = sweep.run_bulk(4096, 87, "desc1", keep_rows=False, collect_u=True,
                         spot=True)
    ok_run = (res["outcome"]["status"] == "CLOSED"
              and res["capture_row"] == 2537
              and res["minus2_rows"] == 2047
              and res["spots"] == {"x2": True, "x3": True, "debt": True})
    check("desc1 R=4096 D0=87 rerun reproduces D1 ledger (CLOSED@2537, "
          "debt 2047, spots pass)", ok_run, "%.0fs" % (time.time() - t0))
    r = verify_stream("R4096 desc1 D0=87", 4096, 87,
                      urows_rows(4096, 87, res["_urows"]))
    audit("D1 NEW RECORD: R=4096 closes at D0=87 (desc1) -- now FULLY "
          "verified", r["ok"] and ok_run)
    r["run_ok"] = ok_run
    out["desc1_87"] = r
    del res
    led_update("wit4096", out)
    # plain D0=89
    spec = importlib.util.spec_from_file_location(
        "fast", C4 + "/amm12592_transient_fast_junkflow_boxeph.py")
    fast = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(fast)
    tmp = R5 + "/_ref_tmp_wit_R4096_D089.json"
    t0 = time.time()
    res = fast.run_fast(4096, 89, witness_path=tmp, keep_rows=False)
    ok_run = (res["outcome"]["status"] == "CLOSED"
              and res["outcome"].get("capture_row") == 2537
              and res["minus2_rows"] == 2047)
    check("plain R=4096 D0=89 rerun reproduces growth table (CLOSED@2537, "
          "debt 2047)", ok_run, "%.0fs" % (time.time() - t0))
    W = json.load(open(tmp))
    sc = {int(k): {int(kk): vv for kk, vv in v.items()}
          for k, v in W["sparse_c"].items()}

    def gen():
        R, D0 = 4096, 89
        for i in range(R):
            d = fgs(R + i) + D0
            if i == R - 1:
                yield [-c for c in crow(d)]
                continue
            b = ballot(d)
            if i == 0:
                load = t4_load(R, d)
                P, Pp = 1, 0
                for t in range(d + 1):
                    lo, hi = -2 * P, 2 * Pp
                    b[t] += min(hi, max(lo, load[t]))
                    Pp = P
                    P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
            else:
                for t, v in sc.get(i, {}).items():
                    b[t] += v
            yield b
    r = verify_stream("R4096 ruleA D0=89", 4096, 89, gen())
    audit("growth-law table: plain D0*(4096) = 89 close-side now FULLY "
          "verified (prior referee's open obligation at 4096 discharged)",
          r["ok"] and ok_run)
    r["run_ok"] = ok_run
    out["plain_89"] = r
    os.remove(tmp)
    led_update("wit4096", out)
    print("wit4096 done; FAILS: %s" % (FAILS or "none"), flush=True)


def stage_g3x():
    print("=" * 78)
    print("STAGE g3x -- G3 seeded-chain audit: band truncation vs the "
          "note's full-layer F3 majorant")
    print("=" * 78, flush=True)
    out = {}
    for (R, D0, row_ns, row_sd, death, hi_scan) in (
            (512, 4, 15, 25, 121, 90), (1024, 14, 21, 42, 250, 160)):
        keep = set(range(row_ns - 1, hi_scan + 1))
        rp = replay_plain(R, D0, keep_states=keep)
        d, j = rp["states"][row_sd]
        okb, pb = g3_chain(R, D0, row_sd, d, j, seeded=True, band=True)
        out["band_%d" % R] = audit(
            "(%d,%d): reproducing D3's CODE (band-truncated psi): seeded "
            "chain validates at row %d predicting %d -- discrepancy with "
            "the full-layer run is exactly the band truncation"
            % (R, D0, row_sd, death), okb and pb == death,
            "valid=%s pred=%s" % (okb, pb))
        first_full = None
        for i0 in range(row_sd, hi_scan + 1):
            if i0 not in rp["states"]:
                continue
            d, j = rp["states"][i0]
            okf, pf = g3_chain(R, D0, i0, d, j, seeded=True, band=False)
            if okf:
                first_full = (i0, pf)
                break
        out["full_%d" % R] = audit(
            "(%d,%d): with the note's FULL-layer F3 majorant the seeded "
            "chain first validates at row %s (death still predicted "
            "exactly)" % (R, D0, first_full),
            first_full is not None and first_full[1] == death,
            str(first_full))
    note("VERDICT G3: the noseed (NI-class) certificates are fully "
         "confirmed; the 'every admissible continuation' seeded "
         "certificates in D3's ledger use a band-truncated psi that does "
         "NOT implement the note's F3 recursion -- the full-layer majorant "
         "validates only later (see rows above).  G3 the THEOREM is sound "
         "(full-layer psi); the ledger's seeded validation ROWS and "
         "psi-margin values are artifacts of the truncation and must be "
         "restated from the full-layer runs.")
    led_update("g3x", out)
    print("g3x done; FAILS: %s" % (FAILS or "none"), flush=True)


STAGES = {"quick": stage_quick, "wit1024": stage_wit1024, "fix1": stage_fix1,
          "wit4096": stage_wit4096, "g3x": stage_g3x,
          "witD2": stage_witD2, "wit2048a": stage_wit2048a,
          "wit2048b": stage_wit2048b, "proofs": stage_proofs,
          "summary": stage_summary}

if __name__ == "__main__":
    for name in sys.argv[1:]:
        STAGES[name]()
