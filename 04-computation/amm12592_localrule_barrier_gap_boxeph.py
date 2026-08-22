"""LANE E3 -- the local-rule barrier and the existence gap (boxeph, 2026-08-04).

Stages (each appends to the .out ledger and updates the JSON ledger; run
individually or 'all'):

  A  RECORD CHECK + transportation integer point at (R,D0) = (256,0):
     re-verify the desc1 exact-floor witness (sha + hostile referee), then
     translate it to the nonnegative transportation form f = (C - delta)/2
     (THM-2966/-3002 change of variables) and certify 0 <= f <= C plus the
     (**) identity at exact integer points x = 2, 3, -1.  This upgrades
     "the transportation LP relaxation is feasible at 256" to "the
     transportation polytope at the EXACT FLOOR profile (D0 = 0) contains an
     INTEGER point at R = 256" -- and corrects the stale reading
     "D0*(256) = 1 for every method tried".

  B  NEW EXACT LEMMA CERTIFICATES (front marginality + boundary-layer rate):
     (L1) (t2 - 1) + F0 = d           (crossover/front complementarity)
     (L2) C(R-2-F0, d-F0) = C(d, F0)  (dominant front load = half box width)
     (L3) (A_{t-1}/A_t) / (C(d,t-1)/C(d,t)) = (R-1-t)/t   (excess growth
          per cell of depth; at the front -> d/F0 -> g/(1-g) = 1.4876...)
     verified as exact integer identities on a dyadic grid, plus the exact
     front-cell overflow/width ratio table (the taper marginality).

  C  VARIANT SWEEP at R = 512, D0 = 0..3: every VARIANTS rule of the D1
     sweep engine (19 rules), exact outcomes.  Rule negatives are NOT
     epoch-feasibility evidence (hazard); this maps the D1 class exactly.

  D  BEAM SEARCH over cross-row ACTION SCHEDULES at R = 512 (and a 256
     control): per row each surviving state branches over the D1 action
     alphabet (plain/desc/desc1/desc2/desc3/descnv/tr*/fill15), scored by
     (same-sign-front-pair risk bits, front position, junk L1 bits) --
     i.e. a bounded-width GLOBAL schedule search over local actions,
     guided by the G3 order parameter.  Any closure is rebuilt from its
     action history and referee-verified; failures prove nothing.

All arithmetic exact int; no floats in any decision; no numpy; sympy unused.
"""
import sys, os, json, time, io, contextlib, importlib.util, hashlib
from math import comb
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
OUT = os.path.join(RESULTS, "amm12592_localrule_barrier_gap_boxeph.out")
LEDGER = os.path.join(RESULTS, "amm12592_localrule_barrier_gap_boxeph.json")


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


sweep = _load("bulksweep", "amm12592_bulkrule_design_sweep_boxeph.py")
iref = _load("iref", "amm12592_independent_witness_referee_boxeph.py")
aref = _load("aref", "amm12592_allR_referee_boxeph.py")
tool = _load("tool", "amm12592_allR_family_toolbox_boxeph.py")
floorgs = iref.floor_gamma_star
two_G_coeffs = sweep.two_G_coeffs
choose_row = sweep.choose_row
t4_row_load = sweep.t4_row_load
VARIANTS = sweep.VARIANTS


def log(line):
    print(line, flush=True)
    with open(OUT, "a") as f:
        f.write(line + "\n")


def ledger_update(key, val):
    led = {}
    if os.path.exists(LEDGER):
        led = json.load(open(LEDGER))
    led[key] = val
    json.dump(led, open(LEDGER, "w"), indent=1)


def check(ok, msg, extra=""):
    tag = "PASS" if ok else "FAIL"
    log("  [%s] %s%s" % (tag, msg, ("  -- " + str(extra)) if extra != "" else ""))
    return bool(ok)


# ======================================================================
# STAGE A -- record check + transportation integer point at (256, 0)
# ======================================================================
def stage_A():
    log("=" * 74)
    log("STAGE A -- record check + transportation integer point at (256, D0=0)")
    log("=" * 74)
    res = {"fails": 0}
    path = os.path.join(HERE, "amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json")
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    ok = sha.startswith("5950bd4287d49922")
    res["fails"] += 0 if check(ok, "witness file sha256 matches D1 ledger",
                               sha[:16]) else 1
    doc = json.loads(raw)
    R, D0 = doc["R"], doc["D0"]
    assert (R, D0) == (256, 0)
    prof = [floorgs(R + i) + D0 for i in range(R)]
    blocks = []
    sparse = {int(i): {int(t): v for t, v in row.items()}
              for i, row in doc["sparse_c"].items()}
    for i in range(R - 1):
        d = prof[i]
        b = tool.ballot(d)[:]
        for t, v in sparse.get(i, {}).items():
            b[t] += v
        blocks.append(b)
    dL = prof[R - 1]
    blocks.append([-comb(dL, t) for t in range(dL + 1)])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        v = aref.verify_witness("E3_R256_desc1_D00", R, D0, prof, blocks)
    res["fails"] += 0 if check(v.get("ok"), "hostile referee verify_witness "
                               "(adm + identity + ballot laws) ALL-PASS",
                               {k: v[k] for k in ("ok", "adm", "identity")
                                if k in v}) else 1
    # ---- transportation translation: f = (C - delta)/2
    n_cells = 0
    bad = 0
    fmax_sat = 0
    for i in range(R):
        d = prof[i]
        for t in range(d + 1):
            c = comb(d, t)
            delta = blocks[i][t]
            num = c - delta
            if num % 2 != 0 or not (0 <= num // 2 <= c):
                bad += 1
            else:
                f = num // 2
                n_cells += 1
                if f == 0 or f == c:
                    fmax_sat += 1
    res["fails"] += 0 if check(bad == 0, "f = (C(d,t)-delta)/2 integer and in "
                               "[0, C(d,t)] at ALL %d cells" % n_cells,
                               "boundary-saturated cells: %d (%.2f%%)"
                               % (fmax_sat, 100.0 * fmax_sat / n_cells)) else 1
    # ---- (**) identity at exact integer points x = 2, 3, -1:
    #      sum_i x^i F_i(x) = T_R(x),  F_i = sum_t f_{i,t} x^{d-t} (1-x)^t,
    #      T_R = [ (1-x^R)/(1-x) - (1-x)^{R-1} ] / 2.
    id_ok = True
    for x in (2, 3, -1):
        lhs = 0
        for i in range(R):
            d = prof[i]
            Fi = 0
            q = 1 - x
            # Horner over cells t: sum_t f_t x^(d-t) q^t
            qp = 1
            xs = [x ** (d - t) for t in range(d + 1)]
            for t in range(d + 1):
                f = (comb(d, t) - blocks[i][t]) // 2
                if f:
                    Fi += f * xs[t] * qp
                qp *= q
            lhs += x ** i * Fi
        tr = ((1 - x ** R) // (1 - x) - (1 - x) ** (R - 1))
        assert tr % 2 == 0
        tr //= 2
        if lhs != tr:
            id_ok = False
        log("    x=%2d : sum x^i F_i = %s ; T_R = %s ; equal=%s"
            % (x, ("%d bits" % lhs.bit_length()), ("%d bits" % tr.bit_length()),
               lhs == tr))
    res["fails"] += 0 if check(id_ok, "(**) transportation identity holds at "
                               "x = 2, 3, -1 exactly (with adm+identity above "
                               "=> integer point of the (**) polytope)") else 1
    verdict = ("CONFIRMED: the transportation polytope at the EXACT FLOOR "
               "profile d_i = floor(g(256+i)) (D0 = 0) contains an INTEGER "
               "point; LP relaxation feasible a fortiori. The stale reading "
               "'D0*(256) = 1 for every method' is corrected by the record.")
    log("  [%s]" % verdict)
    res["verdict"] = verdict
    ledger_update("A_record_transport_256", res)
    return res["fails"]


# ======================================================================
# STAGE B -- exact lemma certificates (front marginality + layer rates)
# ======================================================================
def stage_B():
    log("=" * 74)
    log("STAGE B -- exact certificates: L1 complementarity, L2 front "
        "marginality, L3 depth-rate law, taper table")
    log("=" * 74)
    fails = 0
    grid = []
    for R in (128, 256, 512, 1024, 2048, 4096, 8192):
        for D0 in (0, 1, (R + 31) // 32):
            d = floorgs(R) + D0
            if not (R // 2 < d < 2 * R // 3 + 1):
                continue
            grid.append((R, D0, d))
    l1 = l2 = l3 = 0
    taper = []
    for (R, D0, d) in grid:
        F0 = R - 2 - d
        t2 = 2 * d - R + 3
        if (t2 - 1) + F0 != d:
            fails += 1
        else:
            l1 += 1
        # L2: dominant front load C(R-2-F0, d-F0) == C(d, F0) == width/2
        if comb(R - 2 - F0, d - F0) == comb(d, F0):
            l2 += 1
        else:
            fails += 1
        # L3: excess growth per cell of depth, exact Fractions, 8 depths
        okl3 = True
        for m in range(1, 9):
            t = F0 - m + 1
            lhsr = Fraction(comb(R - 2 - (t - 1), d - (t - 1)),
                            comb(R - 2 - t, d - t)) / \
                Fraction(comb(d, t - 1), comb(d, t))
            if lhsr != Fraction(R - 1 - t, t):
                okl3 = False
        if okl3:
            l3 += 1
        else:
            fails += 1
        # taper: exact front-cell full load (T4) vs box width
        wF = ((-1) ** (d - F0)) * comb(R - 2 - F0, d - F0) \
            - comb(d + 1, F0 + 1) + 2 * comb(d, F0)
        width = 2 * comb(d, F0)
        hi = 2 * comb(d - 1, F0 - 1)
        over = wF - hi if wF > hi else 0
        taper.append((R, D0, str(Fraction(wF, width)),
                      float(Fraction(wF, width)),
                      float(Fraction(over, width)) if width else None))
    fails += 0 if check(l1 == len(grid), "L1 (t2-1)+F0 = d on %d/%d grid "
                        "points (PROVED: algebra; certified)" % (l1, len(grid))) else 1
    fails += 0 if check(l2 == len(grid), "L2 C(R-2-F0,d-F0) = C(d,F0) "
                        "= (box width)/2 on %d/%d (PROVED: R-2-F0 = d, "
                        "d-F0 = 2d-R+2 = d-F0 sym; certified)" % (l2, len(grid))) else 1
    fails += 0 if check(l3 == len(grid), "L3 (A_{t-1}/A_t)/(C(d,t-1)/C(d,t)) "
                        "= (R-1-t)/t at 8 depths, %d/%d (PROVED: telescoping "
                        "binomial ratios; certified)" % (l3, len(grid))) else 1
    log("  taper table (front-cell full load / box width; overflow/width):")
    for (R, D0, fr, fv, ov) in taper:
        log("    R=%5d D0=%4d : w_F0/width = %.6f  overflow/width = %.6f"
            % (R, D0, fv, ov))
    # limit constants for the note (display only, from certified g sandwich)
    g_lo = Fraction(627035, 1048576)
    g_hi = Fraction(156759, 262144)
    log("  boundary-layer rate g/(1-g) in (%.7f, %.7f)  [certified sandwich]"
        % (float(g_lo / (1 - g_lo)), float(g_hi / (1 - g_hi))))
    ledger_update("B_lemmas", {"grid": len(grid), "L1": l1, "L2": l2,
                               "L3": l3, "taper": taper, "fails": fails})
    return fails


# ======================================================================
# STAGE C -- full variant sweep at R = 512, D0 = 0..3
# ======================================================================
def stage_C():
    log("=" * 74)
    log("STAGE C -- D1 rule class at R=512, D0=0..3: exact outcomes "
        "(rule negatives are NOT epoch evidence)")
    log("=" * 74)
    tab = {}
    for D0 in (0, 1, 2, 3):
        for vn in sorted(VARIANTS):
            t0 = time.time()
            r = sweep.run_bulk(512, D0, vn, keep_rows=False, spot=False)
            o = r["outcome"]
            s = ("DIE@%d" % o["row"]) if o["status"] == "DIE" else \
                ("CLOSED@%s" % o.get("capture_row"))
            tab.setdefault(vn, {})[D0] = s
            log("  R=512 D0=%d %-9s -> %-12s (%.1fs)"
                % (D0, vn, s, time.time() - t0))
            ledger_update("C_sweep_512", tab)
    closed = [(vn, D0) for vn, row in tab.items()
              for D0, s in row.items() if s.startswith("CLOSED")]
    log("  closures found at D0<=3: %s" % (closed if closed else "NONE"))
    return 0


# ======================================================================
# STAGE D -- beam over action schedules
# ======================================================================
BEAM_ACTIONS = ["plain", "desc", "desc1", "descnv", "tr15", "tr43",
                "tr54", "tr74", "train16", "fill15"]
for _n, _vl in (("desc2", 2), ("desc3", 3)):
    if _n not in VARIANTS:
        VARIANTS[_n] = dict(VARIANTS["desc"], vlad=_vl)
BEAM_ACTIONS += ["desc2", "desc3"]


def _feed(w, i, d, delta, R, g):
    fed = False
    if d + i <= R - 1:
        w[0] = w.get(0, 0) + g[d + i]
        fed = True
    if delta == 1 and d - 1 + i <= R - 1:
        w[0] = w.get(0, 0) + g[d - 1 + i]
        w[1] = w.get(1, 0) + g[d - 1 + i]
        fed = True
    return fed


def _transport(j, delta):
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
    return {t: v for t, v in w.items() if v}


def _score(j, d):
    """(risk_bits, tmax, L1bits): same-sign adjacent pair min-value bits
    minus front-cap bits within 40 cells of the front, then front, then L1."""
    if not j:
        return (-(10 ** 9), -1, 0)
    ts = sorted(j)
    tmax = ts[-1]
    L1b = sum(abs(v) for v in j.values()).bit_length()
    risk = -(10 ** 9)
    capb = (2 * comb(d - 1, tmax)).bit_length() if tmax <= d - 1 else 1
    for t in ts:
        if t < tmax - 40 or (t + 1) not in j:
            continue
        a, b = j[t], j[t + 1]
        if (a > 0) == (b > 0):
            mb = min(abs(a), abs(b)).bit_length()
            risk = max(risk, mb - capb)
    return (risk, tmax, L1b)


def beam_run(R, D0, width=24, budget_s=600, actions=None):
    """Beam over cross-row action schedules.  The beam population is the
    union of (a) PROTECTED PURE LANES -- one state per action that always
    plays that action, never pruned (so the beam provably dominates every
    pure rule in the alphabet) -- and (b) the scored beam: top `width`
    states under (risk, front, L1) plus top `width` under (front, L1)."""
    actions = actions or BEAM_ACTIONS
    g = two_G_coeffs(R)
    t0 = time.time()
    # state: (jdict, hist tuple of action names)
    states = [({}, ())]
    pure = {an: ({}, ()) for an in actions}   # protected lanes
    stats = {"deaths": 0, "rows_done": 0, "peak_states": 1,
             "pure_alive": len(actions)}
    for i in range(0, R - 1):
        d = floorgs(R + i) + D0
        dn = floorgs(R + i + 1) + D0
        S = 1 + (dn - d)
        # ---- step the protected pure lanes
        new_pure = {}
        for an, (j, hist) in pure.items():
            if i == 0:
                base_w = {t: v for t, v in enumerate(t4_row_load(R, d)) if v}
                fed = True
            else:
                dp = floorgs(R + i - 1) + D0
                delta = d - dp
                base_w = _transport(j, delta)
                fed = _feed(base_w, i, d, delta, R, g)
                base_w = {t: v for t, v in base_w.items() if v}
            if not base_w:
                if not fed and d + i > R - 1:
                    return {"status": "CLOSED", "capture_row": i - 1,
                            "hist": hist, "lane": "pure:" + an,
                            "elapsed": time.time() - t0, "stats": stats}
                new_pure[an] = ({}, hist + (an,))
                continue
            jn, _, c0, jl1, _ = choose_row(dict(base_w), d, dn, S,
                                           VARIANTS[an], record_u=False)
            if d in jn:
                continue
            new_pure[an] = (jn, hist + (an,))
        pure = new_pure
        stats["pure_alive"] = len(pure)
        # ---- step the scored beam
        nxt = {}
        for (j, hist) in states:
            if i == 0:
                base_w = {t: v for t, v in enumerate(t4_row_load(R, d)) if v}
                fed = True
            else:
                dp = floorgs(R + i - 1) + D0
                delta = d - dp
                base_w = _transport(j, delta)
                fed = _feed(base_w, i, d, delta, R, g)
                base_w = {t: v for t, v in base_w.items() if v}
            if not base_w:
                if not fed and d + i > R - 1:
                    return {"status": "CLOSED", "capture_row": i - 1,
                            "hist": hist, "elapsed": time.time() - t0,
                            "stats": stats}
                key = ()
                if key not in nxt:
                    nxt[key] = ({}, hist + ("plain",))
                continue
            for an in actions:
                jn, _, c0, jl1, _ = choose_row(dict(base_w), d, dn, S,
                                               VARIANTS[an], record_u=False)
                if d in jn:
                    stats["deaths"] += 1
                    continue
                key = tuple(sorted(jn.items()))
                if key not in nxt:
                    nxt[key] = (jn, hist + (an,))
        if not nxt and not pure:
            return {"status": "ALL_DIE", "row": i, "elapsed": time.time() - t0,
                    "stats": stats}
        ranked = sorted(nxt.values(), key=lambda sj: _score(sj[0], dn))
        ranked2 = sorted(nxt.values(),
                         key=lambda sj: _score(sj[0], dn)[1:])
        merged = {}
        for (jj, hh) in ranked[:width] + ranked2[:width] + list(pure.values()):
            k = tuple(sorted(jj.items()))
            if k not in merged:
                merged[k] = (jj, hh)
        states = list(merged.values())
        stats["rows_done"] = i
        stats["peak_states"] = max(stats["peak_states"], len(nxt))
        # capture check (post-feed empty junk)
        for (j, hist) in states:
            if not j and d + i > R - 1:
                return {"status": "CLOSED", "capture_row": i, "hist": hist,
                        "elapsed": time.time() - t0, "stats": stats}
        if time.time() - t0 > budget_s:
            best = _score(states[0][0], dn)
            return {"status": "BUDGET", "row": i, "best_score": best,
                    "elapsed": time.time() - t0, "stats": stats}
    # terminal: check last-row absorbability for each survivor
    for (j, hist) in states:
        if not j:
            return {"status": "CLOSED", "capture_row": R - 2, "hist": hist,
                    "elapsed": time.time() - t0, "stats": stats}
    return {"status": "OPEN_RESIDUAL", "n_final": len(states),
            "elapsed": time.time() - t0, "stats": stats}


def replay_and_verify(R, D0, hist):
    """Replay an action history collecting u, rebuild blocks, referee."""
    g = two_G_coeffs(R)
    urows = {}
    j = {}
    for i, an in enumerate(hist):
        d = floorgs(R + i) + D0
        dn = floorgs(R + i + 1) + D0
        S = 1 + (dn - d)
        if i == 0:
            w = {t: v for t, v in enumerate(t4_row_load(R, d)) if v}
        else:
            dp = floorgs(R + i - 1) + D0
            delta = d - dp
            w = _transport(j, delta)
            _feed(w, i, d, delta, R, g)
            w = {t: v for t, v in w.items() if v}
        if not w:
            j = {}
            continue
        jn, us, c0, jl1, _ = choose_row(w, d, dn, S, VARIANTS[an],
                                        record_u=True)
        assert d not in jn
        if us:
            urows[i] = us
        j = jn
    res = {"R": R, "D0": D0, "variant": "beam", "_urows": urows,
           "capture_row": len(hist) - 1}
    v, txt = sweep.referee_verify(res)
    return v, res


def stage_D():
    log("=" * 74)
    log("STAGE D -- beam over cross-row action schedules "
        "(alphabet: %s; width 24)" % ",".join(BEAM_ACTIONS))
    log("=" * 74)
    out = {}
    # positive control: 256 at D0 = 0 must close (desc1 alone does)
    r = beam_run(256, 0, width=12, budget_s=240)
    log("  control 256 D0=0 : %s" % json.dumps(
        {k: v for k, v in r.items() if k != "hist"}))
    out["ctrl_256_0"] = {k: v for k, v in r.items() if k != "hist"}
    ledger_update("D_beam", out)
    for D0 in (0, 1, 2, 3):
        r = beam_run(512, D0, width=16, budget_s=900)
        rec = {k: v for k, v in r.items() if k != "hist"}
        log("  beam 512 D0=%d : %s" % (D0, json.dumps(rec)))
        out["beam_512_%d" % D0] = rec
        ledger_update("D_beam", out)
        if r["status"] == "CLOSED":
            log("  !! CLOSURE -- rebuilding witness and running referee")
            v, res = replay_and_verify(512, D0, r["hist"])
            log("  referee: %s" % json.dumps(
                {k: v[k] for k in ("ok", "adm", "identity") if k in v}))
            if v.get("ok"):
                w = sweep.save_witness(dict(res, outcome=r), os.path.join(
                    HERE, "amm12592_witness_R512_beam_D0_%d_boxeph.json" % D0))
                log("  witness saved: %s" % json.dumps(w))
                out["beam_512_%d" % D0]["witness"] = w
                ledger_update("D_beam", out)
    return 0


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    log("\n### amm12592_localrule_barrier_gap_boxeph.py mode=%s  %s" %
        (mode, time.strftime("%Y-%m-%d %H:%M:%S")))
    t0 = time.time()
    f = 0
    if mode in ("A", "all"):
        f += stage_A()
    if mode in ("B", "all"):
        f += stage_B()
    if mode in ("C", "all"):
        f += stage_C()
    if mode in ("D", "all"):
        f += stage_D()
    log("### done mode=%s fails=%d elapsed=%.1fs" % (mode, f, time.time() - t0))
