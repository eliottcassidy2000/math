"""LANE D1 -- THE BULK RULE: kill eps_inf by design (boxeph, 2026-08-03).

Generalized exact clamp for the T6 clamped-Pascal junk flow.  The proved
conjugacy (T2) holds for ANY admissible in-box choice u_t (even, in the
asymmetric c-box [-2C(d-1,t), +2C(d-1,t-1)]), not only the nearest-point
clamp: the junk passed to the next row is jn_t = w_t - u_t, transported by
the Pascal kernel K = (1,1) (delta=0) or (1,2,1) (delta=1), feed on cells
{0,1} with original 2G_R coefficients, death iff jn_d != 0, capture iff
jn = {} after feed stops.  All decisions exact int; no floats anywhere.

DESIGN (variant "desc", the bulk rule):
  Next-row load at output cell c is  w'_c = jn_c + 2 jn_{c-1} + jn_{c-2}
  (delta'=1; resp. jn_c + jn_{c-1} for delta'=0) + feed' (cells 0,1 only,
  never a finalized output).  Output c is FINALIZED by its lowest
  contributor jn_{c-S}.  Processing cells t descending, when jn_t is chosen
  the contributions from above (jn_{t+1}, jn_{t+2}) are known, so jn_t can
  be chosen inside its exact interval I_t = [w_t - hi_t, w_t - lo_t] to put
  w'_{t+S} INSIDE the next row's box (absorbable => zero junk there), with
  the plain overflow as preference (least intervention) and as fallback
  metric when infeasible (minimize the next overflow).  VENT: two cells
  above the front are allowed nonzero jn to pre-cancel the T9a edge
  transport (the front cell's forced overflow o can only take MORE of the
  same sign; the vent cell above contributes 2*z to the same output, so
  z ~ -o/2 cancels it) -- fires only when the edge transport would overflow
  the next box.  This attacks exactly the march mechanism (T9a/T9'):
  alternation shaping (V1) + kernel lookahead (V2) + venting (V3) in one
  exact local rule.

Certification: variant "plain" (shape off) is asserted trace- and
outcome-identical to the certified fast engine run_fast at
(256,0),(256,1),(512,4),(512,5).  Every run streams two independent
identity spot checks (x=2 and x=3 evaluations of sum x^i Delta_i =
(1-x)^(R-1), computed from the u-ledger, never from the flow state) plus
the ballot-debt count; closures at R <= 1024 are rebuilt into full blocks
and passed to the hostile referee verify_witness (admissibility + full
polynomial identity via t-substitution + ballot laws).

Usage:
  python3 amm12592_bulkrule_design_sweep_boxeph.py regress
  python3 amm12592_bulkrule_design_sweep_boxeph.py run R D0 VARIANT [--verify] [--witness]
  python3 amm12592_bulkrule_design_sweep_boxeph.py scan R D0LIST VARIANT
VARIANT in: plain, desc, descnv (no vent), asc, descw32, descw128, descw512
"""
import sys, os, json, time, io, contextlib, importlib.util, hashlib
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

spec = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(iref)
floorgs = iref.floor_gamma_star

spec2 = importlib.util.spec_from_file_location(
    "fastflow", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fastflow = importlib.util.module_from_spec(spec2)
with contextlib.redirect_stdout(io.StringIO()):
    spec2.loader.exec_module(fastflow)
two_G_coeffs = fastflow.two_G_coeffs

VARIANTS = {
    "plain":   dict(shape=False, vent=0, Wwin=None, asc=False),
    "desc":    dict(shape=True,  vent=2, Wwin=None, asc=False),
    "descnv":  dict(shape=True,  vent=0, Wwin=None, asc=False),
    "asc":     dict(shape=True,  vent=0, Wwin=None, asc=True),
    "descw32": dict(shape=True,  vent=2, Wwin=32,   asc=False),
    "descw128": dict(shape=True, vent=2, Wwin=128,  asc=False),
    "descw512": dict(shape=True, vent=2, Wwin=512,  asc=False),
    "train":   dict(shape=True,  vent=6, Wwin=None, asc=False, guard=48,
                    rho0=2),
    "train16": dict(shape=True,  vent=6, Wwin=None, asc=False, guard=16,
                    rho0=2),
    "trainv2": dict(shape=True,  vent=2, Wwin=None, asc=False, guard=48,
                    rho0=2),
    "train3":  dict(shape=True,  vent=6, Wwin=None, asc=False, guard=48,
                    rho0=3),
    "tr15":    dict(shape=True,  vent=6, Wwin=None, asc=False, guard=16,
                    rho0=(3, 2)),
    "fill":    dict(shape=True,  vent=6, Wwin=None, asc=False, guard=16,
                    rho0=2, fill=True),
    "fill15":  dict(shape=True,  vent=6, Wwin=None, asc=False, guard=16,
                    rho0=(3, 2), fill=True),
    "tr54":    dict(shape=True,  vent=6, Wwin=None, asc=False, guard=16,
                    rho0=(5, 4)),
    "tr43":    dict(shape=True,  vent=6, Wwin=None, asc=False, guard=16,
                    rho0=(4, 3)),
    "tr74":    dict(shape=True,  vent=6, Wwin=None, asc=False, guard=16,
                    rho0=(7, 4)),
    "tr15g48": dict(shape=True,  vent=6, Wwin=None, asc=False, guard=48,
                    rho0=(3, 2)),
    "desc1":   dict(shape=True,  vent=2, Wwin=None, asc=False, vlad=1),
}
for _v in VARIANTS.values():
    _v.setdefault("vlad", None)
for _v in VARIANTS.values():
    _v.setdefault("fill", False)
for _v in VARIANTS.values():
    _v.setdefault("guard", None)
    _v.setdefault("rho0", 2)


def t4_row_load(R, d):
    """Dense row-0 load w_t (T4 closed form), list length d+1. Exact."""
    A = comb(R - 2, d)            # C(R-2-t, d-t)
    B = d + 1                     # C(d+1, t+1)
    Cd = 1                        # C(d, t)
    sgn = 1 if d % 2 == 0 else -1
    w = []
    for t in range(d + 1):
        w.append(sgn * A - B + 2 * Cd)
        if t < d:
            A = A * (d - t) // (R - 2 - t)
            B = B * (d - t) // (t + 2)
            Cd = Cd * (d - t) // (t + 1)
            sgn = -sgn
    assert w[d] == 2
    return w


class DescBinom:
    """C(n,k) with k stepping DOWN by 1; exact; 0 outside [0,n]."""
    __slots__ = ("n", "k", "v")

    def __init__(self, n, k):
        self.n, self.k = n, k
        self.v = comb(n, k) if 0 <= k <= n else 0

    def step(self):
        n, k = self.n, self.k
        self.k = k - 1
        if self.k < 0 or self.k > n:
            self.v = 0
        elif self.k == n:
            self.v = 1
        else:
            self.v = self.v * k // (n - k + 1)
        return self.v


class AscBinom:
    """C(n,k) with k stepping UP by 1; exact; 0 outside [0,n]."""
    __slots__ = ("n", "k", "v")

    def __init__(self, n, k):
        self.n, self.k = n, k
        self.v = comb(n, k) if 0 <= k <= n else 0

    def step(self):
        n, k = self.n, self.k
        self.k = k + 1
        if self.k < 0 or self.k > n:
            self.v = 0
        elif self.k == 0:
            self.v = 1
        else:
            self.v = self.v * (n - k) // (k + 1)
        return self.v


def choose_row(w, d, dn, S, cfg, record_u):
    """Generalized clamp of one row.  w: dict t->even value, support <= d.
    Returns (jn, us, c0, junkL1, diag).  us=None unless record_u."""
    shape, vent, Wwin, asc = cfg["shape"], cfg["vent"], cfg["Wwin"], cfg["asc"]
    ts = sorted(w)
    ta, tb = ts[0], ts[-1]
    assert tb <= d, ("support beyond bottom cell", tb, d)
    thi = min(tb + (vent if shape else 0), d)
    tlo = max(0, ta - (2 if shape else 0)) if shape else ta
    tcut = tlo if Wwin is None else max(tlo, tb - Wwin)
    jn = {}
    us = {} if record_u else None
    c0 = 0
    junkL1 = 0
    nfail = 0
    fail_cells = []
    nshaped = 0
    # vent preference ladder: pre-cancel the front cell's unavoidable edge
    # transport with an alternating halving chain above the front
    vent_prefs = {}
    if shape and vent and thi > tb:
        loF = -2 * comb(d - 1, tb) if tb <= d - 1 else 0
        hiF = 2 * comb(d - 1, tb - 1) if 1 <= tb <= d else 0
        vF = w[tb]
        oF = vF - min(hiF, max(loF, vF))          # plain overflow at front
        out = tb + S                              # its edge-transport output
        loN = -2 * (comb(dn - 1, out) if out <= dn - 1 else 0)
        hiN = 2 * (comb(dn - 1, out - 1) if 1 <= out <= dn else 0)
        excess = oF - min(hiN, max(loN, oF))      # unabsorbable part
        if excess:
            _vl = cfg.get("vlad")
            for k in range(1, min(thi - tb, _vl or (thi - tb)) + 1):
                p = (-1) ** k * (excess if S == 1 else (excess >> k))
                pe = 2 * (p // 2)
                if pe == 0:
                    break
                vent_prefs[tb + k] = pe
    if asc:
        rng = range(tlo, thi + 1)
        Ab = AscBinom(d - 1, tlo)          # C(d-1, t)      -> lo = -2*
        Bb = AscBinom(d - 1, tlo - 1)      # C(d-1, t-1)    -> hi = +2*
        An = AscBinom(dn - 1, tlo)         # out = t
        Bn = AscBinom(dn - 1, tlo - 1)
    else:
        rng = range(thi, tlo - 1, -1)
        Ab = DescBinom(d - 1, thi)
        Bb = DescBinom(d - 1, thi - 1)
        An = DescBinom(dn - 1, thi + S)    # out = t + S
        Bn = DescBinom(dn - 1, thi + S - 1)
    j1 = j2 = 0        # junk at the two previously processed cells
    first = True
    for t in rng:
        if not first:
            Ab.step(); Bb.step(); An.step(); Bn.step()
        first = False
        v = w.get(t, 0)
        lo, hi = -2 * Ab.v, 2 * Bb.v
        up = v if lo <= v <= hi else (lo if v < lo else hi)
        oplain = v - up
        j = oplain
        if shape and t < d and t >= tcut:
            inc = (2 * j1 + j2) if S == 2 else j1
            pref = vent_prefs.get(t, oplain) if (not asc and t > tb) \
                else oplain
            flo = max(v - hi, -2 * An.v - inc)
            fhi = min(v - lo, 2 * Bn.v - inc)
            if flo <= fhi:
                j = flo if pref < flo else (fhi if pref > fhi else pref)
            else:
                nfail += 1
                if len(fail_cells) < 6:
                    fail_cells.append(t)
                if (v - lo) < -2 * An.v - inc:
                    j = v - lo            # I entirely below target: take max
                elif (v - hi) > 2 * Bn.v - inc:
                    j = v - hi            # I entirely above target: take min
            # TRAIN GUARD (front layer): keep the junk train alternating
            # with envelope ratio >= rho0 (plateau/hole kill); boosts only
            # within I_t, never flips a forced sign.
            K = cfg.get("guard")
            if K is not None and t <= tb and t >= tb - K and not asc:
                r0 = cfg["rho0"]
                num, den = r0 if isinstance(r0, tuple) else (r0, 1)
                if j1:
                    tsgn = 1 if j1 < 0 else -1
                    m = 2 * ((num * abs(j1)) // (2 * den))
                elif j2:
                    tsgn = 1 if j2 > 0 else -1
                    m = 2 * ((num * num * abs(j2)) // (2 * den * den))
                else:
                    tsgn, m = 0, 0
                if j and tsgn and ((j > 0) == (tsgn > 0)) and abs(j) < m:
                    jb = min(m, v - lo) if j > 0 else max(-m, v - hi)
                    if jb != j:
                        j = jb
                elif (not j) and tsgn and m and cfg.get("fill"):
                    # hole fill: keep alternation-phased junk instead of
                    # full absorption inside the train
                    jb = tsgn * m
                    jb = min(max(jb, v - hi), v - lo)
                    if (jb > 0) == (tsgn > 0) and jb:
                        j = jb
        if j != oplain:
            nshaped += 1
        u = v - j
        assert lo <= u <= hi and not (u & 1), ("clamp", t, d)
        if j:
            jn[t] = j
            junkL1 += abs(j)
        if u and us is not None:
            us[t] = u
        if t == 0:
            c0 = u
        j2, j1 = j1, j
    return jn, us, c0, junkL1, {"nfail": nfail, "fails": fail_cells,
                                "nshaped": nshaped}


def run_bulk(R, D0, variant="desc", keep_rows=True, collect_u=False,
             spot=True, trace_path=None, flush_every=256):
    """Exact generalized-rule run.  Returns result dict (mirrors run_fast)."""
    cfg = VARIANTS[variant]
    t0 = time.time()
    g = two_G_coeffs(R)
    trace = []
    urows = {}
    SP2 = 0
    SP3 = 0
    minus2 = 0
    j = {}
    d_prev = None
    outcome = None
    capture = None
    front0 = None
    T6b = None

    def rowrec(i, d, jj, junkL1, c0, wb, diag):
        tsj = sorted(jj)
        r = {"i": i, "d": d, "nclamped": len(jj),
             "tmin": (tsj[0] if tsj else None),
             "tmax": (tsj[-1] if tsj else None),
             "junkL1_bits": junkL1.bit_length(), "c0": c0, "e_in0": wb,
             "gap": (d - tsj[-1]) if tsj else None}
        if diag is not None:
            r["nfail"] = diag["nfail"]
            r["nshaped"] = diag["nshaped"]
            if diag["fails"]:
                r["fail_cells"] = diag["fails"]
        return r

    for i in range(0, R - 1):
        d = floorgs(R + i) + D0
        dn = floorgs(R + i + 1) + D0
        S = 1 + (dn - d)
        if i == 0:
            w = {t: v for t, v in enumerate(t4_row_load(R, d)) if v}
            fed = True
        else:
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
            w = {t: v for t, v in w.items() if v}
            fed = False
            if d + i <= R - 1:
                w[0] = w.get(0, 0) + g[d + i]
                fed = True
            if delta == 1 and d - 1 + i <= R - 1:
                w[0] = w.get(0, 0) + g[d - 1 + i]
                w[1] = w.get(1, 0) + g[d - 1 + i]
                fed = True
            w = {t: v for t, v in w.items() if v}
            if not w:
                jn, us, c0, junkL1, diag = {}, ({} if (collect_u or spot)
                                                else None), 0, 0, None
        w_bottom = w.get(d, 0)
        if w:
            jn, us, c0, junkL1, diag = choose_row(
                w, d, dn, S, cfg, record_u=(collect_u or spot))
        if spot and us:
            p2 = 0
            p3 = 0
            for t in sorted(us, reverse=True):
                u = us[t]
                sg2 = -1 if (t & 1) else 1
                p2 += u * sg2 * (1 << (d - t))
                p3 += u * ((-2) ** t) * (3 ** (d - t))
            SP2 += (3 + p2) << i
            SP3 += (5 + p3) * 3 ** i
        elif spot:
            SP2 += 3 << i
            SP3 += 5 * 3 ** i
        if collect_u and us:
            urows[i] = us
        if i == 0:
            front0 = max(jn) if jn else -1
            T6b = d - front0
        died = d in jn
        if died:
            outcome = {"status": "DIE", "row": i,
                       "const_bits": abs(jn[d]).bit_length()}
            r = rowrec(i, d, jn, junkL1, c0, w_bottom, diag)
            r["death_const_bits"] = abs(jn[d]).bit_length()
            trace.append(r)
        else:
            if c0 == -2:
                minus2 += 1
            if keep_rows or (i % 8 == 0) or not jn:
                trace.append(rowrec(i, d, jn, junkL1, c0, w_bottom, diag))
        j = jn
        d_prev = d
        if died:
            break
        if not j and not fed and d + i > R - 1:
            capture = i
            outcome = {"status": "CLOSED", "capture_row": i}
            break
        if trace_path and i % flush_every == 0:
            json.dump({"R": R, "D0": D0, "variant": variant, "partial_row": i,
                       "trace": trace,
                       "elapsed_s": round(time.time() - t0, 1)},
                      open(trace_path, "w"))
    if outcome is None:
        if not j:
            outcome = {"status": "CLOSED", "capture_row": R - 2}
            capture = R - 2
        else:
            from amm12592_allR_family_toolbox_boxeph import (bern_to_poly,
                                                             poly_to_bern)
            P = bern_to_poly([j.get(t, 0) for t in range(d_prev + 1)], d_prev)
            assert P and P[0] == 0
            P = P[1:]
            dL = floorgs(R + R - 1) + D0
            wl = poly_to_bern(P + [0] * (dL + 1 - len(P)), dL) \
                if len(P) <= dL + 1 else None
            ok = wl is not None and all(
                0 <= wl[t] <= 2 * comb(dL, t) and wl[t] % 2 == 0
                for t in range(dL + 1))
            outcome = ({"status": "CLOSED", "capture_row": None,
                        "lastrow_absorbed": True} if ok else
                       {"status": "OPEN_RESIDUAL"})
    # streamed identity spots (only meaningful for CLOSED-with-capture runs)
    spots = None
    if spot and outcome["status"] == "CLOSED" and capture is not None:
        for i in range(capture + 1, R - 1):
            SP2 += 3 << i
            SP3 += 5 * 3 ** i
        SP2 += -(1 << (R - 1))          # last row Delta = -1
        SP3 += -(3 ** (R - 1))
        spots = {"x2": SP2 == (-1) ** (R - 1),
                 "x3": SP3 == (-2) ** (R - 1),
                 "debt": minus2 == (R - 2) // 2}
    res = {"R": R, "D0": D0, "variant": variant, "outcome": outcome,
           "capture_row": capture, "minus2_rows": minus2, "front0": front0,
           "T6b_bound": T6b, "spots": spots,
           "elapsed_s": round(time.time() - t0, 1), "trace": trace}
    if collect_u:
        res["_urows"] = urows
    if trace_path:
        json.dump({k: v for k, v in res.items() if k != "_urows"},
                  open(trace_path, "w"))
    return res


# ------------------------------------------------------------- verification
def build_blocks(R, D0, urows):
    from amm12592_allR_family_toolbox_boxeph import ballot
    prof = [floorgs(R + i) + D0 for i in range(R)]
    blocks = []
    for i in range(R - 1):
        d = prof[i]
        b = ballot(d)
        u = urows.get(i)
        if u:
            b = b[:]
            for t, v in u.items():
                b[t] += v
        blocks.append(b)
    dL = prof[R - 1]
    blocks.append([-comb(dL, t) for t in range(dL + 1)])
    return prof, blocks


def referee_verify(res):
    """Rebuild full blocks from the u-ledger and run the hostile referee."""
    spec3 = importlib.util.spec_from_file_location(
        "aref", os.path.join(HERE, "amm12592_allR_referee_boxeph.py"))
    aref = importlib.util.module_from_spec(spec3)
    with contextlib.redirect_stdout(io.StringIO()):
        spec3.loader.exec_module(aref)
    prof, blocks = build_blocks(res["R"], res["D0"], res["_urows"])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        v = aref.verify_witness("bulk_%s_R%d_D0%d" % (res["variant"],
                                                      res["R"], res["D0"]),
                                res["R"], res["D0"], prof, blocks)
    return v, buf.getvalue()


def save_witness(res, path):
    urows = res["_urows"]
    doc = {"R": res["R"], "D0": res["D0"], "variant": res["variant"],
           "rule": "bulk generalized clamp (see amm12592_bulkrule_design_"
                   "sweep_boxeph.py); Delta_i = ballot(d_i) + sparse_c[i], "
                   "rows absent are pure ballot (2x-1), last row = corner -1",
           "profile": "floor_gamma_star(R+i)+D0",
           "capture_row": res["capture_row"],
           "sparse_c": {str(i): {str(t): v for t, v in row.items()}
                        for i, row in urows.items()}}
    s = json.dumps(doc)
    sha = hashlib.sha256(s.encode()).hexdigest()
    if len(s) <= 40 * 10 ** 6:
        open(path, "w").write(s)
        return {"path": path, "bytes": len(s), "sha256": sha}
    return {"path": None, "bytes": len(s), "sha256": sha,
            "note": "over 40MB cap; sha of canonical json recorded"}


def regress():
    """Variant 'plain' must reproduce the certified engine bit-for-bit."""
    keys = ("d", "nclamped", "tmin", "tmax", "junkL1_bits", "c0", "e_in0")
    allok = True
    for R, D0 in ((256, 0), (256, 1), (512, 4), (512, 5)):
        a = fastflow.run_fast(R, D0)
        b = run_bulk(R, D0, "plain", spot=True)
        atr = {r["i"]: r for r in a["trace"]}
        btr = {r["i"]: r for r in b["trace"]}
        bad = []
        for i in sorted(set(atr) & set(btr)):
            for k in keys:
                if atr[i].get(k) != btr[i].get(k):
                    bad.append((i, k, atr[i].get(k), btr[i].get(k)))
        ok = (a["outcome"] == b["outcome"] and a["minus2_rows"] ==
              b["minus2_rows"] and a["front0"] == b["front0"] and not bad)
        allok = allok and ok
        print("REGRESS R=%d D0=%d: outcome_eq=%s minus2_eq=%s front0_eq=%s "
              "rows=%d mism=%d spots=%s -> %s"
              % (R, D0, a["outcome"] == b["outcome"],
                 a["minus2_rows"] == b["minus2_rows"],
                 a["front0"] == b["front0"], len(set(atr) & set(btr)),
                 len(bad), b["spots"], "OK" if ok else "FAIL"), flush=True)
        if bad:
            print("   first mismatches:", bad[:5], flush=True)
    print("REGRESSION:", "ALL OK" if allok else "FAILED", flush=True)
    return allok


def summarize(res):
    o = res["outcome"]
    tr = res["trace"]
    maxfront = max((r["tmax"] or -1) for r in tr)
    mingap = min((r["gap"] for r in tr if r["gap"] is not None), default=None)
    lastL1 = tr[-1]["junkL1_bits"]
    peakL1 = max(r["junkL1_bits"] for r in tr)
    nfail_tot = sum(r.get("nfail", 0) for r in tr)
    nsh_tot = sum(r.get("nshaped", 0) for r in tr)
    return ("%s R=%d D0=%d: %s  T6b=%s cap=%s minus2=%s front0=%s "
            "maxfront=%d mingap=%s peakL1b=%d failrows=%d shaped=%d "
            "spots=%s (%.1fs)"
            % (res["variant"], res["R"], res["D0"], o, res["T6b_bound"],
               res["capture_row"], res["minus2_rows"], res["front0"],
               maxfront, mingap, peakL1, nfail_tot, nsh_tot, res["spots"],
               res["elapsed_s"]))


OUT = os.path.join(RESULTS, "amm12592_bulkrule_design_sweep_boxeph.out")


def log(line):
    print(line, flush=True)
    with open(OUT, "a") as f:
        f.write(line + "\n")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "regress":
        ok = regress()
        log("regression plain-vs-engine (256/512, D0 0/1/4/5): %s"
            % ("ALL OK" if ok else "FAILED"))
        sys.exit(0 if ok else 1)
    elif mode == "run":
        R, D0, var = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
        verify = "--verify" in sys.argv
        wit = "--witness" in sys.argv
        tp = os.path.join(RESULTS, "amm12592_bulkrule_trace_R%d_D0%d_%s_"
                                   "boxeph.json" % (R, D0, var))
        res = run_bulk(R, D0, var, keep_rows=(R <= 2048),
                       collect_u=(verify or wit), trace_path=tp)
        log(summarize(res))
        if verify and res["outcome"]["status"] == "CLOSED" and \
                res["capture_row"] is not None:
            v, txt = referee_verify(res)
            log("REFEREE %s R=%d D0=%d: ok=%s adm=%s identity=%s debt=%s"
                % (var, R, D0, v.get("ok"), v.get("adm"), v.get("identity"),
                   v.get("debt")))
        if wit and res["outcome"]["status"] == "CLOSED":
            wp = os.path.join(HERE, "amm12592_witness_R%d_bulk_%s_D0_%d_"
                                    "boxeph.json" % (R, var, D0))
            info = save_witness(res, wp)
            log("WITNESS %s" % info)
    elif mode == "scan":
        R = int(sys.argv[2])
        var = sys.argv[4] if len(sys.argv) > 4 else "desc"
        for D0 in [int(v) for v in sys.argv[3].split(",")]:
            tp = os.path.join(RESULTS, "amm12592_bulkrule_trace_R%d_D0%d_%s_"
                                       "boxeph.json" % (R, D0, var))
            res = run_bulk(R, D0, var, keep_rows=False, trace_path=tp)
            log(summarize(res))


def nu_profile(R, D0, variant="plain", rows=None):
    """Exact anatomy: per row, the pairwise-sum residue nu = (1+sigma)j
    (nonalternating content), its band profile, and front-envelope ratios.
    Pure measurement (plain flow unless variant says else)."""
    cfg = VARIANTS[variant]
    g = two_G_coeffs(R)
    j = {}
    d_prev = None
    out = []
    for i in range(0, R - 1):
        d = floorgs(R + i) + D0
        dn = floorgs(R + i + 1) + D0
        S = 1 + (dn - d)
        if i == 0:
            w = {t: v for t, v in enumerate(t4_row_load(R, d)) if v}
        else:
            delta = d - d_prev
            w = {}
            K = (1, 1) if delta == 0 else (1, 2, 1)
            for t, v in j.items():
                for s, ks in enumerate(K):
                    w[t + s] = w.get(t + s, 0) + ks * v
            if d + i <= R - 1:
                w[0] = w.get(0, 0) + g[d + i]
            if delta == 1 and d - 1 + i <= R - 1:
                w[0] = w.get(0, 0) + g[d - 1 + i]
                w[1] = w.get(1, 0) + g[d - 1 + i]
            w = {t: v for t, v in w.items() if v}
        if not w:
            break
        jn, us, c0, junkL1, diag = choose_row(w, d, dn, S, cfg, record_u=False)
        died = d in jn
        if jn and (rows is None or i in rows):
            T = max(jn)
            # nu = pairwise sums (nonalternating content), band L1 bits
            nu = {}
            for t in sorted(jn):
                nu[t] = jn.get(t, 0) + jn.get(t - 1, 0)
                nu[t + 1] = nu.get(t + 1, 0) + 0  # ensure boundary term seen
            for t in list(nu):
                nu[t] = jn.get(t, 0) + jn.get(t - 1, 0)
            nbands = 8
            bl = [0] * nbands
            jl = [0] * nbands
            for t, v in nu.items():
                if 0 <= t <= T:
                    bl[min(nbands - 1, t * nbands // (T + 1))] += abs(v)
            for t, v in jn.items():
                jl[min(nbands - 1, t * nbands // (T + 1))] += abs(v)
            # front envelope: |j| at T, T-1, ..., T-6 (bits) + exact ratios
            env = []
            for k in range(7):
                v = jn.get(T - k, 0)
                env.append(v.bit_length() if v >= 0 else -((-v).bit_length()))
            rec = {"i": i, "d": d, "front": T, "gap": d - T,
                   "junkL1_bits": junkL1.bit_length(),
                   "nu_bits": [b.bit_length() for b in bl],
                   "j_bits": [b.bit_length() for b in jl],
                   "env_sgnbits": env, "died": died}
            out.append(rec)
        j = jn
        d_prev = d
        if died or not j:
            break
    return out


def cmd_nu():
    import sys as _s
    R, D0 = int(_s.argv[2]), int(_s.argv[3])
    var = _s.argv[4] if len(_s.argv) > 4 else "plain"
    rows = set(range(0, 130)) | set(range(130, 400, 5))
    prof = nu_profile(R, D0, var, rows)
    p = os.path.join(RESULTS, "amm12592_bulkrule_nuprof_R%d_D0%d_%s_boxeph.json"
                     % (R, D0, var))
    json.dump(prof, open(p, "w"))
    for r in prof[:0]:
        pass
    print("saved", p, len(prof), "rows", flush=True)


if len(sys.argv) > 1 and sys.argv[1] == "nu":
    cmd_nu()


def block_ledger(R, D0, variant="plain", maxrows=None):
    """Exact per-row measurement for the immortality conjecture: find the
    longest same-sign run within the 40 cells behind the front; record its
    length, leading (deepest) cell, min/max |value| inside the run, and the
    exact geometric cap tail sum_{t' > lead} 2 C(d-1, t') (the total box
    budget strictly beyond the run's lead cell).  Reports first row where
    min|run value| > cap tail (candidate point of no return)."""
    cfg = VARIANTS[variant]
    g = two_G_coeffs(R)
    j = {}
    d_prev = None
    out = {"R": R, "D0": D0, "variant": variant, "rows": [],
           "first_super": None, "outcome": None}
    for i in range(0, R - 1):
        d = floorgs(R + i) + D0
        dn = floorgs(R + i + 1) + D0
        S = 1 + (dn - d)
        if i == 0:
            w = {t: v for t, v in enumerate(t4_row_load(R, d)) if v}
        else:
            delta = d - d_prev
            K = (1, 1) if delta == 0 else (1, 2, 1)
            w = {}
            for t, v in j.items():
                for s2, ks in enumerate(K):
                    w[t + s2] = w.get(t + s2, 0) + ks * v
            if d + i <= R - 1:
                w[0] = w.get(0, 0) + g[d + i]
            if delta == 1 and d - 1 + i <= R - 1:
                w[0] = w.get(0, 0) + g[d - 1 + i]
                w[1] = w.get(1, 0) + g[d - 1 + i]
            w = {t: v for t, v in w.items() if v}
        if not w:
            out["outcome"] = {"CLOSED_capture": i}
            break
        jn, us, c0, junkL1, diag = choose_row(w, d, dn, S, cfg, record_u=False)
        died = d in jn
        if jn:
            T = max(jn)
            # longest same-sign run in cells [T-40, T]
            best = (0, None, None, None)   # len, lead, minv, maxv
            runlen = 0
            runmin = runmax = None
            lead = None
            prev_sign = 0
            for t in range(max(0, T - 40), T + 1):
                v = jn.get(t, 0)
                sg = 0 if v == 0 else (1 if v > 0 else -1)
                if sg != 0 and sg == prev_sign:
                    runlen += 1
                    runmin = min(runmin, abs(v))
                    runmax = max(runmax, abs(v))
                    lead = t
                elif sg != 0:
                    runlen = 1
                    runmin = runmax = abs(v)
                    lead = t
                else:
                    runlen = 0
                prev_sign = sg
                if runlen > best[0]:
                    best = (runlen, lead, runmin, runmax)
            L, lead, mn, mx = best
            rec = {"i": i, "front": T, "gap": d - T, "runlen": L}
            if L >= 2 and lead is not None:
                tail = 0
                c = 2 * comb(d - 1, lead)
                for tt in range(lead + 1, d + 1):
                    c = c * (d - tt) // tt if tt <= d - 1 else 0
                    tail += c
                rec["lead"] = lead
                rec["minv_bits"] = mn.bit_length()
                rec["maxv_bits"] = mx.bit_length()
                rec["captail_bits"] = tail.bit_length()
                rec["super"] = mn > tail
                if rec["super"] and out["first_super"] is None:
                    out["first_super"] = i
            if maxrows is None or i < maxrows or died:
                out["rows"].append(rec)
        j = jn
        d_prev = d
        if died:
            out["outcome"] = {"DIE": i}
            break
        if not j and dn + i > R - 1:
            out["outcome"] = {"CLOSED_capture": i}
            break
    return out


if len(sys.argv) > 1 and sys.argv[1] == "block":
    R, D0 = int(sys.argv[2]), int(sys.argv[3])
    var = sys.argv[4] if len(sys.argv) > 4 else "plain"
    res = block_ledger(R, D0, var)
    pth = os.path.join(RESULTS, "amm12592_bulkrule_blockledger_R%d_D0%d_%s_"
                       "boxeph.json" % (R, D0, var))
    json.dump(res, open(pth, "w"))
    sup = [r for r in res["rows"] if r.get("super")]
    long_runs = [r for r in res["rows"] if r["runlen"] >= 5]
    print("R=%d D0=%d %s: outcome=%s first_super_row=%s n_super_rows=%d "
          "first_run>=5_row=%s" % (R, D0, var, res["outcome"],
          res["first_super"], len(sup),
          long_runs[0]["i"] if long_runs else None), flush=True)
