#!/usr/bin/env python3
"""THM-2365 (10) deep-coloured a != 0 mode lines of the bare tensor H_E on
the canonical typed row: walk endpoints and exact phase-cone feasibility.

EXTENDS (does not rebuild) the gated a = 0 sidecar
    04-computation/lrc14_phase_cone_sidecar_opus_20260727.py
whose chamber machinery, interval machinery, exact cyclotomic sign machinery
and gates are IMPORTED and re-run verbatim (module import; the sidecar's
main() is not executed).

OBJECT.  Canonical typed row THM-2309 (25):
    w = (1,14,27,40,53,66,13,2197,742586), owner j = 1, word sigma = {a},
    grafts k_a = q1 = 14, k_c = q2 = 27 (both = 1 mod 13),
    blocker scales A_j = 1, A_a = 169, A_c = 57122 (THM-2368 (8)).

MODES (THM-2365 (10), exact convention log):
    B(a,b,h) = (1/13^3) sum_{r,s,t} H_E(r,s,t) zeta^{a r + b s + h t},
    zeta = exp(2 pi i/13) (standard embedding; Galois covered as in the
    sidecar).  The walk endpoint computed here is, per mode (a,lam,mu),

      endpoint(a,lam,mu) = sum_{r,s,t} Hnum(r,s,t) zeta^{a r - lam s - mu t}
                         = 13^4 * T_DEN * B(a, -lam, -mu),

    Hnum = 13 * T_DEN * H_E (the sidecar's gated integer tensor), matching
    the sidecar's stated a = 0 identity endpoint = 13^4*T_DEN*B(0,-lam,-mu).
    Per-chamber contribution (common chamber base BEFORE any aggregation,
    MISTAKE-281 discipline; r, s, t all live on ONE chamber partition):

      c_i(a,lam,mu) = len_i * phi_i(a) * V_i(lam,mu)   in Z[zeta_13],
      phi_i(a)      = sum_r dc_i[r] zeta^{a r}         (probe word DFT),
      V_i(lam,mu)   = sum_{s,t} ga_i[s] gc_i[t] N_i[s,t] zeta^{-(lam s+mu t)}.

    For a = 0 this degenerates to the sidecar's scalar weight nd_i = phi_i(0)
    and reproduces its walk exactly (reproduction-gated below).

EXACT INDEX SET (logged; AMBIGUITY in the task count resolved):
    By THM-2365 (13) the target vector of B(a,b,h) is q = (b, a+h); with
    b = -lam, h = -mu the target of endpoint(a,lam,mu) is (-lam, a-mu).
    The ZERO-target mode on line a is therefore (lam,mu) = (0,a), NOT (0,0):
    each a-line has 168 nonzero-target modes {(lam,mu) != (0,a)}.
      * new deep-coloured a != 0 nonzero-target modes: 12 x 168 = 2016;
      * full nonzero-target mode set over all a in F_13: 13 x 168 = 2184
        ( = 12 x 182, the task's count -- it is the TOTAL including the 168
        a = 0 modes already decided by the sidecar; a per-line count of 182
        is not the cardinality of any F_13^2 index subset);
      * zero-target modes (one per line, on the projection-retained line
        b = 0, a+h = 0): 13; these are NOT drift modes (THM-2365 (16)-(17))
        and are decided/report-flagged separately, never counted as
        cone-explainable positivity.
    ALL 13 x 169 = 2197 modes are computed; the a = 0 line must reproduce
    the stored sidecar results (gate); decision value is claimed on the
    2016 new deep-coloured nonzero-target modes.

DECISION RULE (exact arithmetic in every decision):
    A mode's phase-cone certificate HOLDS iff its aggregated primitive
    direction classes in Z[zeta_13] (positive-rational-multiple classes,
    identical reduction to the sidecar's) fit a closed half-plane with at
    least one direction off the boundary line (sidecar classes '<pi' and
    '=pi-cert'); it FAILS on '=pi-fail' / '>pi'.  Floats are used ONLY to
    route which exact certificate is attempted:
      * clear-failure route: an exact positive-spanning triple u,v,w with
        cross(u,v) > 0, cross(v,w) > 0, cross(w,u) > 0 (exact cyclotomic
        cross signs; the 2-D identity cross(v,w)u + cross(w,u)v +
        cross(u,v)w = 0 then writes 0 as a strictly positive combination,
        so no closed half-plane contains the set: class '>pi', exact);
      * all other routes run the sidecar's full exact cone test verbatim.
    First cancelling events are located by float-accelerated scan and then
    EXACTLY certified by bisection with the full exact cone test on prefix
    sets (feasibility is monotone under adding directions, so the minimal
    infeasible prefix is well-defined and the bisection is valid).
    Endpoint zero tests are algebraic (all-coordinates-equal in Z^13).
    Floats appear only in descriptive statistics, marked DESCRIPTIVE.

GATES (all re-run, none weakened):
    (A) full 13^3 bare-table equality against the audited drift-companion
        interval machinery (verbatim, all 169 twists);
    (B) row sums == byte-stored S_E numerators of the drift companion;
    (C) diagonal-plane zero H_E(t,s,t) = 0, all 169;
    (D) per-mode endpoint equality: state-product accumulation == direct
        exact DFT binning of the gated Hnum tensor (all 2197 modes);
    (E) THM-2365 (14) target-null circulation: for every target (b,q2),
        sum_{a in F_13} endpoint(a,-b mod 13,(a-q2) mod 13) = 0 exactly;
    (F) a = 0 reproduction: class counts, endpoint-nonzero count, direction
        statistics and the first cancelling event must equal the stored
        sidecar artifact values (anchors hard-coded from the .out).

SCOPE (MISTAKE-281/283 in force): the row is TYPED, not an asserted scalar
cover; no physical current; no scalar row excluded; no LRC(14)/JC(2)
progress beyond what the cited theorems license; full support is never read
as noncancellation; pairing happens on one common chamber base BEFORE
marginalization; every convention is logged above or printed below.

Script:  04-computation/lrc14_deep_mode_cone_opus_20260728.py
Output:  05-knowledge/results/lrc14_deep_mode_cone_opus_20260728.out
"""

import importlib.util
import os
import sys
import time
from bisect import insort
from fractions import Fraction
from functools import cmp_to_key
from math import atan2, cos, pi, sin

import numpy as np

T0 = time.time()


def log(msg=""):
    print(msg, flush=True)


def el():
    return f"[{time.time() - T0:8.1f}s]"


HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SIDECAR_PATH = os.path.join(HERE, "lrc14_phase_cone_sidecar_opus_20260727.py")
_spec = importlib.util.spec_from_file_location("sidecar_20260727",
                                               SIDECAR_PATH)
sc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sc)

require = sc.require
T_DEN = sc.T_DEN
D2 = sc.D2

# ---------------------------------------------------------------------------
# Sidecar cone machinery (verbatim copies of the closures inside its main();
# all exact-sign primitives are the imported sidecar functions)
# ---------------------------------------------------------------------------


def sort_dirs(dlist):
    def half(d):
        i = sc.im_sign(d)
        if i > 0:
            return 0
        if i < 0:
            return 1
        return 0 if sc.re_sign(d) > 0 else 1

    halves = {d: half(d) for d in dlist}

    def cmp(u, v):
        hu, hv = halves[u], halves[v]
        if hu != hv:
            return -1 if hu < hv else 1
        s = sc.cross_sign(list(u), list(v))
        if s > 0:
            return -1
        if s < 0:
            return 1
        return 0
    return sorted(dlist, key=cmp_to_key(cmp))


def cone_test(dlist):
    """Verbatim sidecar semantics.  Returns (span_class, info)."""
    if not dlist:
        return "empty", None
    if len(dlist) == 1:
        return "<pi", (dlist[0], dlist[0])
    sd = sort_dirs(dlist)
    m = len(sd)
    big_gap_at = None
    pi_gaps = []
    for i in range(m):
        u = sd[i]
        v = sd[(i + 1) % m]
        cs = sc.cross_sign(list(u), list(v))
        if cs < 0:
            big_gap_at = i
            break
        if cs == 0 and sc.dot_sign(list(u), list(v)) < 0:
            pi_gaps.append(i)
    if big_gap_at is not None:
        i = big_gap_at
        return "<pi", (sd[(i + 1) % m], sd[i])
    if pi_gaps:
        i = pi_gaps[0]
        u = sd[i]
        if any(sc.cross_sign(list(u), list(d)) != 0 for d in sd):
            return "=pi-cert", (sd[(i + 1) % m], sd[i])
        return "=pi-fail", None
    return ">pi", None


INFEASIBLE = ("=pi-fail", ">pi")
FEASIBLE = ("<pi", "=pi-cert")

# float embedding (ROUTING + DESCRIPTIVE only; never decides)
EMBC = np.array([cos(2 * pi * j / 13) for j in range(13)])
EMBS = np.array([sin(2 * pi * j / 13) for j in range(13)])
BAND = 1e-6


def float_maxgap(dirs_np):
    ang = np.arctan2(dirs_np @ EMBS, dirs_np @ EMBC)
    ang = np.sort(ang)
    if ang.size == 1:
        return 2 * pi, ang
    gaps = np.diff(ang)
    wrap = ang[0] + 2 * pi - ang[-1]
    return max(float(gaps.max()), float(wrap)), ang


def cert3(dirs_np):
    """Try to certify class '>pi' EXACTLY via a positive-spanning triple.
    Returns True on success (exact), False when no triple was certified
    (caller must fall back to the full exact cone test)."""
    ang = np.arctan2(dirs_np @ EMBS, dirs_np @ EMBC)
    order = np.argsort(ang, kind="stable")
    n = order.size
    sh = ang[order] - ang[order[0]]
    j = int(np.searchsorted(sh, pi)) - 1
    u = [int(t) for t in dirs_np[order[0]]]
    for dj in range(4):
        vi = j - dj
        if vi <= 0:
            break
        v = [int(t) for t in dirs_np[order[vi]]]
        if sc.cross_sign(u, v) <= 0:
            continue
        for dw in range(1, 5):
            wi = j + dw
            if wi >= n:
                break
            w = [int(t) for t in dirs_np[order[wi]]]
            if sc.cross_sign(v, w) > 0 and sc.cross_sign(w, u) > 0:
                return True
    return False


def decide_mode(dirs_np):
    """Exact class of one mode's direction set.
    Returns (cls, info, ndirs, span_float_or_None, route)."""
    n = int(dirs_np.shape[0])
    if n == 0:
        return "empty", None, 0, None, "trivial"
    if n == 1:
        d = tuple(int(v) for v in dirs_np[0])
        return "<pi", (d, d), 1, 0.0, "trivial"
    mg, _ = float_maxgap(dirs_np)
    span = 2 * pi - mg
    if mg < pi - BAND and cert3(dirs_np):
        return ">pi", None, n, span, "cert3"
    dl = [tuple(int(v) for v in row) for row in dirs_np]
    cls, info = cone_test(dl)
    return cls, info, n, span, "full"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log("THM-2365 (10) deep-coloured a != 0 mode lines: walk endpoints and"
        " exact phase-cone feasibility,")
    log("bare tensor H_E, canonical typed row (extension of the gated a = 0"
        " sidecar)")
    log("script=04-computation/lrc14_deep_mode_cone_opus_20260728.py")
    log("machine=opus  date=2026-07-28")
    log("")
    log("[1] conventions and exact mode index set (module docstring is the"
        " full log; files override summaries)")
    log("    w=(1, 14, 27, 40, 53, 66, 13, 2197, 742586); grafts k_a=14,"
        " k_c=27; A_j=1, A_a=169, A_c=57122")
    log("    endpoint(a,lam,mu) = sum Hnum(r,s,t) zeta^(a r - lam s - mu t)"
        " = 13^4*T_DEN*B(a,-lam,-mu), THM-2365 (10)")
    log("    per-chamber contribution c_i = len_i * phi_i(a) * V_i(lam,mu),"
        "  phi_i(a) = sum_r dc_i[r] zeta^(a r)")
    log("    (r,s,t all marginalised on the COMMON chamber base before any"
        " aggregation: MISTAKE-281)")
    log("    target of endpoint(a,lam,mu) is (-lam, a-mu) (THM-2365 (13));"
        " zero-target pair on line a is (lam,mu)=(0,a)")
    log("    AMBIGUITY (logged): task count '2,184 = 12 x 182' = TOTAL"
        " nonzero-target modes over all a in F_13")
    log("    (13 x 168 = 2184, including the 168 a=0 modes already decided);"
        " new a != 0 nonzero-target modes: 12 x 168 = 2016;")
    log("    zero-target modes (13, one per line) are flagged separately and"
        " never counted as drift mechanisms.")
    log("    all 13 x 169 = 2197 modes are computed; a = 0 reproduces the"
        " sidecar (gate F)")
    log("    floats route/describe only; every decision (endpoint zero test,"
        " cone class, event index) is exact")
    log("")

    # ------------------------------------------------------------------
    log("[2] chamber partition (sidecar machinery, re-run)")
    allpts, slowset, counts = sc.gen_breakpoints()
    for k, v in counts.items():
        log(f"    {k}: {v} generated breakpoints")
    n_chambers = int(allpts.size)
    require(n_chambers == 1489966, "chamber count anchor (sidecar [2])")
    log(f"    distinct breakpoints = chambers: {n_chambers} (anchor PASS)")
    from bisect import bisect_left, bisect_right
    LJ = 13 * T_DEN // 14
    RJ = 15 * T_DEN // 14
    pts_list = allpts.tolist()
    lo_i = bisect_right(pts_list, LJ)
    hi_wrap = T_DEN // 14
    inner = [p for p in pts_list[lo_i:] if p < T_DEN] \
        + [p + T_DEN for p in pts_list[:bisect_left(pts_list, hi_wrap)]]
    sweep = [LJ] + inner + [RJ]
    for i in range(len(sweep) - 1):
        require(sweep[i] < sweep[i + 1], "sweep not strictly increasing")
    n_sup = len(sweep) - 1
    require(n_sup == 212845, "contributing-chamber anchor (sidecar [2])")
    log(f"    contributing chambers in supp d(y) = (13/14, 15/14):"
        f" {n_sup} (anchor PASS)")
    log(f"    {el()} partition done")
    log("")

    # ------------------------------------------------------------------
    log("[3] exact 1-D sweep (sidecar machinery, re-run; states indexed)")
    states = {}
    key_index = {}
    chamber_seq = []            # (x0, state_idx) in y-order
    rooted_ids = {}
    rooted_list = []
    cur_rid = -1
    cur_ga = -1
    slow_in_J = set()
    for p in slowset:
        slow_in_J.add(p)
        slow_in_J.add(p + T_DEN)
    first = True
    for i in range(n_sup):
        x0 = sweep[i]
        x1 = sweep[i + 1]
        m = x0 + x1
        if first or x0 in slow_in_J:
            require(sc.owner_d(m), "chamber outside owner support in sweep")
            snap = sc.rooted_snapshot(m)
            rid = rooted_ids.get(snap)
            if rid is None:
                rid = len(rooted_list)
                rooted_ids[snap] = rid
                rooted_list.append(snap)
            cur_rid = rid
            cur_ga = sc.ga_bits(m)
            first = False
        key = (cur_rid, cur_ga, sc.dc_bits(m))
        ln = x1 - x0
        if key in states:
            states[key] += ln
        else:
            key_index[key] = len(states)
            states[key] = ln
        chamber_seq.append((x0, key_index[key]))
    require(sum(states.values()) == T_DEN // 7, "support length mismatch")
    nst = len(states)
    require(nst == 10998 and len(rooted_list) == 19,
            "state/snapshot anchors (sidecar [3])")
    log(f"    chambers walked: {n_sup}; rooted snapshots:"
        f" {len(rooted_list)}; states (rooted,ga,dc): {nst} (anchors PASS)")
    log(f"    {el()} sweep done")
    log("")

    # ------------------------------------------------------------------
    log("[4] rooted N tables, state arrays, per-character V tables")
    N_tabs = []
    for (a_, c_, p_) in rooted_list:
        av = [(a_ >> r) & 1 for r in range(13)]
        cv = [(c_ >> r) & 1 for r in range(13)]
        pv = [(p_ >> h) & 1 for h in range(13)]
        require(0 < sum(av) < 13 and 0 < sum(cv) < 13, "A/C not proper")
        require(3 <= sum(pv) <= 10, "P support outside THM-2368 (14)")
        N = np.zeros((13, 13), dtype=np.int64)
        for s in range(13):
            for t in range(13):
                N[s, t] = sum(av[(h + s) % 13] * cv[(h + t) % 13] * pv[h]
                              for h in range(13))
        N_tabs.append(N)
    log(f"    rooted words proper, |supp P| in [3,10]: PASS"
        f" ({len(rooted_list)} snapshots)")

    LN = np.zeros(nst, dtype=np.int64)
    RID = np.zeros(nst, dtype=np.int64)
    GAV = np.zeros((nst, 13), dtype=np.int64)
    DCV = np.zeros((nst, 13), dtype=np.int64)
    for key, ln in states.items():
        idx = key_index[key]
        rid, ga, dc = key
        LN[idx] = ln
        RID[idx] = rid
        for b in range(13):
            GAV[idx, b] = (ga >> b) & 1
            DCV[idx, b] = (dc >> b) & 1
    NT = np.stack(N_tabs)
    Y = GAV[:, :, None] * (1 - DCV)[:, None, :] * NT[RID]     # (nst,13,13)
    Yflat = Y.reshape(nst, 169)
    require(int(Yflat.max()) <= 1690, "Y bound (state table)")

    Hnum = np.einsum("s,sr,sij->rij", LN, DCV, Y)             # 13*T_DEN scale

    CH169 = [(lam, mu) for lam in range(13) for mu in range(13)]
    SCAT2 = np.zeros((2197, 169), dtype=np.int64)
    for ci, (lam, mu) in enumerate(CH169):
        for j in range(169):
            s, t = divmod(j, 13)
            k = (-(lam * s + mu * t)) % 13
            SCAT2[ci * 13 + k, j] = 1
    Vall = (Yflat @ SCAT2.T).reshape(nst, 169, 13)            # (nst,169,13)
    require(int(np.abs(Vall).max()) <= 1690, "V bound")
    log(f"    state arrays built; Vall (169 characters x 13 zeta-classes per"
        f" state) assembled")
    log(f"    {el()} state processing done")
    log("")

    # ------------------------------------------------------------------
    log("[5] GATES (A)+(B)+(C) (sidecar gates, re-run verbatim)")
    E0 = sc.build_set(sc.PAT_E, sc.ZELL)
    lenE0 = sc.check_intervals(E0, T_DEN)
    require(Fraction(lenE0, T_DEN) == Fraction(1882176, 28589561),
            "E_1 measure anchor")
    log("    E_1 measure anchor 1882176/28589561: PASS")
    HE_tab = {}
    SE_comp = {}
    for s in range(13):
        for t in range(13):
            iv = sc.build_set(sc.PAT_E, sc.twist_shift_vector(s, t))
            lenE = 0
            HE = [0] * 13
            succ = 0
            for A, B in iv:
                lenE += B - A
                succ += sc.accumulate_interval(
                    A, B, sc.PP_T, sc.HPROBE_T, sc.HSUCC_T, HE)
            require(sum(HE) == 2 * lenE - succ,
                    f"(19b) identity fails at {(s, t)}")
            for r in range(13):
                HE_tab[(r, s, t)] = HE[r]
            SE_comp[(s, t)] = sum(HE)
    log("    (19b) successor identity at all 169 twists: PASS")
    gate_a = all(int(Hnum[r, s, t]) == 13 * HE_tab[(r, s, t)]
                 for r in range(13) for s in range(13) for t in range(13))
    log(f"    (A) chamber reconstruction == independent 13^3 table (x13"
        f" scale), all 2197 cells: {'PASS' if gate_a else 'FAIL'}")
    require(gate_a, "GATE (A) failed")
    with open(sc.DRIFT_OUT, "r", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    idx0 = lines.index("    S_E(s,t)*T_DEN:")
    SE_stored = {}
    for s in range(13):
        row = [int(x) for x in lines[idx0 + 1 + s].split()]
        require(len(row) == 13, "stored S_E row malformed")
        for t in range(13):
            SE_stored[(s, t)] = row[t]
    gate_b = all(SE_comp[(s, t)] == SE_stored[(s, t)]
                 and int(Hnum[:, s, t].sum()) == 13 * SE_stored[(s, t)]
                 for s in range(13) for t in range(13))
    log(f"    (B) row sums == byte-stored S_E numerators:"
        f" {'PASS' if gate_b else 'FAIL'}")
    require(gate_b, "GATE (B) failed")
    diag_ok = all(int(Hnum[t, s, t]) == 0
                  for s in range(13) for t in range(13))
    log(f"    (C) diagonal-plane zero H_E(t,s,t)=0 (all 169):"
        f" {'PASS' if diag_ok else 'FAIL'}")
    require(diag_ok, "diagonal zero failed")
    log(f"    {el()} gates done")
    log("")

    # ------------------------------------------------------------------
    log("[6] per-a-line walks, endpoints, exact cone decisions, events")
    log("    aggregation (logged, identical to sidecar stage 2):"
        " contributions merged by primitive direction")
    log("    in Z[zeta_13] (positive-rational-multiple classes)")
    IDXK = np.zeros((13, 13), dtype=np.int64)
    for m in range(13):
        for k in range(13):
            IDXK[m, k] = (k - m) % 13
    R3, S3, T3 = np.indices((13, 13, 13))
    Hrav = Hnum.ravel()

    def alleq(vec):
        return all(int(x) == int(vec[0]) for x in vec)

    END = []                     # END[a][ci] = list of 13 python ints
    line_stats = []
    cone_holders = []            # (a, lam, mu, cls, info, dirs list)
    global_min_excess = None     # (excess, a, lam, mu)  DESCRIPTIVE
    events = []                  # per line dict

    for a in range(13):
        ta = time.time()
        if a == 0:
            PHI = np.zeros((nst, 13), dtype=np.int64)
            PHI[:, 0] = DCV.sum(axis=1)
        else:
            PHI = np.zeros((nst, 13), dtype=np.int64)
            for r in range(13):
                PHI[:, (a * r) % 13] = DCV[:, r]
        Bm = PHI[:, IDXK]                       # (nst,13,13)
        R = np.matmul(Vall, Bm)                 # (nst,169,13) exact int64
        require(int(np.abs(R).max()) <= 21970, "R bound (overflow guard)")
        E = np.einsum("s,sck->ck", LN, R)       # (169,13) exact int64
        require(int(np.abs(E).max()) <= 21970 * (T_DEN // 7),
                "endpoint bound (overflow guard)")
        # GATE (D): direct exact DFT binning of the gated tensor
        for ci, (lam, mu) in enumerate(CH169):
            idxm = ((a * R3 - lam * S3 - mu * T3) % 13).ravel()
            acc = np.zeros(13, dtype=np.int64)
            np.add.at(acc, idxm, Hrav)
            require((acc == E[ci]).all(),
                    f"GATE (D) endpoint mismatch at (a,lam,mu)="
                    f"({a},{lam},{mu})")
        # primitive directions
        Rred = R - R[:, :, 12:13]
        G = np.gcd.reduce(np.abs(Rred), axis=2)
        maskD = G > 0
        prims = Rred // np.where(maskD, G, 1)[:, :, None]

        zt_ci = a                               # zero-target pair (0,a)
        nz = [not alleq(E[ci]) for ci in range(169)]
        cls_counts = {"<pi": 0, "=pi-cert": 0, "=pi-fail": 0, ">pi": 0,
                      "empty": 0}
        ndirs_list = []
        spans = []
        route_counts = {"cert3": 0, "full": 0, "trivial": 0}
        zt_cls = None
        evt_ci = None
        for ci, (lam, mu) in enumerate(CH169):
            dirs_np = np.unique(prims[maskD[:, ci], ci, :], axis=0)
            cls, info, ndirs, span, route = decide_mode(dirs_np)
            if ci == zt_ci:
                zt_cls = (cls, ndirs, nz[ci])
                continue                        # zero-target: flagged only
            cls_counts[cls] += 1
            ndirs_list.append(ndirs)
            route_counts[route] += 1
            if span is not None:
                spans.append((span, lam, mu))
            if cls in FEASIBLE:
                dl = [tuple(int(v) for v in row) for row in dirs_np]
                cone_holders.append((a, lam, mu, cls, info, dl))
            if evt_ci is None and nz[ci] and cls in INFEASIBLE:
                evt_ci = ci
                evt_mask = maskD[:, ci].copy()
                evt_prims = prims[:, ci, :].copy()
        END.append([[int(x) for x in E[ci]] for ci in range(169)])
        nz_count = sum(1 for ci in range(169) if ci != zt_ci and nz[ci])
        if spans:
            smin = min(spans)
            smax = max(spans)
            if a != 0 and (global_min_excess is None
                           or smin[0] - pi < global_min_excess[0]):
                global_min_excess = (smin[0] - pi, a, smin[1], smin[2])
        # first cancelling event for the representative mode
        evt = None
        if evt_ci is not None:
            evt = event_scan(chamber_seq, evt_mask, evt_prims)
        events.append((a, evt_ci, evt))
        line_stats.append((a, nz_count, dict(cls_counts), min(ndirs_list),
                           max(ndirs_list),
                           sum(ndirs_list) / len(ndirs_list),
                           (smin if spans else None),
                           (smax if spans else None), zt_cls,
                           dict(route_counts)))
        lam_e, mu_e = (CH169[evt_ci] if evt_ci is not None else (None, None))
        log(f"    a={a:2d}: nonzero-target modes=168, endpoints nonzero="
            f"{nz_count}/168, classes {cls_counts},")
        log(f"          dirs min={min(ndirs_list)} max={max(ndirs_list)}"
            f" mean={sum(ndirs_list)/len(ndirs_list):.1f}; routes"
            f" {route_counts}")
        if spans:
            log(f"          DESCRIPTIVE float spans: min={smin[0]:.9f} rad"
                f" at (lam,mu)=({smin[1]},{smin[2]}),"
                f" max={smax[0]:.9f} rad; min excess over pi ="
                f" {smin[0]-pi:.9f}")
        log(f"          zero-target mode (lam,mu)=(0,{a}): class={zt_cls[0]}"
            f" dirs={zt_cls[1]} endpoint_nonzero={zt_cls[2]}"
            f" (NOT a drift mode; flagged only)")
        if evt is not None:
            kk, x0e, prime, ecls = evt
            yfrac = Fraction(x0e % T_DEN, T_DEN)
            log(f"          first cancelling event, mode (a,lam,mu)="
                f"({a},{lam_e},{mu_e}): direction #{kk},"
                f" chamber y={yfrac}")
            log(f"          (grid {x0e % T_DEN}/T_DEN), class at event ="
                f" {ecls}, direction {tuple(prime)}")
        else:
            log("          no cancelling event (representative mode stayed"
                " cone-feasible or no candidate)")
        log(f"          {el()} line a={a} done"
            f" ({time.time()-ta:6.1f}s)")
        del R, Rred, prims, G, maskD, Bm, PHI
    log("")

    # ------------------------------------------------------------------
    log("[7] GATE (E): THM-2365 (14) target-null circulation, all 169"
        " target lines")
    for b in range(13):
        for q2 in range(13):
            Ssum = [0] * 13
            lam = (-b) % 13
            for a in range(13):
                mu = (a - q2) % 13
                ci = lam * 13 + mu
                for k in range(13):
                    Ssum[k] += END[a][ci][k]
            require(alleq(Ssum),
                    f"circulation fails at target (b,q2)=({b},{q2})")
    log("    sum_a endpoint(a,-b,(a-q2)) = 0 in Z[zeta_13] for ALL 169"
        " targets (b,q2): PASS")
    log("")

    # ------------------------------------------------------------------
    log("[8] GATE (F): a = 0 line reproduces the stored sidecar artifact")
    a0 = line_stats[0]
    _, nz0, cc0, dmin0, dmax0, dmean0, _, _, _, _ = a0
    require(nz0 == 168, "a=0 endpoint-nonzero count != 168")
    require(cc0[">pi"] == 168 and cc0["<pi"] == 0 and cc0["=pi-cert"] == 0
            and cc0["=pi-fail"] == 0 and cc0["empty"] == 0,
            "a=0 class counts differ from sidecar")
    require(dmin0 == 5740 and dmax0 == 10322
            and f"{dmean0:.1f}" == "9102.9",
            "a=0 direction statistics differ from sidecar")
    evt_ci0 = events[0][1]
    require(evt_ci0 is not None and CH169[evt_ci0] == (0, 1),
            "a=0 representative event mode != (0,1)")
    kk0, x0e0, prim0, _ = events[0][2]
    require(x0e0 % T_DEN == 276563549922660,
            "a=0 event chamber differs from sidecar anchor")
    require(tuple(prim0) == (-6, 0, 0, 0, -1, -2, -2, -2, -2, -2, -1, 0, 0),
            "a=0 event direction differs from sidecar anchor")
    log("    endpoints 168/168 nonzero, classes 168x '>pi', dirs"
        " min=5740 max=10322 mean=9102.9: PASS")
    log("    first event (0,1) at grid 276563549922660/T_DEN with direction"
        " (-6,0,0,0,-1,-2,-2,-2,-2,-2,-1,0,0): PASS")
    log("")

    # ------------------------------------------------------------------
    log("[9] DECISION (the 2016 deep-coloured a != 0 nonzero-target modes)")
    deep_holders = [h for h in cone_holders if h[0] != 0]
    n_nz_deep = sum(st[1] for st in line_stats[1:])
    log(f"    endpoints nonzero: {n_nz_deep} / 2016 deep-coloured"
        f" nonzero-target modes")
    if deep_holders:
        log("    PHASE-CONE CERTIFICATE HOLDS on at least one deep-coloured"
            " a != 0 mode line.")
        log(f"    certified modes ({len(deep_holders)}):"
            f" {[(h[0], h[1], h[2]) for h in deep_holders]}")
        shown = 0
        for (a, lam, mu, cls, info, dl) in deep_holders:
            unit = None
            for k in range(13):
                for sg in (1, -1):
                    sig = [sc.re_sign([sg * x for x in sc.rot13(list(d), k)])
                           for d in dl]
                    if all(s >= 0 for s in sig) and any(s > 0 for s in sig):
                        unit = (sg, k)
                        break
                if unit:
                    break
            start, end = info
            if unit:
                sgs = "" if unit[0] == 1 else "-"
                log(f"    mode (a,lam,mu)=({a},{lam},{mu}): certifying unit"
                    f" u = {sgs}zeta^{unit[1]} (exact 26th-root sector);")
            else:
                log(f"    mode (a,lam,mu)=({a},{lam},{mu}): no 26th-root"
                    f" unit; exact sector = closed arc from")
            log(f"        extreme directions {start} .. {end}"
                f" (arc {'< pi' if cls == '<pi' else '= pi'})")
            shown += 1
            if shown >= 24:
                log(f"    ... ({len(deep_holders)-shown} further certified"
                    f" modes suppressed)")
                break
        log("")
        log("    CONSEQUENCE: this is the first cone-explainable positivity"
            " mechanism on the row: the named")
        log("    mode(s) survive inside one closed half-plane after one unit"
            " rotation -- the anchor the uniform")
        log("    THM-2368 (38) argument needs, now at deep colour a != 0.")
    else:
        log("    PHASE-CONE CERTIFICATE FAILS for ALL 2016 deep-coloured"
            " a != 0 nonzero-target modes.")
        log("    Combined with the sidecar's a = 0 negative (168/168 fail,"
            " gate F) the cone mechanism is DEAD at")
        log("    EVERY mode depth on this row: all 2184 nonzero-target modes"
            " of THM-2365 (10) have phase spans")
        log("    exceeding (or exactly filling) a half turn; the certified"
            " positive drift D_(H_E) > 0 (THM-2550)")
        log("    survives only through incomplete cancellation of"
            " full-rotation phase words at every depth.")
        if global_min_excess is not None:
            ex, ax, lx, mx = global_min_excess
            log(f"    sharpest deep mode (DESCRIPTIVE float): (a,lam,mu)="
                f"({ax},{lx},{mx}) with span excess {ex:.9f} rad over pi")
    log("")

    # ------------------------------------------------------------------
    log("[10] scope (MISTAKE-281/283 guardrails, unchanged)")
    log("    * The row is the canonical TYPED row THM-2309 (25): typed, NOT"
        " an asserted scalar cover.")
    log("    * No physical current is claimed; no scalar row is excluded;"
        " the scalar ledger stays 165;")
    log("      LRC(14) and JC(2) remain OPEN.  Full support is never read as"
        " noncancellation.")
    log("    * A failed cone at every depth excludes ONE mechanism class"
        " (half-plane/cone certificates on")
    log("      THM-2365 (10) modes of the bare tensor); it does NOT bound"
        " the drift, exclude other mechanisms")
    log("      (symmetry quotients, pre-y-integration breaks), or license"
        " any row-uniform statement.")
    log("")
    log(f"{el()} all checks passed")


# ---------------------------------------------------------------------------
# First-cancelling-event scan (float-accelerated, EXACTLY certified)
# ---------------------------------------------------------------------------
def event_scan(chamber_seq, mask_col, prims_col):
    """Scan distinct primitive directions in chamber (y) order; return the
    minimal prefix whose direction set is exactly cone-infeasible, as
    (k, chamber_x0, direction, exact_class), or None if the full set is
    feasible.  Floats only bound the search; the returned index is fixed by
    exact bisection (feasibility is monotone under adding directions)."""
    seen = []
    seen_set = set()
    addx = []
    angs = []
    flagged = None
    for (x0, sidx) in chamber_seq:
        if not mask_col[sidx]:
            continue
        tpl = tuple(int(v) for v in prims_col[sidx])
        if tpl in seen_set:
            continue
        seen.append(tpl)
        seen_set.add(tpl)
        addx.append(x0)
        xx = sum(tpl[j] * EMBC[j] for j in range(13))
        yy = sum(tpl[j] * EMBS[j] for j in range(13))
        insort(angs, atan2(yy, xx))
        n = len(seen)
        if n < 2:
            continue
        if n <= 2000 or n % 25 == 0:
            mg = max(angs[i + 1] - angs[i] for i in range(len(angs) - 1))
            mg = max(mg, angs[0] + 2 * pi - angs[-1])
            if mg <= pi + BAND:
                cls, _ = cone_test(seen)
                if cls in INFEASIBLE:
                    flagged = n
                    break
    if flagged is None:
        cls, _ = cone_test(seen)
        if cls not in INFEASIBLE:
            return None
        flagged = len(seen)
    lo, hi = 2, flagged
    while lo < hi:
        mid = (lo + hi) // 2
        cls, _ = cone_test(seen[:mid])
        if cls in INFEASIBLE:
            hi = mid
        else:
            lo = mid + 1
    cls, _ = cone_test(seen[:lo])
    require(cls in INFEASIBLE, "event bisection lost infeasibility")
    if lo > 2:
        cls_prev, _ = cone_test(seen[:lo - 1])
        require(cls_prev in FEASIBLE, "event bisection not minimal")
    return (lo, addx[lo - 1], seen[lo - 1], cls)


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
