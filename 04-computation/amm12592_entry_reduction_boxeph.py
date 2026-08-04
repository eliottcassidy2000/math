"""LANE F1 -- THE ENTRY REDUCTION: waiting-room theorem E-1 + parameter
certification E-2 + feed-end hypothesis certification (ENTRY-eta).

Session: boxeph multifront, 2026-08-04.  All decisions exact int/Fraction;
floats only in display fields; no numpy; no sympy.

Mathematical content (proofs in
05-knowledge/results/amm12592-entry-proof-boxeph.md):

  W1 (cell-0, unconditional): j_0 <= 0 always; post-feed a_0 drains exactly
     2/row; F2 (a_0 <= d-1) holds at every row >= i_pf + k2 with
     k2 = min{k : (R-2-2k)+ <= d_fe + floor(g_lo k) - 1}  (absorbing).
  W2 (sign decoupling): junk = p - n; one post-feed row obeys, cellwise,
     p' <= max(0, (K p) - capU),  n' <= max(0, (K n) - capL),
     capU_t = 2C(d'-1, t-1), capL_t = 2C(d'-1, t).  (Each side sees only
     its own majorant flow.)
  W3 (layer gain): during the wait, cell 1 gains at most
     G_ub = sum_l 2*(R-1-2l - d_fe_lb - floor(g_lo l))+   (per-row bound,
     valid for every delta word; no drain credit taken).
  W4 (over-cell burn): a core cell whose spill excess over capref is
     eps_t * capref_t stops growing after k_stop(t) rows (cap drift), with
     total growth B_t <= k_stop(t) * eps_t * capref_t; the excesses obey
     the upward recursion eps_{t+1} = eta + (2B_t + B_{t-1})/capref_{t+1},
     explicit for t <= 16, then a scalar fixed-point invariant E.
  E-1 (waiting-room theorem): H-P /\ H-S(eta, theta) at feed-end + Psi(R)
     ==> F1/\F2/\F3/\F4 at i0 = i_pf + k*  ==>  S(R)  (via S-cone-fc,
     two-sided form).
  E-2: Psi(R) certified for dyadic R, eps in {1/32, 1/16}  (this script,
     mode `psi`).

Modes:
  psi        exact-rational sweep of the Psi side conditions, R = 2^9..2^40
  snapcert   H-P / H-S margins on every stored feed-end snapshot
  waitsim    simulate the wait window from each snapshot; verify F1-F3 at
             the E-1 row i0 (consistency of the theorem with reality)
  feedonly R D0   run the feed phase only; write the standard snapshot

Outputs -> 05-knowledge/results/amm12592_entry_{psi,snapcert,waitsim}_boxeph.*
           05-knowledge/results/amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb, isqrt
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

# certified sandwich (fresh re-derivation in amm12592_S_cone_constants_boxeph)
G_LO = Fr(627035, 1 << 20)
G_HI = Fr(627036, 1 << 20)

# E-1 parameters (frozen for the theorem statement)
ETA = Fr(1, 64)          # core tolerance:   r_t <= 1 + eta   on [2, t_c]
THETA = Fr(1, 4)         # zone-2 margin:    r_t <= 1 - theta on (t_c, m+2]
EPS_CAP = Fr(1, 8)       # required uniform bound on all excesses eps_t
T_EXPL = 16              # explicit eps-recursion range; scalar invariant after

CONFIGS = [(128,4),(128,8),(256,8),(256,16),(512,16),(512,32),(1024,32),
           (1024,64),(2048,64),(2048,128),(4096,128),(4096,256),(8192,256),
           (8192,512),(16384,512),(16384,1024),(32768,1024),(32768,2048)]


def tc_of(m):
    return (3 * m) // 4 if m >= 12 else 2


def binom_frac(x, t):
    """C(x, t) for Fraction/int x, integer t >= 0 (falling factorial)."""
    out = Fr(1)
    for s in range(t):
        out *= (x - s)
        out /= (s + 1)
    return out


def psi_one(R, eps_num, eps_den, verbose=False):
    """All Psi side conditions for one (R, eps), exact rationals.
    Returns dict with pass/fail per condition + the derived constants."""
    D0 = -((-R * eps_num) // eps_den)
    rec = {"R": R, "D0": D0, "eps": f"{eps_num}/{eps_den}"}
    # --- profile bounds (Theorem B / S-cone conventions)
    d0_lb = G_LO * R - 1 + D0          # d_0 = floor(gR) + D0
    d0_ub = G_HI * R + D0
    ifeed_lb = (R * (1 - G_HI) - D0) / (1 + G_HI) - 1
    ifeed_ub = (R * (1 - G_LO) - D0) / (1 + G_LO)
    ipf_lb = ifeed_lb + 1
    ipf_ub = ifeed_ub + 2
    dfe_lb = G_LO * (R + ipf_lb - 1) - 1 + D0
    dfe_ub = G_HI * (R + ipf_ub - 1) + D0
    rec["window_ii"] = (Fr(R, 2) < d0_lb) and (d0_ub < Fr(2 * R, 3))
    rec["ipf_le_half"] = ipf_ub <= Fr(R - 2, 2)          # S3 (drain deadline)
    # --- W1: F2 wait k2 (worst case n_0 = R-2, degree staircase floor(g_lo k))
    def f2_ok(k):
        return max(0, R - 2 - 2 * k) <= dfe_lb + (G_LO * k).__floor__() - 1
    lo, hi = 0, R
    while lo < hi:
        mid = (lo + hi) // 2
        if f2_ok(mid):
            hi = mid
        else:
            lo = mid + 1
    k2 = lo
    rec["k2"] = k2
    # --- W3: layer gain G_ub = sum_l 2*(R-1-2l - dfe_lb - floor(g_lo l))+
    # closed staircase sum; term positive while (2+g_lo) l < W', W' = R-1-dfe_lb
    Wp = R - 1 - dfe_lb
    Lmax = int(Wp / (2 + G_LO)) + 2
    G_ub = Fr(0)
    for l in range(max(0, Lmax)):
        term = Wp - 2 * l - (G_LO * l).__floor__()
        if term > 0:
            G_ub += 2 * term
    rec["G_ub_over_R2"] = float(G_ub / R ** 2)
    # delta=0 layer rows gain nothing: n_0 <= R-2 < 2(dfe_lb - 1) ?
    rec["delta0_layer_free"] = (R - 2) < 2 * (dfe_lb - 1)
    # --- W4: eps recursion, explicit for t = 2..T_EXPL, capref at dfe bounds.
    # capref_t within [capL(t), capU(t)] from the dfe sandwich:
    def capref_lb(t):
        return 2 * binom_frac(dfe_lb - 1, t)
    def capref_ub(t):
        return 2 * binom_frac(dfe_ub - 1, t)
    def kstop(eps_t, t):
        # smallest k with floor(g_lo k) >= eps_t * dfe_ub / t   (cap-growth
        # lower bound C(D+s-1,t)/C(D-1,t) >= 1 + t*s/D_ub, s = floor(g_lo k))
        need = eps_t * dfe_ub / t
        k = int(((need + 1) / G_LO)) + 2
        assert (G_LO * k).__floor__() >= need
        return k
    eps = {1: None}                      # eps_t; cell 1 handled via G_ub
    B = {0: Fr(0), 1: G_ub}              # total-growth bounds B_t
    ks = {}
    eps[2] = ETA + (2 * B[1] + 2 * (R - 2)) / capref_lb(2)   # n_0 term <= R-2
    ks[2] = kstop(eps[2], 2)
    B[2] = ks[2] * eps[2] * capref_ub(2)
    ok_eps = eps[2] <= EPS_CAP
    for t in range(3, T_EXPL + 1):
        eps[t] = ETA + (2 * B[t - 1] + B[t - 2]) / capref_lb(t)
        ks[t] = kstop(eps[t], t)
        B[t] = ks[t] * eps[t] * capref_ub(t)
        ok_eps = ok_eps and eps[t] <= EPS_CAP
    rec["eps2"] = float(eps[2]); rec["eps3"] = float(eps[3])
    rec["epsT"] = float(eps[T_EXPL])
    rec["kstop2"] = ks[2]
    # scalar tail invariant: for t in (T_EXPL, t_c_max]:  eps <= E with
    # E = eta + cmax * E^2 (+ two-step term c2max * E^2), where
    # c(t) = 2 k_stop(t) capref(t)/capref(t+1) growth folded:
    #   bump_{t+1} <= 2 B_t / capref_{t+1} <= 2 kstop(t) eps_t capU(t)/capL(t+1)
    # with kstop(t) <= (eps_t dfe_ub/t + 1)/g_lo + 2:
    #   c(t) ~ 2 dfe_ub (t+1) / (g_lo t (dfe_lb - 1 - t)) * capU/capL slop.
    # t_c_max <= 3 m_ub/4; m_ub from T6a:
    m_ub_frac = R - 1 - 2 * d0_lb + dfe_ub + ipf_ub
    m_ub = int(m_ub_frac) + 1
    rec["m_ub_over_R"] = float(Fr(m_ub, R))
    tc_max = (3 * m_ub) // 4
    slop = capref_ub(T_EXPL) / capref_lb(T_EXPL)   # dfe-sandwich slop, >= 1
    def ccoef(t):
        return 2 * dfe_ub * (t + 1) / (G_LO * t * (dfe_lb - 1 - t)) * slop \
            + 2 * Fr(2 + t, 1) / (G_LO * t)  # +2 rows slack in kstop, coarse
    cmax = max(ccoef(T_EXPL), ccoef(tc_max))       # c is convex in t; ends
    # two-step term B_{t-1}/capref_{t+1}: coefficient <= cmax * (t+1)/(dfe-t)
    c2max = cmax * Fr(tc_max + 1, 1) / (dfe_lb - 1 - tc_max)
    # find rational E in (eta, EPS_CAP] with  eta + (cmax + c2max) E^2 <= E
    ctot = cmax + c2max
    E = None
    cand = 2 * ETA
    for _ in range(40):
        if ETA + ctot * cand * cand <= cand:
            E = cand
            break
        cand *= Fr(11, 10)
        if cand > EPS_CAP:
            break
    rec["cmax"] = float(cmax); rec["E_tail"] = float(E) if E else None
    ok_tail = E is not None and eps[T_EXPL] <= E and E <= EPS_CAP
    # kstop decreasing beyond T_EXPL: kstop(t) <= (E dfe_ub/t)/g_lo + 3
    kstop_tail_ub = int((E * dfe_ub / (T_EXPL + 1)) / G_LO) + 3 if E else None
    kstar = max(k2, max(ks.values()), kstop_tail_ub or 0)
    rec["kstar"] = kstar
    rec["ok_eps_expl"] = ok_eps
    rec["ok_tail"] = ok_tail
    # --- support / no-death: T6a bound at i_pf, frozen during wait
    rec["support_ok"] = m_ub + 2 * 0 + 4 < dfe_lb   # frozen: no +2k term
    # theta budget at the zone boundary (t_c+1, t_c+2):
    # bump <= (2 B_tc + B_{tc-1})/capref_{tc+1} <= (cmax+c2max) E^2 <= theta?
    rec["theta_budget_ok"] = E is not None and ctot * E * E <= THETA
    # --- F4' budget: shifted-window a-priori clocks (a_t(i0) <= C(D-1,t+1))
    i0_ub = ipf_ub + kstar
    Dub = G_HI * (R + i0_ub) + D0 + 1          # degree upper bound at i0
    Dlb = dfe_lb                               # degree lower bound at i0
    need1 = Dub * (Dub - 1) / 2                # a_1 <= C(Dub-1, 2)
    need2 = Dub * (Dub - 1) * (Dub - 2) / 6    # a_2 <= C(Dub-1, 3)
    def k1ok(K):        # cell-1 drain-assisted clock (a-priori corollary)
        return (G_LO + 2) * K * K - 4 * K >= need1
    def k2ok(K):        # t=2 cap-drift clock
        return (G_LO * K * K / 2 - K) * (2 * Dlb - 3) >= need2
    def bs(ok):
        lo, hi = 1, 8 * R
        if not ok(hi):
            return None
        while lo < hi:
            mid = (lo + hi) // 2
            if ok(mid):
                hi = mid
            else:
                lo = mid + 1
        return lo
    K1c, K2c = bs(k1ok), bs(k2ok)
    # p-part clocks: T_p(t) <= (Dub - t)/(2t) max at t = 2; p_1 drains 2/row
    Tp = (Dub - 2) / 4
    Tp1 = Dub / 4    # p_1 <= capU_2/4 = (Dub-1)/2; /2 per row
    rec["K1c"] = K1c; rec["K2c"] = K2c
    budget = i0_ub + max(K1c or 10**9, K2c or 10**9, Tp, Tp1)
    rec["budget_used_over_R"] = float(budget / R)
    rec["F4_ok"] = (K1c is not None and K2c is not None and
                    budget <= R - 2)
    # drain clock is free:  i0 + n_0(i0)/2 <= i_pf + (R-2)/2 <= R-2  (S3)
    rec["drain_free"] = rec["ipf_le_half"]
    rec["PSI_PASS"] = all(rec[k] for k in
                          ("window_ii", "ipf_le_half", "delta0_layer_free",
                           "ok_eps_expl", "ok_tail", "support_ok",
                           "theta_budget_ok", "F4_ok", "drain_free"))
    if verbose:
        print(f"  R=2^{R.bit_length()-1} D0={D0}: PASS={rec['PSI_PASS']} "
              f"k2={k2} kstar={kstar} eps2={rec['eps2']:.4f} "
              f"E={rec['E_tail']} m_ub/R={rec['m_ub_over_R']:.3f} "
              f"budget/R={rec['budget_used_over_R']:.3f}", flush=True)
    return rec


def psi_sweep():
    out = {"params": {"eta": [1, 64], "theta": [1, 4], "eps_cap": [1, 8],
                      "T_expl": T_EXPL, "tc": "floor(3m/4) if m>=12 else 2"},
           "sweeps": {}}
    for en, ed in ((1, 32), (1, 16)):
        print(f"[psi] eps = {en}/{ed}:", flush=True)
        recs = []
        first_pass = None
        for kk in range(9, 41):
            r = psi_one(1 << kk, en, ed, verbose=(kk <= 20 or kk % 5 == 0))
            recs.append(r)
            if r["PSI_PASS"] and first_pass is None:
                first_pass = 1 << kk
            if not r["PSI_PASS"]:
                first_pass = None
        allpass_from = None
        for r in recs:
            if r["PSI_PASS"] and allpass_from is None:
                allpass_from = r["R"]
            if not r["PSI_PASS"]:
                allpass_from = None
        out["sweeps"][f"{en}_{ed}"] = {"records": recs,
                                       "all_pass_from": allpass_from}
        print(f"[psi] eps={en}/{ed}: PSI holds for all dyadic R >= "
              f"{allpass_from}", flush=True)
    json.dump(out, open(os.path.join(
        RESULTS, "amm12592_entry_psi_boxeph.json"), "w"), indent=1,
        default=str)
    return out


# ---------------------------------------------------------------- snapshots
def load_snap(R, D0):
    p = os.path.join(RESULTS,
                     f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json")
    S = json.load(open(p))
    j = {int(t): int(v) for t, v in S["junk"].items()}
    return S["i_postfeed"], S["d_feedend"], j


def snapcert():
    """ENTRY-eta certificate on every stored feed-end snapshot: H-P (p = 0,
    or p under the shifted half-cap cone) and H-S (two-zone surface)."""
    rows = []
    for R, D0 in CONFIGS + [(65536, 2048), (65536, 4096), (131072, 4096)]:
        p = os.path.join(RESULTS,
                         f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json")
        if not os.path.exists(p):
            continue
        ipf, dfe, j = load_snap(R, D0)
        pos = {t: v for t, v in j.items() if v > 0}
        n = {t: -v for t, v in j.items() if v < 0}
        m = max(j)
        tc = tc_of(m)
        capref = [2 * comb(dfe - 1, t) for t in range(m + 3)]
        # H-P: positive part under (1/2) * shifted caps  (capU_t = 2C(d-1,t-1))
        hp_ok = (0 not in pos) and (1 not in pos)
        for t in range(1, m + 3):
            spl = 2 * pos.get(t - 1, 0) + pos.get(t - 2, 0)
            if 2 * spl > 2 * comb(dfe - 1, t - 1):
                hp_ok = False
        # H-S: two zones
        eta_need = Fr(0)
        z2_max = Fr(0)
        core_ok = z2_ok = True
        for t in range(2, m + 3):
            r = Fr(2 * n.get(t - 1, 0) + n.get(t - 2, 0), capref[t])
            if t <= tc:
                eta_need = max(eta_need, r - 1)
                if r > 1 + ETA:
                    core_ok = False
            else:
                z2_max = max(z2_max, r)
                if r > 1 - THETA:
                    z2_ok = False
        rec = {"R": R, "D0": D0, "i_pf": ipf, "d_fe": dfe, "m": m, "tc": tc,
               "npos": len(pos), "HP_ok": hp_ok,
               "a0_minus_d": n.get(0, 0) - dfe,
               "eta_needed": float(max(eta_need, 0)),
               "eta_margin": float(ETA - eta_need),
               "z2_max": float(z2_max), "z2_margin": float(1 - THETA - z2_max),
               "core_ok": core_ok, "z2_ok": z2_ok,
               "ENTRY_ETA_PASS": hp_ok and core_ok and z2_ok}
        rows.append(rec)
        print(f"R={R:6d} D0={D0:5d}: PASS={rec['ENTRY_ETA_PASS']} "
              f"m={m:4d} tc={tc:3d} p=0:{not pos} "
              f"eta_need={rec['eta_needed']:.5f} (<= {float(ETA):.5f}) "
              f"z2max={rec['z2_max']:.4f} (<= {float(1-THETA):.4f})",
              flush=True)
    npass = sum(r["ENTRY_ETA_PASS"] for r in rows)
    print(f"[snapcert] {npass}/{len(rows)} snapshots satisfy "
          f"ENTRY-eta(1/64, 1/4)", flush=True)
    json.dump({"eta": [1, 64], "theta": [1, 4], "rows": rows,
               "n_pass": npass, "n_total": len(rows)},
              open(os.path.join(RESULTS,
                   "amm12592_entry_snapcert_boxeph.json"), "w"), indent=1)
    return rows


# ------------------------------------------------------------------ waitsim
_spec = importlib.util.spec_from_file_location(
    "fastflow", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fastflow = importlib.util.module_from_spec(_spec)
with contextlib.redirect_stdout(io.StringIO()):
    _spec.loader.exec_module(fastflow)
two_G_coeffs = fastflow.two_G_coeffs
initial_junk = fastflow.initial_junk
floor_gamma_star = fastflow.floor_gamma_star


def waitsim():
    """From each snapshot: simulate the (autonomous) wait window and verify
    the E-1 conclusion by direct inspection: support frozen, no death,
    F1 /\ F2 /\ F3 hold at i0 = i_pf + kstar (and record the FIRST row where
    they hold).  Engine semantics verbatim."""
    out = []
    for R, D0 in CONFIGS:
        ipf, dfe, j = load_snap(R, D0)
        en, ed = (1, 32) if D0 * 32 <= R + 31 else (1, 16)
        ps = psi_one(R, en, ed)
        kstar = ps["kstar"]
        i0 = ipf + kstar
        d_prev = dfe
        m_fe = max(j)
        first_ok = None
        support_max = m_fe
        died = False
        for i in range(ipf, min(i0 + 1, R - 1)):
            d = floor_gamma_star(R + i) + D0
            delta = d - d_prev
            assert d + i > R, "wait window must be post-feed"
            # check F1/F2/F3 on the state entering row i
            if j:
                n = {t: -v for t, v in j.items()}
                m = max(j)
                allneg = all(v < 0 for v in j.values())
                F1 = allneg and m + 2 < d
                F2 = n.get(0, 0) <= d - 1
                cap = [2 * comb(d - 1, t) for t in range(m + 3)]
                F3 = all(2 * n.get(t - 1, 0) + n.get(t - 2, 0) <= cap[t]
                         for t in range(2, m + 3))
                if F1 and F2 and F3 and first_ok is None:
                    first_ok = i
            # transport + clamp (verbatim)
            w = {}
            K = (1, 1) if delta == 0 else (1, 2, 1)
            for t, v in j.items():
                for s, ks in enumerate(K):
                    w[t + s] = w.get(t + s, 0) + ks * v
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
                died = True
                break
            j = jn
            d_prev = d
            if j:
                support_max = max(support_max, max(j))
            if not j:
                break
        ok_at_i0 = (first_ok is not None and first_ok <= i0 and not died)
        rec = {"R": R, "D0": D0, "i_pf": ipf, "kstar": kstar, "i0": i0,
               "first_F123_row": first_ok, "offset_pf": (first_ok - ipf)
               if first_ok is not None else None,
               "support_fe": m_fe, "support_max": support_max,
               "support_frozen": support_max <= m_fe,
               "died": died, "E1_consistent": ok_at_i0}
        out.append(rec)
        print(f"R={R:6d} D0={D0:5d}: first F1^F2^F3 at i_pf+"
              f"{rec['offset_pf']}, i0=i_pf+{kstar}, "
              f"support {m_fe}->{support_max} frozen={rec['support_frozen']} "
              f"E1_consistent={ok_at_i0}", flush=True)
    npass = sum(r["E1_consistent"] for r in out)
    print(f"[waitsim] {npass}/{len(out)} consistent", flush=True)
    json.dump(out, open(os.path.join(
        RESULTS, "amm12592_entry_waitsim_boxeph.json"), "w"), indent=1)
    return out


# ----------------------------------------------------------------- feedonly
def feedonly(R, D0, flush_every=256):
    """Run rows 1..i_pf only (the feed phase); write the standard feed-end
    snapshot + a progress ledger.  Engine semantics verbatim
    (certified: amm12592_S_cone_certify_vs_engine_boxeph.json)."""
    t0 = time.time()
    g = two_G_coeffs(R)
    d_prev = floor_gamma_star(R) + D0
    j, junkL1, c0 = initial_junk(R, d_prev)
    minus2 = 1 if c0 == -2 else 0
    led_path = os.path.join(
        RESULTS, f"amm12592_entry_feedonly_R{R}_D0{D0}_boxeph.json")
    for i in range(1, R - 1):
        d = floor_gamma_star(R + i) + D0
        delta = d - d_prev
        if d + i > R:                      # post-feed: snapshot and stop
            snap = {"R": R, "D0": D0, "i_postfeed": i, "d_feedend": d_prev,
                    "entry": {"note": "written by entry_reduction feedonly"},
                    "junk": {str(t): str(v) for t, v in j.items()}}
            json.dump(snap, open(os.path.join(RESULTS,
                f"amm12592_S_cone_feedend_R{R}_D0{D0}_boxeph.json"), "w"))
            el = round(time.time() - t0, 1)
            m = max(j) if j else -1
            print(f"[feedonly] R={R} D0={D0}: i_pf={i} d_fe={d_prev} "
                  f"m={m} nlive={len(j)} minus2={minus2} ({el}s)", flush=True)
            json.dump({"R": R, "D0": D0, "i_postfeed": i, "d_feedend": d_prev,
                       "m": m, "nlive": len(j), "minus2_rows": minus2,
                       "elapsed_s": el, "done": True},
                      open(led_path, "w"))
            return snap
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
        if d + i <= R - 1:
            w[0] = w.get(0, 0) + g[d + i]
        if delta == 1 and d - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[d - 1 + i]
            w[1] = w.get(1, 0) + g[d - 1 + i]
        jn = {}
        c0 = 0
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
        assert d not in jn, ("DEATH in feed phase", i)
        if c0 == -2:
            minus2 += 1
        j = jn
        d_prev = d
        if i % flush_every == 0:
            ts = sorted(j)
            json.dump({"R": R, "D0": D0, "partial_row": i,
                       "front": ts[-1] if ts else None, "nlive": len(j),
                       "elapsed_s": round(time.time() - t0, 1),
                       "done": False}, open(led_path, "w"))
    raise RuntimeError("feed never ended (impossible)")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "psi":
        psi_sweep()
    elif mode == "snapcert":
        snapcert()
    elif mode == "waitsim":
        waitsim()
    elif mode == "feedonly":
        feedonly(int(sys.argv[2]), int(sys.argv[3]))
    else:
        raise SystemExit(f"unknown mode {mode}")
