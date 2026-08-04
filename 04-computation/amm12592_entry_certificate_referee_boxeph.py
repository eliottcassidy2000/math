"""LANE F1 -- INDEPENDENT referee for the new-scale ENTRY certificates.

For a given (R, D0): load the scan ledger
`amm12592_entry_scan_R{R}_D0{D0}_boxeph.json`, take the SAVED feed-end
state (junk entering row i_pf), and, with an implementation independent
of both scan engines (fresh clamp via math.comb per cell; fresh floor
profile from the imported certified Lucas/Fibonacci engine only):

  1. re-evolve the post-feed rows i_pf .. i_fc (all feed-free: d + i > R
     asserted per row);
  2. at i_fc re-verify F1 (negativity, support, m + 2 < d), F2, F3
     exactly;
  3. recompute the F4 clocks by the mul//div staircase (independent of
     the scanner's Pascal-add staircase) and compare (drain, Tmax,
     capture_ub, budget_margin) against the ledger's fc record;
  4. verify the T8 identity a_0 = (R-2) - 2*minus2 at i_pf is consistent
     with the ledger, and that the certificate row implies S(R) via
     Theorem S-cone-fc (capture_ub <= R - 2).

Output -> 05-knowledge/results/amm12592_entry_certreferee_R{R}_D0{D0}_boxeph.json
All exact ints; no floats in any decision.
"""
import sys, os, json, time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

import importlib.util as _ilu, io as _io, contextlib as _ctx
_spec = _ilu.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
_iref = _ilu.module_from_spec(_spec)
with _ctx.redirect_stdout(_io.StringIO()):
    _spec.loader.exec_module(_iref)
floor_gamma_star = _iref.floor_gamma_star


def clamp_row(j, dp):
    """One post-feed row at new degree dp (junk = load - clamp), fresh
    comb per cell (independent implementation)."""
    w = {}
    for t, v in j.items():
        w[t] = w.get(t, 0) + v
    # kernel breadth from delta encoded by caller: caller passes K
    return w


def evolve_row(j, dp, delta):
    K = (1, 1) if delta == 0 else (1, 2, 1)
    w = {}
    for t, v in j.items():
        for s, c in enumerate(K):
            w[t + s] = w.get(t + s, 0) + c * v
    jn = {}
    for t, v in sorted(w.items()):
        if t == 0:
            lo, hi = -2, 0
        else:
            lo, hi = -2 * comb(dp - 1, t), 2 * comb(dp - 1, t - 1)
        u = min(hi, max(lo, v))
        assert v % 2 == 0 and (v - u) % 2 == 0, ("parity", t)
        if v != u:
            jn[t] = v - u
    return jn


def referee(R, D0):
    t0 = time.time()
    led = json.load(open(os.path.join(
        RESULTS, f"amm12592_entry_scan_R{R}_D0{D0}_boxeph.json")))
    fe, fc = led["feedend"], led["fc"]
    assert fe["junk"] is not None, "no snapshot saved (m > 400?)"
    i_pf, i_fc = led["i_pf"], led["i_fc"]
    j = {int(t): int(v) for t, v in fe["junk"].items()}
    checks = {}
    # T8 consistency at i_pf
    a0_pf = -j.get(0, 0)
    checks["T8_identity"] = (a0_pf == (R - 2) - 2 * fe["minus2"])
    # evolve i_pf .. i_fc - 1 (state entering row i is j; certificate is
    # the state entering row i_fc)
    d_prev = fe["d_fe"]
    for i in range(i_pf, i_fc):
        d = floor_gamma_star(R + i) + D0
        assert d + i > R, ("row not post-feed", i)
        j = evolve_row(j, d, d - d_prev)
        d_prev = d
    d = floor_gamma_star(R + i_fc) + D0
    assert d + i_fc > R
    # F1/F2/F3 at i_fc
    m = max(j)
    a = {t: -v for t, v in j.items()}
    checks["F1_neg"] = all(v < 0 for v in j.values())
    checks["F1_front"] = (m + 2 < d and m == fc["m"])
    a0 = a.get(0, 0)
    checks["F2"] = (a0 <= d - 1 and a0 == fc["a0"])
    capfc = {t: 2 * comb(d - 1, t) for t in range(m + 3)}
    checks["F3"] = all(2 * a.get(t - 1, 0) + a.get(t - 2, 0) <= capfc[t]
                       for t in range(2, m + 3))
    # independent clocks (mul//div staircase)
    drain = -((-a0) // 2)
    alive = [t for t in range(1, m + 1) if a.get(t, 0)]
    need = {t: a[t] for t in alive}
    cum = {t: 0 for t in alive}
    Tt = {}
    dd = d
    capk = dict(capfc)
    k = 0
    while need and k < R - 2 - i_fc:
        k += 1
        dn = floor_gamma_star(R + i_fc + k) + D0
        dl = dn - dd
        if dl:
            for t in range(m + 3):
                capk[t] = capk[t] * dd // (dd - t)
        a0k = max(0, a0 - 2 * (k - 1))
        for t in list(need):
            if t == 1:
                dec = 2 * (dn - 1) - (1 + dl) * a0k
                if dec > 0:
                    cum[1] += dec
            else:
                dec = capk[t] - capfc[t]
                if dec > 0:
                    cum[t] += dec
            if cum[t] >= need[t]:
                Tt[t] = k
                del need[t]
        dd = dn
    checks["clocks_finish"] = not need
    Tmax = max(Tt.values()) if Tt else 0
    ub = i_fc + max(drain, Tmax)
    checks["drain_eq"] = (drain == fc["drain"])
    checks["Tmax_eq"] = (Tmax == fc["Tmax"])
    checks["capture_ub_eq"] = (ub == fc["capture_ub"])
    checks["budget_margin_eq"] = ((R - 2) - ub == fc["budget_margin"])
    checks["F4_budget"] = (ub <= R - 2)
    ok = all(checks.values())
    rec = {"R": R, "D0": D0, "i_pf": i_pf, "i_fc": i_fc, "checks": checks,
           "S_R_established": ok,
           "elapsed_s": round(time.time() - t0, 1)}
    json.dump(rec, open(os.path.join(
        RESULTS, f"amm12592_entry_certreferee_R{R}_D0{D0}_boxeph.json"),
        "w"), indent=1)
    print(f"[cert-referee] R={R} D0={D0}: "
          f"{'ALL-PASS -- S(R) established' if ok else 'FAIL ' + str(checks)}"
          f" ({rec['elapsed_s']}s)", flush=True)
    return ok


if __name__ == "__main__":
    referee(int(sys.argv[1]), int(sys.argv[2]))
