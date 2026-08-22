"""LANE D3 -- THE GOLDEN TRANSIENT BOUND, part 2: the super-pair death
certificate (boxeph, 2026-08-03).  Exact int arithmetic; no floats anywhere.

Machine side of Lemma G3 (see amm12592-golden-transient-bound-boxeph.md):

  G3 (pair-march death lemma).  Suppose after the clamp of row i0 the junk j
  (degree d) satisfies: L = max supp j >= 4, j_{L-1} and j_L nonzero of the
  same sign s, values b0 = |j_{L-1}|, a0 = |j_L|.  Consider ANY admissible
  continuation (any in-box even choices at every later row).  Along the
  forced diagonal L_{k+1} = L_k + 1 + delta_k, d_{k+1} = d_k + delta_k
  (delta = Beatty word of gamma*), define lower bounds (a_k, b_k) and the
  above-front disturbance majorant psi_k by the recursions

    delta=1:  fl = a_k - 2 psi_k(L+1) - psi_k(L+2)   [load at cell L+2]
              ml = 2 a_k + b_k - psi_k(L+1)          [load at cell L+1]
    delta=0:  fl = a_k - psi_k(L+1)                  [load at cell L+1]
              ml = a_k + b_k                         [load at cell L]
    ends (next degree d' = d + delta, sign s):
              end_s(c) = 2 C(d'-1, c)   if s < 0   (lower box end)
                       = 2 C(d'-1, c-1) if s > 0   (upper box end)
    step OK iff fl > end_s(front cell) and ml > end_s(mate cell); then
              a_{k+1} = fl - end_s(front cell),  b_{k+1} = ml - end_s(mate)
    psi update (new layer cells c > new front):
              psi'(c) = sum_s K_s psi(c - s) + seed(c)
              seed(c) = 2 C(d', c)  [box width: any admissible rule can
                        inject at most this at an in-box cell]  -- seeds ON
                        for the all-continuations theorem; seeds OFF gives
                        the no-injection class (plain rule A provably in it).

  If every step is OK until the front cell reaches d', junk is FORCED at the
  bottom cell there for EVERY admissible continuation: death at row
  i0 + (d_{i0} - L_{i0}) exactly.  (Proof in the note: sign forcing at
  overflow cells + one-directional transport + the layer majorant.)

Modes:
  deaths:   replay near-threshold dying runs (plain/desc1), find the first
            row where the certificate validates (seeds off and on), compare
            predicted vs actual death row, record margins.
  controls: replay adjacent CLOSING runs; verify the closed-form
            precondition (same-sign pair with m > Theta' = sum_{t>=L-1}
            2C(d,t)) NEVER fires (hostile self-test of the lemma).
  grid:     exact threshold table for the closed-form corollary:
            L*(d) = min{L : C(d,L-1) >= sum_{t=L}^{d} C(d,t)} vs 2d/3,
            plus (star-0) <=> 2L >= d+1 identity check.

Outputs: 05-knowledge/results/amm12592_golden_superblock_certificate_boxeph.out
         05-knowledge/results/amm12592_golden_superblock_certificate_boxeph.json
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
OUT = os.path.join(RESULTS,
                   "amm12592_golden_superblock_certificate_boxeph.out")
JS = os.path.join(RESULTS,
                  "amm12592_golden_superblock_certificate_boxeph.json")

_argv = sys.argv[:]
sys.argv = [_argv[0], "noop"]
spec = importlib.util.spec_from_file_location(
    "sweep", os.path.join(HERE, "amm12592_bulkrule_design_sweep_boxeph.py"))
sweep = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(sweep)
sys.argv = _argv
floorgs = sweep.floorgs
two_G_coeffs = sweep.two_G_coeffs
t4_row_load = sweep.t4_row_load
choose_row = sweep.choose_row
VARIANTS = sweep.VARIANTS

RES = {}
if os.path.exists(JS):
    try:
        RES = json.load(open(JS))
    except Exception:
        RES = {}


def log(line):
    print(line, flush=True)
    with open(OUT, "a") as f:
        f.write(line + "\n")


def save():
    json.dump(RES, open(JS, "w"), indent=1)


def end_s(dprime, c, s):
    """Absorption box end on the sign-s side at next degree dprime, cell c."""
    if s < 0:
        return 2 * comb(dprime - 1, c) if c <= dprime - 1 else 0
    return 2 * comb(dprime - 1, c - 1) if 1 <= c <= dprime else 0


def chain_certificate(R, D0, i0, d, L, a0, b0, s, seeds):
    """Forward death certificate from the row-i0 post-clamp state.
    Returns dict: validated (bool), predicted_death_row, steps, margins,
    fail info if not validated.  Everything exact."""
    a, b = a0, b0
    psi = {}
    i = i0
    steps = 0
    min_fmargin_bits = None     # min over steps of bits(fl - Ef) - bits(Ef)
    max_psi_vs_a = None         # max over steps of bits(psi(L+1)) - bits(a)
    while True:
        dn = floorgs(R + i + 1) + D0
        delta = dn - d
        assert delta in (0, 1)
        if delta == 1:
            cf, cm = L + 2, L + 1
            fl = a - 2 * psi.get(L + 1, 0) - psi.get(L + 2, 0)
            ml = 2 * a + b - psi.get(L + 1, 0)
        else:
            cf, cm = L + 1, L
            fl = a - psi.get(L + 1, 0)
            ml = a + b
        if psi:
            p1 = psi.get(L + 1, 0)
            if p1:
                v = p1.bit_length() - a.bit_length()
                max_psi_vs_a = v if max_psi_vs_a is None \
                    else max(max_psi_vs_a, v)
        Ef = end_s(dn, cf, s)
        if not fl > Ef:
            return {"validated": False, "fail_row": i + 1, "steps": steps,
                    "fail": "front", "fl_bits": fl.bit_length() if fl > 0
                    else -((-fl).bit_length() if fl else 0),
                    "Ef_bits": Ef.bit_length()}
        fm = (fl - Ef).bit_length() - max(1, Ef.bit_length())
        min_fmargin_bits = fm if min_fmargin_bits is None \
            else min(min_fmargin_bits, fm)
        if cf >= dn:
            # junk forced at the bottom cell: DEATH at row i+1
            return {"validated": True, "predicted_death_row": i + 1,
                    "steps": steps + 1,
                    "min_front_margin_bits": min_fmargin_bits,
                    "max_psi_vs_a_bits": max_psi_vs_a}
        Em = end_s(dn, cm, s)
        if not ml > Em:
            return {"validated": False, "fail_row": i + 1, "steps": steps,
                    "fail": "mate"}
        a2 = fl - Ef
        b2 = ml - Em
        # layer majorant update: new layer cells c in [cf+1, dn]
        newpsi = {}
        hi_old = max(psi) if psi else L        # old layer top
        hi_new = min(dn, max(hi_old + 1 + delta, cf + 1 + delta))
        for c in range(cf + 1, hi_new + 1):
            if delta == 1:
                v = psi.get(c, 0) + 2 * psi.get(c - 1, 0) + psi.get(c - 2, 0)
            else:
                v = psi.get(c, 0) + psi.get(c - 1, 0)
            if seeds:
                v += 2 * comb(dn, c)
            if v:
                newpsi[c] = v
        psi = newpsi
        a, b, L, d = a2, b2, cf, dn
        i += 1
        steps += 1


def theta_prime(d, L):
    """Theta' = sum_{t=L-1}^{d} 2 C(d,t)  (exact)."""
    tot = 0
    c = comb(d, L - 1)
    for t in range(L - 1, d + 1):
        tot += 2 * c
        c = c * (d - t) // (t + 1) if t < d else 0
    return tot


def replay(R, D0, variant, certify_deaths=True, max_attempts=None):
    """Replay via the D1 generalized clamp; per-row precondition check and
    (for dying runs) chain-certificate validation."""
    cfg = VARIANTS[variant]
    g = two_G_coeffs(R)
    j = {}
    d_prev = None
    rec = {"R": R, "D0": D0, "variant": variant, "outcome": None,
           "first_precond_row": None, "precond_rows": 0,
           "first_valid_noseed": None, "first_valid_seeded": None,
           "attempts_noseed": 0, "attempts_seeded": 0}
    t0 = time.time()
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
            rec["outcome"] = {"CLOSED_capture": i}
            break
        jn, us, c0, junkL1, diag = choose_row(w, d, dn, S, cfg,
                                              record_u=False)
        died = d in jn
        if died:
            rec["outcome"] = {"DIE": i}
            rec["death_row"] = i
        # ---- per-row analysis of the post-clamp state
        if jn and not died:
            L = max(jn)
            if L >= 4 and (L - 1) in jn and \
                    (jn[L] > 0) == (jn[L - 1] > 0):
                a0, b0 = abs(jn[L]), abs(jn[L - 1])
                m = min(a0, b0)
                th = theta_prime(d, L)
                if m > th:
                    rec["precond_rows"] += 1
                    if rec["first_precond_row"] is None:
                        rec["first_precond_row"] = i
                        rec["precond_state"] = {
                            "i": i, "d": d, "L": L, "gap": d - L,
                            "m_bits": m.bit_length(),
                            "theta_bits": th.bit_length()}
                if certify_deaths and m > th:
                    s = 1 if jn[L] > 0 else -1
                    if rec["first_valid_noseed"] is None and \
                            (max_attempts is None or
                             rec["attempts_noseed"] < max_attempts):
                        rec["attempts_noseed"] += 1
                        r = chain_certificate(R, D0, i, d, L, a0, b0, s,
                                              seeds=False)
                        if r["validated"]:
                            rec["first_valid_noseed"] = i
                            rec["chain_noseed"] = r
                            rec["chain_noseed"]["gap_at_valid"] = d - L
                    if rec["first_valid_seeded"] is None and \
                            (max_attempts is None or
                             rec["attempts_seeded"] < max_attempts):
                        rec["attempts_seeded"] += 1
                        r = chain_certificate(R, D0, i, d, L, a0, b0, s,
                                              seeds=True)
                        if r["validated"]:
                            rec["first_valid_seeded"] = i
                            rec["chain_seeded"] = r
                            rec["chain_seeded"]["gap_at_valid"] = d - L
        j = jn
        d_prev = d
        if died:
            break
        if not j and d + i > R - 1:
            rec["outcome"] = {"CLOSED_capture": i}
            break
    rec["elapsed_s"] = round(time.time() - t0, 1)
    return rec


def mode_deaths(runs):
    out = []
    for R, D0, var in runs:
        r = replay(R, D0, var, certify_deaths=True)
        out.append(r)
        cn = r.get("chain_noseed") or {}
        cs = r.get("chain_seeded") or {}
        actual = r.get("death_row")
        log(f"DEATH {var} R={R} D0={D0}: outcome={r['outcome']} "
            f"first_precond={r['first_precond_row']} "
            f"({r.get('precond_state')})")
        log(f"   noseed: first_valid={r['first_valid_noseed']} "
            f"predicted_death={cn.get('predicted_death_row')} "
            f"actual={actual} match={cn.get('predicted_death_row') == actual} "
            f"min_front_margin_bits={cn.get('min_front_margin_bits')} "
            f"attempts={r['attempts_noseed']}")
        log(f"   seeded: first_valid={r['first_valid_seeded']} "
            f"predicted_death={cs.get('predicted_death_row')} "
            f"match={cs.get('predicted_death_row') == actual} "
            f"min_front_margin_bits={cs.get('min_front_margin_bits')} "
            f"max_psi_vs_a_bits={cs.get('max_psi_vs_a_bits')} "
            f"attempts={r['attempts_seeded']}  ({r['elapsed_s']}s)")
        RES.setdefault("deaths", []).append(r)
        save()
    return out


def mode_controls(runs):
    for R, D0, var in runs:
        r = replay(R, D0, var, certify_deaths=False)
        ok = (r["first_precond_row"] is None and
              r["outcome"] is not None and
              "CLOSED_capture" in r["outcome"])
        log(f"CONTROL {var} R={R} D0={D0}: outcome={r['outcome']} "
            f"precond_rows={r['precond_rows']} "
            f"first_precond={r['first_precond_row']} -> "
            f"{'OK (never fires)' if ok else 'ATTENTION'} "
            f"({r['elapsed_s']}s)")
        RES.setdefault("controls", []).append(r)
        save()


def mode_grid():
    """L*(d) table + (star-0) identity check."""
    recs = []
    for d in (60, 100, 153, 306, 612, 1262, 2538, 5090):
        # L*(d) = min L with C(d, L-1) >= sum_{t=L}^{d} C(d,t)
        Ls = None
        for L in range(d // 2, d):
            c = comb(d, L - 1)
            tail = 0
            cc = comb(d, L)
            failed = False
            for t in range(L, d + 1):
                tail += cc
                if tail > c:
                    failed = True          # full tail can only be larger
                    break
                cc = cc * (d - t) // (t + 1) if t < d else 0
            if not failed and c >= tail:
                Ls = L
                break
        # (star-0): C(d,L-1) >= C(d,L)  <=>  2L >= d+1  (exact check on window)
        s0 = all((comb(d, L - 1) >= comb(d, L)) == (2 * L >= d + 1)
                 for L in range(max(1, d // 2 - 3), d // 2 + 5))
        recs.append({"d": d, "Lstar": Ls, "Lstar_minus_2d3": Ls - (2 * d) // 3,
                     "star0_iff_2L_ge_d+1": s0})
        log(f"GRID d={d}: L*={Ls}  L*-2d/3={Ls - (2*d)//3}  "
            f"tau*={Ls}/{d}={Ls/d:.5f}  star0_identity={s0}")
    RES["grid"] = recs
    save()


if __name__ == "__main__":
    t0 = time.time()
    log("=" * 72)
    log("LANE D3 script B: super-pair death certificate (exact ints)")
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("grid", "all"):
        mode_grid()
    if mode in ("deaths", "all"):
        runs = [(512, 4, "plain"), (1024, 14, "plain"), (1024, 13, "desc1"),
                (2048, 37, "plain")]
        if mode == "deaths" and len(sys.argv) > 2:
            R, D0, var = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
            runs = [(R, D0, var)]
        mode_deaths(runs)
    if mode in ("controls", "all"):
        runs = [(512, 5, "plain"), (1024, 15, "plain"), (1024, 14, "desc1"),
                (2048, 38, "plain")]
        if mode == "controls" and len(sys.argv) > 2:
            R, D0, var = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
            runs = [(R, D0, var)]
        mode_controls(runs)
    log(f"done in {round(time.time()-t0,1)}s")
    save()
