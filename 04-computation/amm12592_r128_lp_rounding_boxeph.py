"""AMM 12592 / ANGLE 3: transportation-LP form + parity rounding at R = 128.

Nonnegative transportation form (05-knowledge/results/
amm12592-epoch-closure-nonnegative-transportation-form-boxeph.md):
substituting delta_{i,k} = binom(d_i,k) - 2 f_{i,k} turns THM-3002 (*)

    sum_{i<R} x^i Delta_i = (1-x)^{R-1},  |delta| <= binom, delta == binom (2)

into the capacitated transportation problem

    (**)  sum_{i<R} x^i F_i = T_R := [ (1-x^R)/(1-x) - (1-x)^{R-1} ] / 2,
          F_i = sum_k f_{i,k} B_{d_i,k},   f_{i,k} integer in [0, binom(d_i,k)],

where PARITY IS FREE (any integer f gives delta == binom mod 2 automatically)
and the Lucas box |delta| <= binom is exactly the capacity 0 <= f <= binom.

This file:
  PHASE A (claim check, run FIRST): translate the verified R=8 witness (and
    R=16/32/64/128-direct) to f-coordinates via f = (binom - delta)/2; confirm
    integrality + bounds + the exact (**) identity against T_R in Z[x]; confirm
    T_R in Z[x] at R = 128 (all binom(127,j) odd -- Lucas, 127 = 2^7 - 1).
  PHASE B (relaxation): the translated R=128 point is an EXACT rational
    (indeed integer) point of the real relaxation -- feasibility certificate.
    Emit the box-saturation / slack profile (which cells sit on capacity).
  PHASE C (LP-guided rounding search): independent construction pipeline.
    Per row of the residual recursion  p sigma_i = sigma_{i-1} - Delta_i,
    sigma_{-1} = q^{R-1}  (delta-coordinates; the affine image of (**)), the
    per-cell REAL LP optimum under the triangular Bernstein elimination is the
    box clamp  w* = clip(res_j, +-cap).  We take w* and round to the forced
    parity class, branching over the rounding direction / small even offsets
    on the first `nbranch` cells (rounding debt absorbed by the residual
    recursion = carry scheme), deterministic nearest-parity clamp elsewhere,
    beam ranked by (residual degree, residual L1) ('l1deg', the R=64 winner),
    with a wide 2-row endgame (row R-2 offsets up to +-6 on 3 cells,
    row R-1 exact Bernstein decode).  Every landed witness is verified by the
    search-independent exact referee before saving.

Exact integer arithmetic only; floats never touch any verification path.

Usage:
  python3 amm12592_r128_lp_rounding_boxeph.py --phase ab      # claims + slack
  python3 amm12592_r128_lp_rounding_boxeph.py --phase ladder  # R=8..64 pipeline validation
  python3 amm12592_r128_lp_rounding_boxeph.py --phase r128    # R=128 hunt
"""
import argparse, json, os, sys, time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(HERE, "..", "05-knowledge", "results"))

# gamma* = log(phi)/log(sqrt 5); 16-digit proxy, floors equal true gamma*
# floors for all m <= 512 (proven upstream).
GS = (5979874356654402, 10**16)

def profile(R):
    return [(GS[0] * (R + i)) // GS[1] for i in range(R)]

_QK = [[1]]
def qk(k):
    while len(_QK) <= k:
        r = _QK[-1]
        nr = [0] * (len(r) + 1)
        for i, c in enumerate(r):
            nr[i] += c
            nr[i + 1] -= c
        _QK.append(nr)
    return _QK[k]

def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a

# ---------------------------------------------------------------------------
# exact verification (delta side), independent of the search
# ---------------------------------------------------------------------------
def admissible(delta, dd):
    if len(delta) - 1 > dd:
        return False
    return all(abs(v) <= comb(dd, k) and (v - comb(dd, k)) % 2 == 0
               for k, v in enumerate(delta))

def epoch_identity(R, blocks, d):
    acc = [0] * (R + max(d) + 2)
    for i, de in enumerate(blocks):
        for k, v in enumerate(de):
            if v:
                t = qk(k)
                off = i + d[i] - k
                for s in range(k + 1):
                    acc[off + s] += v * t[s]
    return trim(acc) == trim(qk(R - 1))

# ---------------------------------------------------------------------------
# PHASE A: f-coordinates, T_R, exact (**) verification
# ---------------------------------------------------------------------------
def T_poly(R):
    """T_R = [ sum_{j<R} x^j - (1-x)^{R-1} ] / 2 as an exact integer vector.
    Returns (coeffs, integral) -- integral fails iff some binom(R-1,j) even."""
    q = qk(R - 1)
    co, ok = [], True
    for j in range(R):
        num = 1 - (q[j] if j < len(q) else 0)
        if num % 2:
            ok = False
        co.append(num // 2)
    return co, ok

def to_f(blocks, d):
    """delta -> f = (binom - delta)/2. Returns (fblocks, integral, bounded)."""
    fb, integral, bounded = [], True, True
    for i, de in enumerate(blocks):
        dd = d[i]
        row = []
        for k in range(dd + 1):
            b = comb(dd, k)
            v = de[k] if k < len(de) else 0
            num = b - v
            if num % 2:
                integral = False
            f = num // 2
            if f < 0 or f > b:
                bounded = False
            row.append(f)
        fb.append(row)
    return fb, integral, bounded

def star2_identity(R, fblocks, d):
    """Exact check: sum_i x^i F_i == T_R in Z[x]."""
    acc = [0] * (R + max(d) + 2)
    for i, fr in enumerate(fblocks):
        for k, v in enumerate(fr):
            if v:
                t = qk(k)
                off = i + d[i] - k
                for s in range(k + 1):
                    acc[off + s] += v * t[s]
    T, ok = T_poly(R)
    return ok and trim(acc) == trim(T)

def saturation_profile(fblocks, d):
    """Per-row counts of capacity-saturated cells (f=0 or f=binom) in the
    transportation polytope; the slack anatomy of the integer point."""
    rows = []
    tot = sat = 0
    for i, fr in enumerate(fblocks):
        dd = d[i]
        n0 = sum(1 for k, f in enumerate(fr) if f == 0)
        nC = sum(1 for k, f in enumerate(fr) if f == comb(dd, k))
        # min slack among interior-capable cells
        ms = min(min(f, comb(dd, k) - f) for k, f in enumerate(fr))
        rows.append({"i": i, "d": dd, "cells": dd + 1, "at0": n0,
                     "atcap": nC, "minslack": ms})
        tot += dd + 1
        sat += n0 + nC
    return rows, tot, sat

def phase_ab(out):
    P = lambda *a: print(*a, file=out, flush=True)
    P("=" * 74)
    P("PHASE A: f-coordinate (transportation) claim checks")
    P("=" * 74)
    small = json.load(open(os.path.join(
        HERE, "amm12592_floor_witnesses_R8_R16_R32.json")))
    wit = {w["R"]: w for w in small}
    for name, fn in [("R64 slim", "amm12592_floor_witness_R64.json"),
                     ("R128 direct", "amm12592_floor_witness_R128_direct.json")]:
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            w = json.load(open(p))
            wit[w["R"]] = w
    f128 = None
    for R in sorted(wit):
        w = wit[R]
        d = w["profile"]
        assert d == profile(R), f"profile mismatch at R={R}"
        blocks = w["blocks"]
        a_ok = all(admissible(blocks[i], d[i]) for i in range(R))
        i_ok = epoch_identity(R, blocks, d)
        fb, f_int, f_bnd = to_f(blocks, d)
        T, T_int = T_poly(R)
        s_ok = star2_identity(R, fb, d)
        P(f"R={R:4d}: delta-witness admissible={a_ok} identity={i_ok} | "
          f"f=(binom-delta)/2 integral={f_int} in-bounds={f_bnd} | "
          f"T_R integral={T_int} | (**) sum x^i F_i == T_R exact: {s_ok}")
        assert a_ok and i_ok and f_int and f_bnd and T_int and s_ok
        if R == 128:
            f128 = (fb, d)
    P("")
    P("CLAIM CONFIRMED: the parity condition is absorbed exactly -- every")
    P("verified delta-witness maps to an INTEGER point of the capacitated")
    P("transportation polytope 0 <= f <= binom with sum_i x^i F_i = T_R, and")
    P("back.  Any integer f in the polytope is automatically parity-correct.")
    P("")
    q = qk(127)
    P(f"T_128 integrality: all binom(127,j) odd (127 = 2^7 - 1, Lucas): "
      f"{all(c % 2 for c in q)}")
    if f128 is None:
        P("R128 direct witness absent; PHASE B skipped")
        return None
    P("")
    P("=" * 74)
    P("PHASE B: real relaxation at R = 128 -- exact feasibility + slack")
    P("=" * 74)
    fb, d = f128
    rows, tot, sat = saturation_profile(fb, d)
    P("The translated R=128 point is itself an exact integer point of the")
    P("real relaxation: RELAXATION FEASIBLE (certificate = the point; no")
    P("floating point anywhere).  Box-saturation anatomy of that point:")
    P(f"  total cells {tot}, saturated (f=0 or f=cap) {sat} "
      f"({100.0*sat/tot:.2f}%)")
    worst = sorted(rows, key=lambda r: (r["at0"] + r["atcap"]), reverse=True)[:8]
    P("  most saturated rows (i, d_i, cells, at0, atcap): " +
      ", ".join(f"({r['i']},{r['d']},{r['cells']},{r['at0']},{r['atcap']})"
                for r in worst))
    nz_slack = sorted(rows, key=lambda r: r["minslack"])[:8]
    P("  min per-row interior slack (i, minslack): " +
      ", ".join(f"({r['i']},{r['minslack']})" for r in nz_slack))
    unit_sat = sum(1 for r in rows if r["minslack"] == 0)
    P(f"  rows with some cell ON its bound: {unit_sat}/128 (the forced top")
    P("  cell f_(i,d_i) in {0,1} = Delta_i(0) = +-1 saturates EVERY row --")
    P("  the transportation polytope has no interior in that direction;")
    P("  integer infeasibility can never be excluded by margin alone).")
    pt = {"R": 128, "profile": d, "f_blocks": fb,
          "form": "sum_i x^i sum_k f_{i,k} B_{d_i,k} = T_128, "
                  "0<=f<=binom(d_i,k)",
          "T_definition": "[ (1-x^128)/(1-x) - (1-x)^127 ] / 2",
          "verified": True,
          "provenance": "f = (binom - delta)/2 from "
                        "amm12592_floor_witness_R128_direct.json; exact (**) "
                        "identity + bounds re-verified here",
          "saturation": {"total_cells": tot, "saturated": sat,
                         "rows": rows}}
    fp = os.path.join(HERE, "amm12592_r128_transportation_point_boxeph.json")
    json.dump(pt, open(fp, "w"))
    P(f"  transportation-form integer point saved: {fp}")
    return fp

# ---------------------------------------------------------------------------
# PHASE C: LP-clamp + parity-rounding beam
# ---------------------------------------------------------------------------
def lp_step_candidates(sigma, d, nbranch, maxoff, sign_branch=False):
    """One residual row at block degree d.  Returns list of
    (deltas, sigma_next, debt, overshoot).  Cell processing order
    k = d-1 .. 0 handles res positions j = d-k = 1..d triangularly; the
    real per-cell LP optimum is w* = clip(res_j, +-cap); cell j=1 is forced
    to leave a unit constant (sigma_next[0] = +-1, else the next forced top
    cell breaks parity); cells j = 2..nbranch branch over the parity
    roundings of w* +- even offsets <= maxoff; the rest take the
    deterministic nearest-parity clamp (downward bias, as in the R=64
    winner).  debt = sum |w - w*| (rounding mass absorbed by the residual);
    overshoot = sum max(0, |res_j| - cap) (capacity clipping = the
    redistribution obstruction)."""
    if not sigma or abs(sigma[0]) != 1:
        return []
    L = max(len(sigma), d + 1)
    res0 = list(sigma) + [0] * (L - len(sigma))
    v = sigma[0]
    t = qk(d)
    for s in range(min(d + 1, L)):
        res0[s] -= v * t[s]
    de0 = [0] * (d + 1)
    de0[d] = v
    # branch prefix: cells cellidx = 0 .. nbranch-1  (k = d-1 .. d-nbranch)
    prefix = [(res0, de0, 0, 0)]
    for cellidx in range(min(nbranch, d)):
        k = d - 1 - cellidx
        j = d - k
        cap = comb(d, k)
        tq = qk(k)
        nxt = []
        for res, de, debt, over in prefix:
            rj = res[j]
            wstar = -cap if rj < -cap else (cap if rj > cap else rj)
            ov = abs(rj) - cap if abs(rj) > cap else 0
            if cellidx == 0:
                cands = [w for w in (rj - 1, rj + 1)
                         if -cap <= w <= cap and (w - cap) % 2 == 0]
            else:
                if (wstar - cap) % 2 == 0:
                    base = (wstar,)
                else:
                    base = (wstar - 1, wstar + 1)
                cs = set()
                for b in base:
                    for o in range(-maxoff, maxoff + 1, 2):
                        w = b + o
                        if -cap <= w <= cap and (w - cap) % 2 == 0:
                            cs.add(w)
                cands = sorted(cs)
            for w in cands:
                r2 = list(res)
                if w:
                    off = d - k
                    for s in range(k + 1):
                        r2[off + s] -= w * tq[s]
                d2 = list(de)
                d2[k] = w
                nxt.append((r2, d2, debt + abs(w - wstar), over + ov))
        prefix = nxt
        if not prefix:
            return []
    # deterministic tail down to cell k=1; cell k=0 (the sign cell
    # delta_{i,0} = Delta_i(1) = +-1, cap 1, parity odd) is BRANCHED BOTH
    # ways: it alone drives the evaluation-at-1 ballot walk
    # sigma_i(1) = sigma_{i-1}(1) -+ 1, whose budget |sigma_i(1)| <= R-1-i
    # is a hard LP cut (see solve_lp_round) that the local clamp cannot see.
    out = []
    for res, de, debt, over in prefix:
        for k in range(d - 1 - nbranch, 0, -1):
            j = d - k
            cap = comb(d, k)
            rj = res[j]
            w = -cap if rj < -cap else (cap if rj > cap else rj)
            if abs(rj) > cap:
                over += abs(rj) - cap
            if (w - cap) % 2:
                w2 = w - 1 if w - 1 >= -cap else w + 1
                debt += abs(w2 - (rj if -cap <= rj <= cap else w))
                w = w2
            de[k] = w
            if w:
                tq = qk(k)
                off = d - k
                for s in range(k + 1):
                    res[off + s] -= w * tq[s]
        if res[0] != 0:
            continue
        if sign_branch:
            signs = (1, -1)
        else:
            # deterministic residual clamp (the R=64 winner's rule):
            # w0 = parity-adjusted clip of res[d] to +-1, downward bias
            rd = res[d]
            w = -1 if rd < -1 else (1 if rd > 1 else rd)
            if (w - 1) % 2:
                w = w - 1 if w - 1 >= -1 else w + 1
            signs = (w,)
        for w0 in signs:
            r2 = list(res)
            r2[d] -= w0
            d2 = list(de)
            d2[0] = w0
            out.append((d2, trim(r2[1:]),
                        debt + abs(w0 - max(-1, min(1, res[d]))), over))
    return out

def final_decode(sigma, d):
    if len(sigma) - 1 > d:
        return None
    res = list(sigma) + [0] * (d + 1 - len(sigma))
    de = [0] * (d + 1)
    for k in range(d, -1, -1):
        cap = comb(d, k)
        want = res[d - k]
        if abs(want) > cap or (want - cap) % 2:
            return None
        de[k] = want
        if want:
            t = qk(k)
            off = d - k
            for s in range(k + 1):
                res[off + s] -= want * t[s]
    return de if not trim(res) else None

def solve_lp_round(d, beam=400, nbranch=2, maxoff=2, end_nbranch=3,
                   end_maxoff=6, log=None, ckpt_path=None, ckpt_every=12):
    """LP-clamp + parity-rounding beam for sum_i x^i Delta_i == (1-x)^{R-1}
    at profile d.  Returns (blocks, msg, stats)."""
    R = len(d)
    t0 = time.time()
    states = [((), tuple(qk(R - 1)), (), ())]   # (acc, sigma, debts, overs)
    stats = {"rows": []}
    for i in range(R - 2):
        nxt = []
        budget = R - 1 - i          # |sigma_i(1)| <= # remaining rows
        for acc, sig, debts, overs in states:
            for de, ns, debt, over in lp_step_candidates(
                    list(sig), d[i], nbranch, maxoff,
                    sign_branch=(i >= R - 32)):
                if not ns or abs(ns[0]) != 1:
                    continue
                if abs(sum(ns)) > budget:   # evaluation-at-1 ballot cut
                    continue
                nxt.append((acc + (tuple(de),), tuple(ns),
                            debts + (debt,), overs + (over,)))
        if not nxt:
            return None, f"died at row {i}", stats
        nxt.sort(key=lambda st: (len(st[1]), sum(abs(x) for x in st[1]),
                                 abs(sum(st[1]))))
        # dedup + diversity cap: the k=0 sign branch produces "twins" that
        # differ only in the top few coefficients; cap each low-part group
        # at 2 so twins cannot flood the beam and starve trajectory
        # diversity (the row-7 death mode of the all-branch run).
        seen, gcount, uniq = set(), {}, []
        for st in nxt:
            if st[1] in seen:
                continue
            seen.add(st[1])
            g = (len(st[1]), st[1][:len(st[1]) - 4])
            c = gcount.get(g, 0)
            if c >= 2:
                continue
            gcount[g] = c + 1
            uniq.append(st)
        states = uniq[:beam]
        b = states[0][1]
        row = {"i": i, "states": len(states), "bestdeg": len(b) - 1,
               "bestL1": sum(abs(x) for x in b), "bestev1": sum(b),
               "bestdebt": states[0][2][-1], "bestover": states[0][3][-1],
               "t": round(time.time() - t0, 1)}
        stats["rows"].append(row)
        if log:
            print(f"  row {i:3d}: states={row['states']:5d} "
                  f"deg={row['bestdeg']:3d} L1={row['bestL1']} "
                  f"ev1={row['bestev1']} "
                  f"debt={row['bestdebt']} over={row['bestover']} "
                  f"t={row['t']}s", file=log, flush=True)
        if ckpt_path and (i % ckpt_every == 0 or i >= R - 16):
            json.dump({"row": i,
                       "states": [[list(map(list, st[0])), list(st[1]),
                                   list(st[2]), list(st[3])]
                                  for st in states[:100]],
                       "stats": stats},
                      open(ckpt_path, "w"))
    # wide 2-row endgame
    da, db = d[R - 2], d[R - 1]
    tried = 0
    for acc, sig, debts, overs in states:
        for de, ns, debt, over in lp_step_candidates(
                list(sig), da, end_nbranch, end_maxoff, sign_branch=True):
            if not ns or abs(ns[0]) != 1:
                continue
            tried += 1
            fin = final_decode(list(ns), db)
            if fin is not None:
                stats["endgame_decodes"] = tried
                stats["winner_debts"] = list(debts) + [debt]
                stats["winner_overs"] = list(overs) + [over]
                return [list(b) for b in acc] + [de, fin], "SOLVED", stats
    stats["endgame_decodes"] = tried
    return None, f"beam {beam} exhausted at 2-row endgame", stats

def run_ladder(out, beams=(200, 200, 300, 400)):
    P = lambda *a: print(*a, file=out, flush=True)
    P("")
    P("=" * 74)
    P("PHASE C-ladder: pipeline validation at R = 8, 16, 32, 64")
    P("=" * 74)
    ok_all = True
    for R, bm in zip((8, 16, 32, 64), beams):
        d = profile(R)
        t0 = time.time()
        sol, msg, stats = solve_lp_round(d, beam=bm)
        dt = time.time() - t0
        if sol is None:
            P(f"R={R:3d} beam={bm}: {msg}  ({dt:.1f}s)")
            ok_all = False
            continue
        a_ok = all(admissible(sol[i], d[i]) for i in range(R))
        i_ok = epoch_identity(R, sol, d)
        fb, fi, fbnd = to_f(sol, d)
        s_ok = star2_identity(R, fb, d)
        wd = stats.get("winner_debts", [])
        wo = stats.get("winner_overs", [])
        P(f"R={R:3d} beam={bm}: {msg} in {dt:.1f}s; verify adm={a_ok} "
          f"id={i_ok} (**)={s_ok}; winner debt total={sum(wd)} "
          f"max={max(wd) if wd else 0}, overshoot total={sum(wo)} "
          f"max={max(wo) if wo else 0}")
        assert a_ok and i_ok and s_ok
    return ok_all

def run_r128(out, beam, tag=""):
    P = lambda *a: print(*a, file=out, flush=True)
    P("")
    P("=" * 74)
    P(f"PHASE C-r128: LP-rounding hunt at R = 128, beam={beam} {tag}")
    P("=" * 74)
    R = 128
    d = profile(R)
    ck = os.path.join(HERE, "amm12592_r128_lp_ckpt_boxeph.json")
    t0 = time.time()
    sol, msg, stats = solve_lp_round(d, beam=beam, log=out, ckpt_path=ck)
    dt = time.time() - t0
    P(f"  -> {msg}  ({dt:.1f}s, endgame decodes={stats.get('endgame_decodes')})")
    if sol is None:
        return False
    a_ok = all(admissible(sol[i], d[i]) for i in range(R))
    i_ok = epoch_identity(R, sol, d)
    fb, fi, fbnd = to_f(sol, d)
    s_ok = star2_identity(R, fb, d)
    P(f"  EXACT VERIFY: admissible={a_ok} epoch-identity={i_ok} "
      f"f-integral={fi} f-bounded={fbnd} (**)-identity={s_ok}")
    assert a_ok and i_ok and fi and fbnd and s_ok
    wd, wo = stats["winner_debts"], stats["winner_overs"]
    P(f"  winner rounding-debt profile: total={sum(wd)} max={max(wd)} "
      f"rows>0: {sum(1 for x in wd if x)}")
    P(f"  winner overshoot profile: total={sum(wo)} max={max(wo)} "
      f"last row>0: {max((j for j, x in enumerate(wo) if x), default=-1)}")
    wp = os.path.join(HERE, "amm12592_floor_witness_R128_lp.json")
    json.dump({"R": R, "profile": d, "blocks": sol, "verified": True,
               "H": [1],
               "source_label": "gamma* floor (LP-clamp + parity-rounding "
                               "beam, Angle 3, independent of the direct "
                               "beam witness)",
               "search": {"beam": beam, "nbranch": 2, "maxoff": 2,
                          "end_nbranch": 3, "end_maxoff": 6,
                          "rank": "l1deg", "deterministic": True},
               "debt_profile": wd, "overshoot_profile": wo},
              open(wp, "w"))
    P(f"  WITNESS WRITTEN: {wp}")
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="ab",
                    choices=["ab", "ladder", "r128", "all"])
    ap.add_argument("--beam", type=int, default=400)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = open(args.out, "a") if args.out else sys.stdout
    if args.phase in ("ab", "all"):
        phase_ab(out)
    if args.phase in ("ladder", "all"):
        run_ladder(out)
    if args.phase in ("r128", "all"):
        ok = run_r128(out, args.beam)
        sys.exit(0 if ok else 3)

if __name__ == "__main__":
    main()
