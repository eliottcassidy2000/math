"""AMM 12592 / ANGLE 3 endgame widener: replay the 2-row completion from the
row-125 beam checkpoint of amm12592_r128_lp_rounding_boxeph.py with
(a) wider head-cell offsets on row 126 and, decisively,
(b) TAIL-CELL BRANCHING on row 126 (cells k = 0..tailn-1, the ones whose
    leading positions are the small-cap tail levels of the final Bernstein
    decode at d = 152) -- the deterministic clamp cannot steer those levels,
    and the final decode is cap-tight at both ends of position space.
Exact integers only; any hit is fully re-verified before saving."""
import importlib.util, json, os, sys, time
from itertools import product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "lpr", os.path.join(HERE, "amm12592_r128_lp_rounding_boxeph.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

R = 128
d = m.profile(R)
da, db = d[R - 2], d[R - 1]        # 151, 152

def endgame_candidates(sigma, dd, headn, headoff, tailn, tailoff):
    """Row-(R-2) step: forced top cell, forced unit continuation on cell
    d-1; head cells d-2..d-headn with even offsets <= headoff around the
    parity-rounded clamp; TAIL cells k = 0..tailn-1 with even offsets
    <= tailoff (processed last, branching on their leading positions);
    deterministic nearest-parity clamp elsewhere.  Yields (de, sigma_next)."""
    if not sigma or abs(sigma[0]) != 1:
        return
    L = max(len(sigma), dd + 1)
    res0 = list(sigma) + [0] * (L - len(sigma))
    v = sigma[0]
    t = m.qk(dd)
    for s in range(min(dd + 1, L)):
        res0[s] -= v * t[s]
    de0 = [0] * (dd + 1)
    de0[dd] = v
    states = [(res0, de0)]
    # cells k = dd-1 down to tailn, branching on the first headn cellidx
    for cellidx in range(dd - tailn):
        k = dd - 1 - cellidx
        j = dd - k
        cap = comb(dd, k)
        tq = m.qk(k)
        nxt = []
        for res, de in states:
            rj = res[j]
            wstar = -cap if rj < -cap else (cap if rj > cap else rj)
            if cellidx == 0:
                cands = [w for w in (rj - 1, rj + 1)
                         if -cap <= w <= cap and (w - cap) % 2 == 0]
            elif cellidx < headn:
                base = ((wstar,) if (wstar - cap) % 2 == 0
                        else (wstar - 1, wstar + 1))
                cs = set()
                for b in base:
                    for o in range(-headoff, headoff + 1, 2):
                        w = b + o
                        if -cap <= w <= cap and (w - cap) % 2 == 0:
                            cs.add(w)
                cands = sorted(cs)
            else:
                w = wstar
                if (w - cap) % 2:
                    w = w - 1 if w - 1 >= -cap else w + 1
                cands = [w]
            for w in cands:
                if len(cands) == 1:
                    r2, d2 = res, de
                else:
                    r2, d2 = list(res), list(de)
                d2[k] = w
                if w:
                    off = dd - k
                    for s in range(k + 1):
                        r2[off + s] -= w * tq[s]
                nxt.append((r2, d2))
        states = nxt
        if not states:
            return
    # tail cells k = tailn-1 .. 0, each with even offsets around the clamp
    for k in range(tailn - 1, -1, -1):
        j = dd - k
        cap = comb(dd, k)
        tq = m.qk(k)
        nxt = []
        for res, de in states:
            rj = res[j]
            wstar = -cap if rj < -cap else (cap if rj > cap else rj)
            base = ((wstar,) if (wstar - cap) % 2 == 0
                    else (wstar - 1, wstar + 1))
            cs = set()
            for b in base:
                for o in range(-tailoff, tailoff + 1, 2):
                    w = b + o
                    if -cap <= w <= cap and (w - cap) % 2 == 0:
                        cs.add(w)
            for w in sorted(cs, key=lambda x: abs(x - wstar)):
                r2, d2 = list(res), list(de)
                d2[k] = w
                if w:
                    off = dd - k
                    for s in range(k + 1):
                        r2[off + s] -= w * tq[s]
                nxt.append((r2, d2))
        states = nxt
    for res, de in states:
        if res[0] != 0:
            continue
        yield de, m.trim(res[1:])

def main():
    headn = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    headoff = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    tailn = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    tailoff = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    ck = json.load(open(os.path.join(HERE,
                        "amm12592_r128_lp_ckpt_boxeph.json")))
    assert ck["row"] == R - 3
    sts = sorted(ck["states"], key=lambda s: sum(abs(x) for x in s[1]))
    print(f"endgame widener: {len(sts)} states, headn={headn} "
          f"headoff={headoff} tailn={tailn} tailoff={tailoff}", flush=True)
    t0 = time.time()
    tried = 0
    for si, (acc, sig, debts, overs) in enumerate(sts):
        L1 = sum(abs(x) for x in sig)
        for de, ns in endgame_candidates(list(sig), da, headn, headoff,
                                         tailn, tailoff):
            if not ns or abs(ns[0]) != 1:
                continue
            tried += 1
            fin = m.final_decode(list(ns), db)
            if fin is None:
                continue
            sol = [list(b) for b in acc] + [de, fin]
            a_ok = all(m.admissible(sol[i], d[i]) for i in range(R))
            i_ok = m.epoch_identity(R, sol, d)
            fb, fi, fbnd = m.to_f(sol, d)
            s_ok = m.star2_identity(R, fb, d)
            print(f"HIT at state {si} (L1={L1}) after {tried} decodes, "
                  f"{time.time()-t0:.1f}s", flush=True)
            print(f"EXACT VERIFY: admissible={a_ok} identity={i_ok} "
                  f"f-int={fi} f-bnd={fbnd} (**)={s_ok}", flush=True)
            assert a_ok and i_ok and fi and fbnd and s_ok
            wp = os.path.join(HERE, "amm12592_floor_witness_R128_lp.json")
            json.dump({"R": R, "profile": d, "blocks": sol,
                       "verified": True, "H": [1],
                       "source_label": "gamma* floor (LP-clamp + parity-"
                       "rounding beam 400 rows 0..125 + tail-branched "
                       "endgame, Angle 3, independent construction)",
                       "search": {"beam": 400, "nbranch": 2, "maxoff": 2,
                                  "rank": "l1deg", "endgame":
                                  {"headn": headn, "headoff": headoff,
                                   "tailn": tailn, "tailoff": tailoff},
                                  "deterministic": True},
                       "winner_row125_L1": L1},
                      open(wp, "w"))
            print(f"WITNESS WRITTEN: {wp}", flush=True)
            return 0
        if si % 10 == 9:
            print(f"  ... state {si+1}/{len(sts)} done, decodes={tried}, "
                  f"{time.time()-t0:.1f}s", flush=True)
    print(f"no completion: {tried} decodes over {len(sts)} states, "
          f"{time.time()-t0:.1f}s  (search negative, NOT evidence)",
          flush=True)
    return 3

if __name__ == "__main__":
    sys.exit(main())
