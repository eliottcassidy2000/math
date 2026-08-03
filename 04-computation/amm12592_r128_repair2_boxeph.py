"""AMM 12592 R=128 exact banded 2-row repair (v2, complete within ansatz).

State sigma at row R-3 (sigma_125).  Ansatz: Delta_126 = ride(sgn) + even
corrections 2-adjusted on cells k (a_k, |ride_k + a_k| <= C(151,k), a_k even),
Delta_127 = decode((sigma - Delta_126_poly)/p, 152), required box-admissible.

Cell algebra: adding a_k to row-126 cell k shifts the final decode cells by
    de'_j = de_j - (a_j + 2 a_{j-1} + a_{j-2})
(the (x+q)^2 band from the degree gap d127 - (d126 - 1) = 2).
DFS over k = 0..K with exact interval feasibility per cell; beyond K the
final caps C(152,k) dwarf any band contribution.  The assembled witness is
verified exactly (independent identity + admissibility) before saving.
"""
import json, os, sys, time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r128_floor_solver_boxeph as eng

def decode_cells(poly, d):
    if len(poly) - 1 > d:
        return None
    res = list(poly) + [0] * (d + 1 - len(poly))
    de = [0] * (d + 1)
    for k in range(d, -1, -1):
        want = res[d - k]
        de[k] = want
        if want:
            t = eng.qk(k)
            off = d - k
            for s in range(k + 1):
                res[off + s] -= want * t[s]
    return de if not eng.trim(res) else None

def box_ok(cells, d):
    return len(cells) == d + 1 and all(
        abs(v) <= comb(d, k) and (v - comb(d, k)) % 2 == 0
        for k, v in enumerate(cells))

def poly_of_cells(cells, d):
    poly = [0] * (d + 1)
    for k, v in enumerate(cells):
        if v:
            t = eng.qk(k)
            off = d - k
            for s in range(k + 1):
                poly[off + s] += v * t[s]
    return poly

def repair2(sigma, d1, d2, sgn, K=14, branch=7):
    """Exact banded 2-row repair; returns (de126, de127) or None."""
    ride = eng.final_decode([-sgn, 2 * sgn], d1)
    if ride is None:
        return None
    # base final decode
    y = [-sgn, 2 * sgn]
    n = max(len(sigma), 2)
    rem = [(sigma[t] if t < len(sigma) else 0) -
           (y[t] if t < len(y) else 0) for t in range(n)]
    if not rem or rem[0] != 0:
        return None
    z = eng.trim(rem[1:])
    de = decode_cells(z, d2)
    if de is None:
        return None
    # parity screen (a_k even cannot change parity)
    for k in range(min(K + 6, d2 + 1)):
        if (de[k] - comb(d2, k)) % 2:
            return None
    c2 = [comb(d2, k) for k in range(d2 + 1)]
    c1 = [comb(d1, k) for k in range(d1 + 1)]

    sol = {}

    def cands(lo, hi, want):
        """few even candidates in [lo,hi], preferring near `want`."""
        if lo > hi:
            return []
        w = max(lo, min(hi, want))
        w -= (w - lo) % 2 if (w % 2 != lo % 2) else 0
        # normalize to even values with the parity of lo (lo,hi are even-adjusted)
        out = []
        for v in (w, lo, hi, (lo + hi) // 2, w - 2, w + 2, 0):
            if lo <= v <= hi and v % 2 == 0 and v not in out:
                out.append(v)
        return out[:branch]

    def dfs(k, a_prev1, a_prev2):
        if k > K:
            return True
        band = 2 * a_prev1 + a_prev2
        # a_k must satisfy |de_k - a_k - band| <= c2[k]
        lo1 = de[k] - band - c2[k]
        hi1 = de[k] - band + c2[k]
        # and |ride_k + a_k| <= c1[k]
        rk = ride[k] if k < len(ride) else 0
        lo2, hi2 = -c1[k] - rk, c1[k] - rk
        lo, hi = max(lo1, lo2), min(hi1, hi2)
        # even alignment
        lo += lo % 2
        hi -= hi % 2
        if lo > hi:
            return False
        # prefer a_k that centers the next cell
        want = (de[k + 1] - a_prev1) // 2 if k + 1 <= K else 0
        want -= want % 2
        for a in cands(lo, hi, want):
            sol[k] = a
            if dfs(k + 1, a, a_prev1):
                return True
        sol.pop(k, None)
        return False

    if not dfs(0, 0, 0):
        return None
    de126 = list(ride)
    for k, a in sol.items():
        de126[k] = ride[k] + a
    # exact reconstruction of Delta_127 from the corrected Delta_126
    y2 = poly_of_cells(de126, d1)
    n = max(len(sigma), len(y2))
    rem = [(sigma[t] if t < len(sigma) else 0) -
           (y2[t] if t < len(y2) else 0) for t in range(n)]
    if not rem or rem[0] != 0:
        return None
    z2 = eng.trim(rem[1:])
    de127 = decode_cells(z2, d2)
    if de127 is None or not box_ok(de127, d2):
        return None
    if not box_ok(de126, d1):
        return None
    return de126, de127

def repair_states(states, R, d, out, label, log=print):
    """states: list of (acc, sigma) at row R-3."""
    t0 = time.time()
    for si, (acc, sig) in enumerate(states):
        for sgn in (1, -1):
            r = repair2(sig, d[R - 2], d[R - 1], sgn)
            if r is None:
                continue
            de126, de127 = r
            blocks = list(acc) + [de126, de127]
            a, b = eng.verify_witness(R, blocks, d)
            log(f"  repair2 state {si} sign {sgn}: verify adm={a} id={b}")
            if a and b:
                eng.write_json(out, {
                    "R": R, "profile": d, "blocks": blocks, "verified": True,
                    "H": [1],
                    "source_label": f"gamma* floor ({label} + exact banded "
                                    f"2-row repair)",
                    "search": {"mode": "repair2", "state": si, "sign": sgn}})
                log(f"  WITNESS WRITTEN: {out}")
                return blocks
    log(f"  repair2: no closure over {len(states)} states "
        f"({time.time()-t0:.1f}s)")
    return None

def main():
    R = 128
    d = eng.prof(R)
    out = os.path.join(HERE, "amm12592_floor_witness_R128_direct.json")
    ck_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, "amm12592_r128_beam_checkpoint.json")
    with open(ck_path) as f:
        ck = json.load(f)
    print(f"repair2 on {os.path.basename(ck_path)} (phase {ck.get('phase')}, "
          f"row {ck['row']}, {len(ck['states'])} states)", flush=True)
    assert ck["row"] == R - 3
    states = [(st["acc"], st["sigma"]) for st in ck["states"]]
    r = repair_states(states, R, d, out, f"beam ckpt row {ck['row']}")
    return 0 if r is not None else 1

if __name__ == "__main__":
    sys.exit(main())
