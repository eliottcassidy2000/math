"""AMM 12592 R=128 exact 2-row BOX-PUSH repair (fast, int-only).

For a row-125 state sigma with sigma = sgn*E_1 + J (J even):
    Delta_126 := ride  (decode of sgn*(2x-1) at d126),
    Delta_127 := -sgn*fullbox + W,   W = decode(J/p, d127)
works iff p | J and each cell -sgn*C(152,k) + W_k lies in the box with parity.
Parity is automatic (J even => W even; box cell parity = C(152,k)).
The full box sits AT the cap, so mid-k cells absorb ~2*C(152,k) junk; only
small-k cells (caps 1, 2*152, 2*11476, ...) truly constrain.

Fallback within the same ansatz: move a small even correction 2u into the
ride row at low-order cells (k = 0..KREP) to fix small-k violations of W:
Delta_126 = ride + sum_k 2 a_k B_{151,k} shifts W by -decode(2u/p ...) --
implemented as a small search over a_k in the row-126 margins for the cells
k = 1..KREP, adjusting W accordingly (exact, brute force over tiny windows).

Every witness is verified exactly before saving.  Usage:
    python3 amm12592_r128_boxpush_repair_boxeph.py [checkpoint.json]
"""
import json, os, sys, time
from math import comb
from itertools import product as iproduct

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r128_floor_solver_boxeph as eng

def decode_cells(poly, d):
    """Exact Bernstein-d cell decode of poly (must have deg <= d), no
    admissibility test; returns cells or None if deg too high / leftover."""
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

def box_admissible(cells, d):
    return all(abs(v) <= comb(d, k) and (v - comb(d, k)) % 2 == 0
               for k, v in enumerate(cells)) and len(cells) == d + 1

def poly_of_cells(cells, d):
    poly = [0] * (d + 1)
    for k, v in enumerate(cells):
        if v:
            t = eng.qk(k)
            off = d - k
            for s in range(k + 1):
                poly[off + s] += v * t[s]
    return poly

def boxpush(sigma, d126, d127, sgn, KREP=4, WIN=40):
    """Try Delta_126 = ride + small low-cell correction, Delta_127 =
    -sgn*box + decode((sigma - Delta_126_poly)/p).  Exact; returns
    (de126, de127) or None."""
    ride = eng.final_decode([-sgn, 2 * sgn], d126)
    if ride is None:
        return None
    box = [-sgn * comb(d127, k) for k in range(d127 + 1)]
    base126 = [-sgn, 2 * sgn]
    # J = sigma - ride_poly (as polynomials); Delta_127 = (J)/p + box-part:
    # sigma = Delta126 + p*Delta127; with Delta126 = ride: p*Delta127 =
    # sigma - (-sgn + 2sgn x); Delta127 = that / p; then subtract box base
    # implicitly by requiring box-admissibility of box + W where
    # W = Delta127_cells - box... simpler: decode Delta127 directly.
    def attempt(corr):
        """corr: dict k -> even correction to ride cell k (row 126)."""
        de126 = list(ride)
        for k, a in corr.items():
            de126[k] = ride[k] + a
            if abs(de126[k]) > comb(d126, k):
                return None
        y = poly_of_cells(de126, d126) if corr else list(base126) + [0] * (d126 - 1)
        n = max(len(sigma), len(y))
        rem = [(sigma[t] if t < len(sigma) else 0) -
               (y[t] if t < len(y) else 0) for t in range(n)]
        if rem and rem[0] != 0:
            return None
        z = eng.trim(rem[1:])
        de127 = decode_cells(z, d127)
        if de127 is None:
            return None
        if not box_admissible(de127, d127):
            return None
        return de126, de127
    r = attempt({})
    if r is not None:
        return r
    # small corrections on low-order row-126 cells (affect top positions):
    # enumerate even windows on cells 1..KREP one and two at a time
    cand = []
    for k in range(1, KREP + 1):
        m = comb(d126, k) - abs(ride[k])
        w = min(WIN, m)
        cand.append([2 * t for t in range(-w // 2, w // 2 + 1)])
    for k in range(1, KREP + 1):
        for a in cand[k - 1]:
            if a == 0:
                continue
            r = attempt({k: a})
            if r is not None:
                return r
    for k1 in range(1, KREP + 1):
        for k2 in range(k1 + 1, KREP + 1):
            for a1 in cand[k1 - 1]:
                if a1 == 0:
                    continue
                for a2 in cand[k2 - 1]:
                    if a2 == 0:
                        continue
                    r = attempt({k1: a1, k2: a2})
                    if r is not None:
                        return r
    return None

def main():
    R = 128
    d = eng.prof(R)
    out = os.path.join(HERE, "amm12592_floor_witness_R128_direct.json")
    ck_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, "amm12592_r128_beam_checkpoint.json")
    with open(ck_path) as f:
        ck = json.load(f)
    row = ck["row"]
    print(f"box-push repair on {os.path.basename(ck_path)} "
          f"(phase {ck.get('phase')}, row {row}, {len(ck['states'])} states)",
          flush=True)
    assert row == R - 3
    t0 = time.time()
    for si, st in enumerate(ck["states"]):
        sig = st["sigma"]
        for sgn in (1, -1):
            r = boxpush(sig, d[R - 2], d[R - 1], sgn)
            if r is None:
                continue
            de126, de127 = r
            blocks = list(st["acc"]) + [de126, de127]
            a, b = eng.verify_witness(R, blocks, d)
            print(f"  state {si} sign {sgn}: repair found; verify adm={a} "
                  f"id={b}", flush=True)
            if a and b:
                eng.write_json(out, {
                    "R": R, "profile": d, "blocks": blocks, "verified": True,
                    "H": [1],
                    "source_label": "gamma* floor (beam1000 + exact box-push "
                                    "2-row tail repair)",
                    "search": {"mode": "box-push", "state": si, "sign": sgn,
                               "ckpt": os.path.basename(ck_path)}})
                print(f"  WITNESS WRITTEN: {out}", flush=True)
                print(f"CLOSED ({time.time()-t0:.1f}s)")
                return 0
        if si % 10 == 9:
            print(f"  ..{si+1}/{len(ck['states'])} ({time.time()-t0:.1f}s)",
                  flush=True)
    print(f"no box-push repair on {len(ck['states'])} states "
          f"({time.time()-t0:.1f}s)")
    return 1

if __name__ == "__main__":
    sys.exit(main())
