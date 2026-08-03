"""AMM 12592 R=128 exact 2-row TAIL REPAIR.

Input: a beam state sigma_125 = sgn*E_1 + J where the junk J is supported on
high positions (the diagnosed beam1000 population: head = attractor, interior
zero, 3-term tail junk ~6e4 at positions 148..150, which overflows the small
top-cell caps of a single-row decode).

Closure needs Delta_126 (admissible, d=151) and Delta_127 (admissible, d=152)
with  sigma_125 = Delta_126 + p*Delta_127.
Write Delta_126 = ride(2x-1 decode at 151)*sgn + D1, Delta_127 = full-box
(-sgn) + D2; then the requirement is exactly

    D1 + p*D2 = J,   D1 = sum_k 2 a_k B_{151,k},  D2 = sum_k 2 b_k B_{152,k},

with box margins |ride_k + 2 a_k| <= C(151,k) (parity preserved since the
corrections are even) and |(-C(152,k)) + 2 b_k| <= C(152,k) i.e.
0 <= b_k <= C(152,k).

With J supported on positions >= 148 and deg constraints (positions 152,153
must vanish), only small-k cells participate: a_k (k=0..KA), b_k (k=2..KB).
The linear system is solved exactly over Z (sympy), then the affine kernel is
enumerated over a bounded window, checking all box margins.  Every produced
witness is verified independently (exact identity + admissibility).
"""
import json, os, sys, time
from math import comb
from itertools import product as iproduct

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r128_floor_solver_boxeph as eng
import amm12592_r128_hook_solver_boxeph as hk

import sympy

def cell_poly(d, k):
    """dense coeffs of B_{d,k} = x^{d-k}(1-x)^k (length d+1)."""
    v = [0] * (d + 1)
    t = eng.qk(k)
    for s, c in enumerate(t):
        v[d - k + s] = c
    return v

def tail_repair(sigma, d126, d127, sgn, KA=6, KB=8, win=3):
    """Try to write sigma = Delta_126 + p*Delta_127 with the ride/full-box
    ansatz + even tail corrections.  Returns (de126, de127) or None."""
    # base blocks
    ride = eng.final_decode([-sgn, 2 * sgn], d126)
    if ride is None:
        return None
    box = [-sgn * comb(d127, k) for k in range(d127 + 1)]  # full box -sgn
    # J = sigma - (ride_poly + p*box_poly); ride_poly = -sgn+2sgn x, p*box = -sgn x
    # base sum = -sgn + 2sgn x - sgn x = -sgn + sgn x = sgn*E_1. J = sigma - sgn*E_1.
    n = max(len(sigma), 2)
    J = list(sigma) + [0] * (n - len(sigma))
    J[0] -= -sgn
    J[1] -= sgn
    J = eng.trim(J)
    if not J:
        return ride, box  # pure attractor
    if any(c % 2 for c in J):
        return None
    # unknown cells: a_k at d126 (k=0..KA), b_k at d127 (k=0..KB), b as poly
    # shifted by 1 (multiplied by p).  Build equations on positions
    # pos0 .. d127+1 restricted to where anything is nonzero.
    N = d127 + 2  # positions 0..d127+1 for p*B_{d127,k}
    Avars = sympy.symbols(f"a0:{KA+1}", integer=True)
    Bvars = sympy.symbols(f"b0:{KB+1}", integer=True)
    coefvec = [sympy.Integer(0)] * N
    for k in range(KA + 1):
        t = eng.qk(k)
        for s, c in enumerate(t):
            coefvec[d126 - k + s] += 2 * c * Avars[k]
    for k in range(KB + 1):
        t = eng.qk(k)
        for s, c in enumerate(t):
            coefvec[1 + d127 - k + s] += 2 * c * Bvars[k]
    Jv = list(J) + [0] * (N - len(J))
    eqs = []
    for pos in range(N):
        e = coefvec[pos] - Jv[pos]
        if e != 0:
            eqs.append(e)
    sol = sympy.solve(eqs, list(Avars) + list(Bvars), dict=True)
    if not sol:
        return None
    sol = sol[0]
    # free variables = those not in sol
    allv = list(Avars) + list(Bvars)
    free = [v for v in allv if v not in sol]
    def caps_ok(vals):
        de126 = list(ride)
        de127 = list(box)
        for k in range(KA + 1):
            de126[k] = ride[k] + 2 * vals[Avars[k]]
            if abs(de126[k]) > comb(d126, k):
                return None
        for k in range(KB + 1):
            de127[k] = box[k] + 2 * vals[Bvars[k]]
            if abs(de127[k]) > comb(d127, k):
                return None
        return de126, de127
    rng = range(-win, win + 1)
    for combo in iproduct(rng, repeat=len(free)):
        subs = dict(zip(free, [sympy.Integer(c) for c in combo]))
        vals = {}
        ok = True
        for v in allv:
            val = sympy.Integer(subs[v]) if v in subs else sol[v].subs(subs)
            if not val.is_integer:
                ok = False
                break
            vals[v] = int(val)
        if not ok:
            continue
        r = caps_ok(vals)
        if r is not None:
            return r
    return None

def try_state(R, d, acc, sigma, out, label):
    for sgn in (1, -1):
        r = tail_repair(sigma, d[R - 2], d[R - 1], sgn)
        if r is None:
            continue
        de126, de127 = r
        blocks = list(acc) + [de126, de127]
        a, b = eng.verify_witness(R, blocks, d)
        print(f"    repair candidate sign {sgn}: verify adm={a} id={b}")
        if a and b:
            eng.write_json(out, {
                "R": R, "profile": d, "blocks": blocks, "verified": True,
                "H": [1],
                "source_label": f"gamma* floor ({label} + exact 2-row tail repair)",
                "search": {"mode": "tail-repair", "sign": sgn}})
            print(f"    WITNESS WRITTEN: {out}")
            return True
    return False

def main():
    R = 128
    d = eng.prof(R)
    out = os.path.join(HERE, "amm12592_floor_witness_R128_direct.json")
    ck_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, "amm12592_r128_beam_checkpoint.json")
    with open(ck_path) as f:
        ck = json.load(f)
    row = ck["row"]
    print(f"tail repair on checkpoint {os.path.basename(ck_path)} "
          f"(phase {ck.get('phase')}, row {row}, {len(ck['states'])} states)")
    assert row == R - 3, f"need row R-3={R-3} states (sigma_125), got {row}"
    t0 = time.time()
    for si, st in enumerate(ck["states"]):
        sig = st["sigma"]
        # quick screen: head must be +-E_1-like and interior small
        if try_state(R, d, st["acc"], sig, out, f"beam1000 ckpt state {si}"):
            print(f"CLOSED from state {si} ({time.time()-t0:.1f}s)")
            return 0
        if si % 10 == 9:
            print(f"  ..{si+1}/{len(ck['states'])} no repair yet "
                  f"({time.time()-t0:.1f}s)", flush=True)
    print(f"no tail repair on {len(ck['states'])} states "
          f"({time.time()-t0:.1f}s)")
    return 1

if __name__ == "__main__":
    sys.exit(main())
