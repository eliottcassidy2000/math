#!/usr/bin/env python3
"""Wavefront probe (angle B1 die anatomy): run rule A (= rule B policy none)
at (R, D0) and dump, every `step` rows, the bit-length profile of the lowest
`W` coefficients of sigma_i plus the E-form deviation depth: the smallest j
such that sigma_i[j] != E_{R-2-i}[j].  Exact ints only."""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amm12592_ruleB_ballot_scheduled_boxeph import (
    crow, qpow, Em, ptrim, psub, floor_gamma_star, bern_to_poly, clamp_cellwise, RES)

def probe(R, D0, W=48, step=4):
    prof = [floor_gamma_star(R + i) + D0 for i in range(R)]
    sig = qpow(R - 1)
    rows = []
    t0 = time.time()
    for i in range(R):
        d = prof[i]
        m = R - 2 - i
        ideal = list(sig)
        if m >= 0:
            need = m + 2
            if len(ideal) < need:
                ideal += [0] * (need - len(ideal))
            for j, c in enumerate(Em(m)):
                ideal[j + 1] -= c
        delta, junk = clamp_cellwise(ideal[:d + 1], d, 'tozero')
        Dp = bern_to_poly(delta, d)
        t = psub(sig, Dp)
        if t and t[0] != 0:
            rows.append(dict(i=i, die=True, const_bits=abs(t[0]).bit_length(),
                             head=[abs(c).bit_length() for c in sig[:W]]))
            print(f"R={R} D0={D0}: DIE row {i} const 2^{abs(t[0]).bit_length()} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            return rows, None
        sig = t[1:] if t else []
        if i % step == 0 or i >= R - 3:
            E = Em(m - 1) if m - 1 >= -1 else []
            EE = list(E) + [0] * (max(0, len(sig) - len(E)))
            dev = next((j for j in range(len(sig)) if sig[j] != EE[j]), -1)
            rows.append(dict(i=i, head=[abs(c).bit_length() for c in sig[:W]],
                             efdev=dev, junkbits=junk.bit_length(),
                             deg=len(sig) - 1))
    print(f"R={R} D0={D0}: survived all rows, residual L1="
          f"{sum(abs(c) for c in sig)} ({time.time()-t0:.0f}s)", flush=True)
    return rows, sig

if __name__ == "__main__":
    R = int(sys.argv[1]); D0 = int(sys.argv[2])
    rows, sig = probe(R, D0)
    out = os.path.join(RES, f"amm12592_ruleB_wavefront_R{R}_D{D0}_boxeph.json")
    with open(out, "w") as f:
        json.dump({"R": R, "D0": D0, "rows": rows}, f)
    print("wrote", out)
    # console summary: wavefront = first coefficient with bits > 1
    for r in rows:
        if r.get("die") or (r["i"] % 16 == 0):
            hb = r["head"]
            wf = next((j for j, b in enumerate(hb) if b > 1), -1)
            amp = max(hb) if hb else 0
            print(f"  i={r['i']:4d} wf_pos={wf:3d} head_amp_bits={amp:4d} "
                  f"efdev={r.get('efdev', 'DIE')} junkbits={r.get('junkbits', -1)}")
