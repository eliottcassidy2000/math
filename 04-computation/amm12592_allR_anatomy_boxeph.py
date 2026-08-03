"""Witness anatomy in ballot coordinates for R = 8,16,32,64,128."""
import json, os
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

W = load_witnesses()
for R in sorted(W):
    w = W[R]; prof, blocks = w["profile"], w["blocks"]
    assert epoch_sum(R, blocks, prof) == qpow(R-1), f"R={R} identity FAIL"
    assert all(admissible(blocks[i], prof[i]) for i in range(R)), f"R={R} adm FAIL"
    sig = residuals(R, blocks, prof)
    print(f"\n================ R = {R} ================  (identity + admissibility re-verified)")
    # last row
    last = blocks[-1]; d = prof[-1]
    is_minus_box = all(last[k] == -comb(d,k) for k in range(d+1))
    is_plus_box  = all(last[k] ==  comb(d,k) for k in range(d+1))
    print(f"row R-1: -fullbox={is_minus_box} +fullbox={is_plus_box}")
    # attractor entry: sigma_i == E_m = -1 + x + ... + x^m
    def isEm(P):
        if not P or P[0] != -1: return None
        if all(v == 1 for v in P[1:]): return len(P)-1
        return None
    entry = None
    for i in range(R):
        m = isEm(sig[i])
        if m is not None and entry is None:
            entry = (i, m)
    print(f"attractor entry: sigma_i = E_m first at i={entry[0] if entry else None} m={entry[1] if entry else None}")
    # ballot deviation per row
    dev_rows = []
    for i in range(R):
        d = prof[i]
        c = [blocks[i][k] - b for k, b in enumerate(ballot(d))]
        if any(c):
            nz = [k for k, v in enumerate(c) if v]
            dev_rows.append((i, min(nz), max(nz), max(abs(v) for v in c)))
    print(f"rows with c_i != 0: {len(dev_rows)} of {R}")
    for i, k0, k1, mx in dev_rows[:40]:
        print(f"  row {i:3d} d={prof[i]:3d}: cell-support k=[{k0},{k1}] max|c|={mx}")
    if len(dev_rows) > 40: print(f"  ... ({len(dev_rows)-40} more)")
    # C cross-check for slim/backbone witnesses: Delta_0 + sum_{1..R-2} x^i c_i == closed_C(R)?
    Cw = bern_to_poly(blocks[0], prof[0])
    for i in range(1, R-1):
        ci = psub(bern_to_poly(blocks[i], prof[i]), [-1, 2])
        Cw = padd(Cw, pshift(ci, i))
    match = (Cw == closed_C(R)) if is_minus_box else None
    print(f"C(witness) == closed_C(R): {match}" + ("" if is_minus_box else "  [row R-1 not -fullbox; C-form needs Delta_(R-1)=-1]"))
