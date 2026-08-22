"""Sparse polynomial-side corrections c_i = Delta_i - (p-q) for slim witnesses.
Print as {deg: coeff} sparse dicts; align windows against d_i and R."""
import json, os
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

W = load_witnesses()
out = {}
for R in (8, 16, 32, 64):
    w = W[R]; prof, blocks = w["profile"], w["blocks"]
    print(f"\n===== R={R}: c_i as sparse polynomials (deg:coeff), backbone (p-q) removed =====")
    rows = []
    for i in range(R):
        d = prof[i]
        P = bern_to_poly(blocks[i], d)
        c = psub(P, [-1, 2])   # subtract p-q
        sp = {j: v for j, v in enumerate(c) if v}
        rows.append({"i": i, "d": d, "c": sp})
        if sp:
            nterms = len(sp); mx = max(abs(v) for v in sp.values())
            lo, hi = min(sp), max(sp)
            tag = f"win=[{lo},{hi}] hi-d={hi-d} terms={nterms} max={mx}"
            print(f"row {i:3d} d={d:3d} {tag}: {sp if nterms<=14 else dict(list(sorted(sp.items()))[:14])}")
    out[R] = rows
json.dump(out, open(os.path.join(HERE, "amm12592_allR_sparse_corrections_boxeph.json"), "w"), indent=0)
print("\nwrote amm12592_allR_sparse_corrections_boxeph.json")
