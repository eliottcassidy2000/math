"""Decode slim-witness blocks into closed-form atoms (R = 8, 16)."""
import json, os
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

W = load_witnesses()
for R in (8, 16):
    w = W[R]; prof, blocks = w["profile"], w["blocks"]
    print(f"\n===== R={R} blocks =====")
    for i in range(R):
        d = prof[i]
        P = bern_to_poly(blocks[i], d)
        b = ballot(d)
        c = [blocks[i][k]-b[k] for k in range(d+1)]
        print(f"row {i:2d} d={d:2d} cells={blocks[i]}")
        print(f"          poly={P}")
        if any(c): print(f"          c-cells={c}")
