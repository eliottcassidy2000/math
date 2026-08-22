"""Dump residual trajectories sigma_i for slim witnesses; hunt closed forms."""
import json, os
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

W = load_witnesses()
for R in (8, 16, 32):
    w = W[R]; prof, blocks = w["profile"], w["blocks"]
    sig = residuals(R, blocks, prof)
    print(f"\n===== R={R} residual trajectory =====")
    print(f"sigma_-1 = q^{R-1}")
    for i in range(R):
        P = sig[i]
        L1 = sum(abs(v) for v in P)
        if R <= 16 or L1 <= 60 or i < 6:
            print(f"i={i:2d} deg={len(P)-1 if P else -1:3d} L1={L1:6d}  {P if len(P)<=40 else '(long)'}")
        else:
            print(f"i={i:2d} deg={len(P)-1 if P else -1:3d} L1={L1}")
