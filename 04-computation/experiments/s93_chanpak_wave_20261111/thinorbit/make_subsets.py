#!/usr/bin/env python3
"""Lane thinorbit (S93): write sub-alphabets of the 56 (or 264) base-2 block matrices by determinant / rank."""
import sys
import sympy as sp

src = sys.argv[1]
tag = sys.argv[2]
rows = [list(map(int, l.split())) for l in open(src)]
sel = {"uni": [], "inv": [], "sing": [], "rank2": [], "rank1": [], "nonuni_inv": []}
for r in rows:
    M = sp.Matrix(3, 3, r)
    d = int(M.det())
    rk = M.rank()
    if abs(d) == 1:
        sel["uni"].append(r)
    if d != 0:
        sel["inv"].append(r)
        if abs(d) != 1:
            sel["nonuni_inv"].append(r)
    else:
        sel["sing"].append(r)
        if rk == 2:
            sel["rank2"].append(r)
        else:
            sel["rank1"].append(r)
for k, v in sel.items():
    with open(f"sub_{tag}_{k}.txt", "w") as f:
        for r in v:
            f.write(" ".join(map(str, r)) + "\n")
    print(k, len(v))
